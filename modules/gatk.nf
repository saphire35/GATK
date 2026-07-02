process QCClean {

    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val cpus_per_task_ch
    val galore_quality 
    val galore_length

    output:
    tuple val(sample_id), path("*_val_*.fq"), emit: clean_reads_ch
    path "QC_reports/*"
    tuple val(sample_id), path("*.log"), emit: QCClean_logs

    script:
    """
    mkdir -p fastqc_raw fastqc_clean trimmed

    read1=${reads[0]}
    read2=${reads[1]}

    echo "Processing $sample_id"
    echo "R1: \$read1"
    echo "R2: \$read2"
#-q 30 -e 0.1 n 3 -O 6 m 3
    # FastQC raw + Trim Galore
    fastqc \$read1 \$read2 -o fastqc_raw -t ${cpus_per_task_ch}  >> ${sample_id}_fastqc_raw.log 2>&1
    trim_galore \
        -q ${galore_quality} --length ${galore_length} \
        --illumina --paired \$read1 \$read2 -o trimmed  \
        >> ${sample_id}_fastqc_raw.log 2>&1

    # FastQC trimmed
    fastqc trimmed/*_val_1.fq trimmed/*_val_2.fq -o fastqc_clean -t ${cpus_per_task_ch} >> ${sample_id}_fastqc_raw.log 2>&1
    mv trimmed/*_val_*.fq .
    
    echo "1. QCClean" >> ${sample_id}_gatk_qc.log
    echo "============================================================================== " >> ${sample_id}_gatk_qc.log
    grep -E "Total reads processed|Reads written|Number of sequence pairs removed"   ${sample_id}_fastqc_raw.log >> ${sample_id}_gatk_qc.log
    grep -A2 "Processing " ${sample_id}_fastqc_raw.log  | grep -v -A1  "Processing reads " >> ${sample_id}_gatk_qc.log


    mkdir QC_reports
    cp trimmed/*report.txt fastqc_clean/*.html fastqc_raw/*.html QC_reports/
    """
}


process MergeFastq {

    input:
    path(reads)

    output:
    tuple val("merged"), path("merged_*.fq"), emit: merged_reads_ch

    script:
    """
    cat *_1_val_1.fq > merged_1.fq
    cat *_2_val_2.fq > merged_2.fq
    """
}

process IndexReference {

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: reference_with_index
    path "*.dict", emit: reference_dict
    tuple path(reference), path("${reference}.*"), path("*.dict"), emit: reference_with_index2
    path "${reference}.fai", emit: reference_fai
    path "*log", emit: IndexReference_logs

    script:
    """
    bwa index $reference
    samtools faidx $reference
    java -jar  /usr/local/bin/picard.jar CreateSequenceDictionary \
        R=$reference \
        O=${reference.baseName}.dict

    echo "2. IndexReference " >> gatk_reference.log
    echo "============================================================================== \n" >> gatk_reference.log
    echo "number of bases: \$(awk '{sum += \$2} END {print sum}' *.fai)" >> gatk_reference.log
    echo "number of contigs: \$(wc -l < *.fai) \n" >> gatk_reference.log
    
    """
}

process ReadAlign {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(index_files)
    val cpus_per_task_ch
    val F
    val f
    val q

    output:
    tuple val(sample_id), path ("${sample_id}.sorted.raw.bam*"), path ("${sample_id}.mappingstats.txt"), emit: aligned_bams_ch
    tuple path ("*sorted.raw.bam*"), path ("*sorted.raw.bam.bai"), emit: reference_index_ch
    tuple val(sample_id), path("*.log"), emit: ReadAlign_logs 

    script:
    """
    read1=${reads[0]}
    read2=${reads[1]}   

    #temp alignment without filter
    bwa mem -t "$cpus_per_task_ch" "$reference" "\$read1" "\$read2" > tmp.sam

    total=\$(samtools view -c tmp.sam)
    mapped=\$(samtools view -c -F 4 tmp.sam)
    mapq=\$(samtools view -c -q "$q" tmp.sam)
    minus_F=\$(samtools view -c -q "$q" -F "$F" tmp.sam)
    minus_f=\$(samtools view -c -q "$q" -F "$F" -f "$f" tmp.sam)

    #filter and sort
    samtools view -h -q "$q" -F "$F" -f "$f" -b tmp.sam \
    | samtools sort -@ "$cpus_per_task_ch" -o "${sample_id}.sorted.raw.bam"
    samtools index "${sample_id}.sorted.raw.bam"

    kept=\$(samtools view -c "${sample_id}.sorted.raw.bam")
    removed=\$((total - kept))

    samtools flagstat "${sample_id}.sorted.raw.bam" > "${sample_id}.mappingstats.txt"
    samtools quickcheck "${sample_id}.sorted.raw.bam"

    echo "3. ReadAlign " >> ${sample_id}_gatk_align.log
    echo "==============================================================================" >> ${sample_id}_gatk_align.log
    echo "Filtering summary:" >> ${sample_id}_gatk_align.log
    echo "  Total alignments before filtering  for ${sample_id}  : \$total" >> ${sample_id}_gatk_align.log
    echo "  Mapped alignments (-F 4 only)       : \$mapped" >> ${sample_id}_gatk_align.log
    echo "  After MAPQ filter (-q $q)            : \$mapq" >> ${sample_id}_gatk_align.log
    echo "  After excluded flags (-F $F)         : \$minus_F" >> ${sample_id}_gatk_align.log
    echo "  Final kept alignments (-f $f)        : \$kept" >> ${sample_id}_gatk_align.log
    echo "  Total removed alignments  for ${sample_id}  : \$removed" >> ${sample_id}_gatk_align.log
    echo "" >> ${sample_id}_gatk_align.log

    cat *.mappingstats.txt >> ${sample_id}_gatk_align.log
    rm tmp.sam
    """
}

process AddReadGroups {
    tag "$sample_id"

    input:
    tuple val(sample_id), path (bam_files), path (mapstat_file)
    path metadata


    output:
    //tuple val(sample_id), path ("*sorted.rg*"), emit: sorted_bams_ch
    tuple val(sample_id), path("*.sorted.rg.bam"), emit: sorted_bams_ch
    tuple val(sample_id), path("*log"), emit: AddReadGroups_logs

    script:
       def rgsm_value = params.merged_vcf ? sample_id : '$simpleID'
    def rglb_value = params.merged_vcf ? sample_id : '$simpleID'
    """
    #retrieve metadata for the sample
    #order: simple_ID sample_ID read1 read2 instrument flowcell lane barcode sex run_num seq_num
    line=\$(awk -v id="${sample_id}" 'BEGIN{FS="\t"} \$2==id {print; exit}' $metadata)

    if [ -z "\$line" ]; then
        line=\$(awk -v id="${sample_id}" 'BEGIN{FS="\t"} \$1==id {print; exit}' $metadata)
    fi

    if [ -z "\$line" ]; then
        echo "[ERROR] No metadata line found for sample_id=${sample_id}" >&2
        echo "[ERROR] Metadata columns expected: simple_ID sample_ID read1 read2 instrument flowcell lane barcode sex run_num seq_num" >&2
        exit 1
    fi

    simpleID=\$(echo "\$line" | cut -f1)
    instrument=\$(echo "\$line" | cut -f5 )
    flowcell=\$(echo "\$line" | cut -f6 )
    lane=\$(echo "\$line" | cut -f7 )
    barcode=\$(echo "\$line" | cut -f8 )
    seqnum=\$(echo "\$line" | cut -f11 )
    bam_file=${bam_files[0]}
    echo "launching AddOrReplaceReadGroups for $sample_id with bam file \$bam_file"
    java -Xmx10g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups \
        I=\$bam_file \
        O=${sample_id}.sorted.rg.bam \
        RGSM=${rgsm_value} \
        RGLB=${rglb_value}.\${seqnum} \
        RGID=\${flowcell}.\${lane} \
        RGPU=\${flowcell}\${lane}.\${barcode} \
        RGPL=\$instrument >> ${sample_id}_add.log 2>&1

    # Index
    java -Xmx10g -jar /usr/local/bin/picard.jar BuildBamIndex I=${sample_id}.sorted.rg.bam >> ${sample_id}_add.log 2>&1


    echo "4. AddReadGroups" >> ${sample_id}_gatk_addrg.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_addrg.log
    grep -E "AddOrReplaceReadGroups.*Created"  ${sample_id}_add.log >> ${sample_id}_gatk_addrg.log
    echo "\n" >> ${sample_id}_gatk_addrg.log

    """
}

process Dedup {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file)
    path metadata

    output:
    tuple val(sample_id), path("${sample_id}.sorted.dedup.bam"), path("${sample_id}.sorted.dedup.bai"), emit: dedup_bams_ch
    tuple val(sample_id), path("*.log"), emit: Dedup_logs

    script:
    """
    echo "launching MarkDuplicates for $sample_id with bam file $bam_file"

    java -Xmx10g -jar /usr/local/bin/picard.jar MarkDuplicates \
        I=$bam_file \
        O=${sample_id}.sorted.dedup.bam \
        M=${sample_id}.dedup.metrics.txt \
        REMOVE_DUPLICATES=true

    java -Xmx10g -jar /usr/local/bin/picard.jar BuildBamIndex \
        I=${sample_id}.sorted.dedup.bam


    echo "5. Dedup " >> ${sample_id}_gatk_dedup.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_dedup.log
    grep -A100 "## METRICS CLASS" ${sample_id}.dedup.metrics.txt >> ${sample_id}_gatk_dedup.log
    echo "\n" >> ${sample_id}_gatk_dedup.log
    """
}

process MergeBams {

    input:
    path bam_files
    path metadata

    output:
    path "merged.bam", emit: merged_bam

    script:
    //build input string
    def inputs = bam_files.collect { bam -> "I=${bam}" }.join(' ')
    """
    java -Xmx20g -jar /usr/local/bin/picard.jar MergeSamFiles \
        ${inputs} \
        O=merged.bam \
        TMP_DIR=/gpfs/ts0/scratch/mv323/tmp \
        VALIDATION_STRINGENCY=LENIENT
    """
}



process  HaplotypeCaller{
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*.vcf*"


    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    tuple path(reference), path(reference_idx)
    path reference_dict
    val cpus_per_task_ch

    output:
    tuple val(sample_id),
          path("${sample_id}.g.vcf.gz"),
          path("${sample_id}.g.vcf.gz.tbi"),
          emit: gatk_gvcf_ch


    script:
    def total_cpus = cpus_per_task_ch as int
    def gatk_cpus = Math.max(1, Math.min(total_cpus.intdiv(2), 8))

    """
    gatk --java-options "-Xmx10g" HaplotypeCaller \
        -R "$reference" \
        -I "$bam_file" \
        -O "${sample_id}.g.vcf.gz" \
        -ERC GVCF \
        --sample-ploidy 1 \
        --native-pair-hmm-threads ${gatk_cpus}


    mkdir -p "${params.outdir}"
    cp "${sample_id}.g.vcf.gz" "${sample_id}.g.vcf.gz.tbi" "${params.outdir}/"
    """

}

process GenotypeInterval {
    tag "${sample_id}_${interval}"

    input:
    tuple   val(interval),
            val(sample_id),
            path(gvcfs),
            path(tbis),
            path(reference_files)

    output:
    tuple val(sample_id),
          path("${sample_id}_${interval}.vcf.gz"),
          path("${sample_id}_${interval}.vcf.gz.tbi")


    script:
    def reference = reference_files.find { file -> file.name.endsWith(".fasta") }
 def gvcf_list = gvcfs instanceof List ? gvcfs : [gvcfs]
    def variants = gvcf_list.collect { gvcf -> "--variant ${gvcf}" }.join(' ')

    """
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR

    gatk --java-options "-Djava.io.tmpdir=\$TMPDIR" GenomicsDBImport \\
        ${variants} \\
        --genomicsdb-workspace-path INTERVAL_${sample_id}_${interval}_db \\
        --intervals ${interval}

    gatk --java-options "-Djava.io.tmpdir=\$TMPDIR" GenotypeGVCFs \\
        -R ${reference} \\
        -V gendb://INTERVAL_${sample_id}_${interval}_db \\
        -O ${sample_id}_${interval}.vcf.gz

    tabix -f -p vcf ${sample_id}_${interval}.vcf.gz

    rm -rf INTERVAL_${sample_id}_${interval}_db
    """
}


process ConcatVCFs {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(vcfs),
          path(tbis)

    output:
    tuple val(sample_id),
          path("${sample_id}_genotyped.vcf.gz"),
          path("${sample_id}_genotyped.vcf.gz.tbi")

    script:
    def sorted_vcfs = vcfs.sort { it.name }
    def inputs = sorted_vcfs.collect { "-I ${it}" }.join(' ')

    """
    gatk --java-options "-Xmx10g" GatherVcfs \
        ${inputs} \
        -O ${sample_id}_genotyped.vcf.gz


    tabix -f -p vcf ${sample_id}_genotyped.vcf.gz
    """
}





process HaplotypeCallerFreeBayes {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*vcf*"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    tuple path(reference), path(reference_idx)
    path reference_dict
    val cpus_per_task_ch

    output:
    path("${sample_id}.freebayes.vcf"), emit: freebayes_vcf_ch
    tuple val(sample_id), path("*log"), emit: HaplotypeCallerFreeBayes_logs

    script:
    """
    #index file
    FAI_FILE=\$(ls *.fai)
    
    #generate regions(contigs) of 100.000 to speed things up by //
    fasta_generate_regions.py "\$FAI_FILE" 100000 > regions.txt
    #standard values
    freebayes-parallel regions.txt ${cpus_per_task_ch} \
        -f "$reference" \
        -p 1 \
        --min-alternate-count 3 \
        --min-alternate-fraction 0.2 \
        "$bam_file" \
        > "${sample_id}.freebayes.vcf"
#removed --pooled-discrete \
    echo "7. freebayes variant number for $sample_id:" >> ${sample_id}_gatk_hcFB.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_hcFB.log
    grep -v "^#" "${sample_id}.freebayes.vcf" | wc -l >> ${sample_id}_gatk_hcFB.log
    mkdir -p "${params.outdir}"
    cp "${sample_id}.freebayes.vcf" "${params.outdir}/${sample_id}.freebayes.vcf"
    """
}



process CompareCallers {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*vcf*"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(gatk_gvcf)
    path freebayes_vcf

    output:
    tuple val(sample_id), path("${sample_id}.gatk.raw.vcf")
    tuple val(sample_id), path("${sample_id}.freebayes.sorted.vcf.gz"), emit: freebayes_vcf
    tuple val(sample_id), path("${sample_id}.gatk.sorted.vcf.gz"), emit: gatk_vcf
    tuple val(sample_id), path("${sample_id}.intersect.position.vcf"), emit: intersect_position_vcf
    tuple val(sample_id), path("${sample_id}.intersect.allele.strict.vcf"), emit: intersect_allele_strict_vcf
    tuple val(sample_id), path("*.log"), emit: CompareCallers_logs

    script:
    """
    bcftools view -v snps,indels "$gatk_gvcf" > ${sample_id}.gatk.raw.vcf

    bcftools sort ${sample_id}.gatk.raw.vcf -o ${sample_id}.gatk.sorted.vcf
    bcftools sort "$freebayes_vcf" -o ${sample_id}.freebayes.sorted.vcf

    bgzip -f ${sample_id}.gatk.sorted.vcf
    bgzip -f ${sample_id}.freebayes.sorted.vcf

    tabix -f -p vcf ${sample_id}.gatk.sorted.vcf.gz
    tabix -f -p vcf ${sample_id}.freebayes.sorted.vcf.gz

    bedtools intersect -header \
        -a ${sample_id}.gatk.sorted.vcf.gz \
        -b ${sample_id}.freebayes.sorted.vcf.gz \
        > ${sample_id}.intersect.position.vcf

    bcftools isec -n=2 -w1 \
        ${sample_id}.gatk.sorted.vcf.gz \
        ${sample_id}.freebayes.sorted.vcf.gz \
        -o ${sample_id}.intersect.allele.strict.vcf


    echo "9. CompareCallers (${sample_id}): number of intersected SNPs:" >> ${sample_id}_gatk_comparison.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_comparison.log
    echo "Position-based intersection: \n" >> ${sample_id}_gatk_comparison.log
    grep -v "^#" ${sample_id}.intersect.allele.strict.vcf | wc -l >> ${sample_id}_gatk_comparison.log
    echo "Allele-based strict intersection: \n" >> ${sample_id}_gatk_comparison.log
    grep -v "^#" ${sample_id}.intersect.position.vcf | wc -l >> ${sample_id}_gatk_comparison.log
    """
}

process IntersectRunVCFs {

    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*.vcf*"

    input:
    tuple val(sample_id), path(vcf_files)

    output:
    tuple val(sample_id), path("${sample_id}.common_runs.allele.vcf"), emit: common_runs_allele_vcf_ch
    tuple val(sample_id), path("${sample_id}.common_runs.position.vcf"), emit: common_runs_position_vcf_ch
    tuple val(sample_id), path("*.log"), emit: IntersectRunVCFs_logs


    script:
    def vcf_list = vcf_files instanceof List    ? vcf_files 
                                                : [vcf_files]
    def n_vcfs = vcf_list.size()
    def input_vcfs = vcf_list.collect { it.toString() }.join(' ')

    """
    echo "IntersectRunVCFs" > ${sample_id}_intersect_runs.log
    echo "==============================================================================" >> ${sample_id}_intersect_runs.log
    echo "Sample: ${sample_id}" >> ${sample_id}_intersect_runs.log
    echo "Number of run VCFs: ${n_vcfs}" >> ${sample_id}_intersect_runs.log
    echo "" >> ${sample_id}_intersect_runs.log

    mkdir -p normalized_vcfs

    i=0
    for vcf in ${input_vcfs}
    do
        i=\$((i+1))
        out="normalized_vcfs/run_\${i}.vcf.gz"

        echo "Preparing \$vcf -> \$out" >> ${sample_id}_intersect_runs.log

        if [[ "\$vcf" == *.vcf.gz ]]; then
            cp "\$vcf" "\$out"
        else
            bgzip -c "\$vcf" > "\$out"
        fi

        tabix -f -p vcf "\$out"
    done

    if [ ${n_vcfs} -eq 1 ]; then
        echo "Only one VCF found; no run intersection performed." >> ${sample_id}_intersect_runs.log
        cp normalized_vcfs/run_1.vcf.gz ${sample_id}.common_runs.allele.vcf.gz
    else
        bcftools isec \
            -n=${n_vcfs} \
            -w1 \
            -O z \
            -o ${sample_id}.common_runs.allele.vcf.gz \
            normalized_vcfs/*.vcf.gz >> ${sample_id}_intersect_runs.log 2>&1
    fi

    tabix -f -p vcf ${sample_id}.common_runs.allele.vcf.gz

    bcftools view \
        -O v \
        -o ${sample_id}.common_runs.allele.vcf \
        ${sample_id}.common_runs.allele.vcf.gz

    #copy to keep signature of RefineFilter
    cp ${sample_id}.common_runs.allele.vcf ${sample_id}.common_runs.position.vcf

    echo "" >> ${sample_id}_intersect_runs.log
    echo -n "Variants common to all runs: " >> ${sample_id}_intersect_runs.log
    bcftools view -H ${sample_id}.common_runs.allele.vcf | wc -l >> ${sample_id}_intersect_runs.log

   """
}

process RefineFilter {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*.vcf*"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(intersect_allele_strict_vcf)
    tuple path(reference),     path(reference_idx)
    path (reference_dict)

    output:
    tuple val(sample_id), path("${sample_id}_SNP_final_strict.vcf.gz"), emit: final_vcf
    path "${sample_id}_SNP_final_strict.vcf", emit: final_vcf_uncompressed
    tuple val(sample_id), path("*.log"), emit: RefineFilter_logs


    script:
    """
    PREFIX="${sample_id}"

#keep only SNPs
    gatk --java-options "-Xmx20g" SelectVariants \
        -R $reference \
        -V $intersect_allele_strict_vcf \
        --select-type-to-include SNP \
        -O \${PREFIX}_tmp_snps.vcf >> \${PREFIX}_filter.log 2>&1


#filtering:
#QD: quality by depth
#FS: fisher strand, biais of strand(if only seen in forward, strange)
#MQ: mapping quality: if enough reads mapped on it.
#MappingQualityRankSum: if the mapping quality of the alt allele is much worse than ref
    gatk --java-options "-Xmx20g" VariantFiltration \
        -R $reference \
        -V \${PREFIX}_tmp_snps.vcf \
        -O \${PREFIX}_tmp_gatk_flagged.vcf \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || DP < 10.0 || MQRankSum < -12.5" \
        --filter-name "gatk_hard_filter" >> \${PREFIX}_filter.log 2>&1

#haploid filtering: ratio between ref and alt is > 0.9
#sequencing depth too high (repeat region, duplication)
#GT genotype=heterozygote
#|| GT == "het"


    bcftools filter \
        -e 'FILTER != "PASS" || FORMAT/DP < 4 || FORMAT/DP > 200 ||  FORMAT/AD[0:1]/FORMAT/DP < 0.9' \
        -O z \
        -o \${PREFIX}_SNP_final_strict.vcf.gz \
        \${PREFIX}_tmp_gatk_flagged.vcf >> \${PREFIX}_filter.log 2>&1

    tabix -p vcf \${PREFIX}_SNP_final_strict.vcf.gz

    echo "10. RefineFilter: number of SNPs after filtering:" >> \${PREFIX}_gatk_filter.log
    echo "==============================================================================" >> \${PREFIX}_gatk_filter.log
    echo -n "Total SNPs before filters: " >> \${PREFIX}_gatk_filter.log
    bcftools view -H \${PREFIX}_tmp_snps.vcf | wc -l >> \${PREFIX}_gatk_filter.log
    
    echo -n "Total SNPs conserved (pure haploid + DP valids): " >> \${PREFIX}_gatk_filter.log
    bcftools view -H \${PREFIX}_SNP_final_strict.vcf.gz | wc -l >> \${PREFIX}_gatk_filter.log
    rm -f \${PREFIX}_tmp_snps.vcf \${PREFIX}_tmp_gatk_flagged.vcf
    gunzip -c \${PREFIX}_SNP_final_strict.vcf.gz > \${PREFIX}_SNP_final_strict.vcf    
    """
}


process Logs {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/logs", mode: 'copy'

    input:
    tuple val(sample_id), path(log_files)

    output:
    path("${sample_id}_gatk_full.log"), optional: true
    path("reference_global.log"), optional: true

    script:
    """
    echo "=== Processing Logs for Sample: ${sample_id} ==="
    
    if [ "${sample_id}" = "REF" ]; then
        ls *reference.log >/dev/null 2>&1 && cat *reference.log > reference_global.log || touch reference_global.log
    else
        ls *qc.log >/dev/null 2>&1 && cat *qc.log > entries_qc.log || touch entries_qc.log
        ls *align.log >/dev/null 2>&1 && cat *align.log > entries_align.log || touch entries_align.log
        ls *addrg.log >/dev/null 2>&1 && cat *addrg.log > entries_addrg.log || touch entries_addrg.log
        ls *dedup.log >/dev/null 2>&1 && cat *dedup.log > entries_dedup.log || touch entries_dedup.log

        echo "==========================================================================" > "${sample_id}_gatk_full.log"
        echo "GATK COMPLETE PIPELINE LOG FOR SAMPLE: ${sample_id}" >> "${sample_id}_gatk_full.log"
        echo "==========================================================================" >> "${sample_id}_gatk_full.log"
        
        echo -e "\n--- STEP 1: QC CLEAN ---" >> "${sample_id}_gatk_full.log"
        [ -s entries_qc.log ] && cat entries_qc.log >> "${sample_id}_gatk_full.log" || echo "No QC logs found for this run mode." >> "${sample_id}_gatk_full.log"

        echo -e "\n--- STEP 2: BWA ALIGN ---" >> "${sample_id}_gatk_full.log"
        [ -s entries_align.log ] && cat entries_align.log >> "${sample_id}_gatk_full.log" || echo "No alignment logs found for this run mode." >> "${sample_id}_gatk_full.log"

        echo -e "\n--- STEP 3: ADD READ GROUPS ---" >> "${sample_id}_gatk_full.log"
        [ -s entries_addrg.log ] && cat entries_addrg.log >> "${sample_id}_gatk_full.log" || echo "No read group logs found for this run mode." >> "${sample_id}_gatk_full.log"

        echo -e "\n--- STEP 4: MARK DUPLICATES ---" >> "${sample_id}_gatk_full.log"
        [ -s entries_dedup.log ] && cat entries_dedup.log >> "${sample_id}_gatk_full.log" || echo "No duplication logs found for this run mode." >> "${sample_id}_gatk_full.log"

        echo -e "\n--- STEP 5: VARIANT CALLING ---" >> "${sample_id}_gatk_full.log"
        ls *hcFB.log *fb.log *hc.log *comparison.log >/dev/null 2>&1 && cat *hcFB.log *fb.log *hc.log *comparison.log >> "${sample_id}_gatk_full.log" 2>/dev/null || echo "No variant calling logs found." >> "${sample_id}_gatk_full.log"

        echo -e "\n--- STEP 6: REFINE FILTER ---" >> "${sample_id}_gatk_full.log"
        ls *filter.log >/dev/null 2>&1 && cat *filter.log >> "${sample_id}_gatk_full.log" || echo "No filter logs found." >> "${sample_id}_gatk_full.log"
        
        rm -f entries_qc.log entries_align.log entries_addrg.log entries_dedup.log
    fi
    """
}

