
process QCClean {

    tag "$sample_id"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val cpus_per_task_ch
    val galore_quality 
    val galore_length

    output:
    tuple val(sample_id), path("trimmed/*_val_*.fq"), emit: clean_reads_ch
    path "QC_reports/*"
    path "*.log", emit: QCClean_logs

    script:
    """
    mkdir -p fastqc_raw fastqc_clean trimmed

    read1=${reads[0]}
    read2=${reads[1]}

    echo "Processing $sample_id"
    echo "R1: \$read1"
    echo "R2: \$read2"

    # FastQC raw + Trim Galore
    fastqc \$read1 \$read2 -o fastqc_raw -t ${cpus_per_task_ch} &
    trim_galore -q ${galore_quality} --length ${galore_length} --illumina --paired \$read1 \$read2 -o trimmed 

    # FastQC trimmed
    fastqc trimmed/*_val_1.fq trimmed/*_val_2.fq -o fastqc_clean -t ${cpus_per_task_ch}
    
    
    echo "1. QCClean \n" >> ${sample_id}_gatk.log
    echo "============================================================================== \n" >> ${sample_id}_gatk.log
    grep -E "Total reads processed|Reads written|Number of sequence pairs removed"   .command.log >> ${sample_id}_gatk.log
    grep -A2 "Processing " .command.log  | grep -v -A1  "Processing reads" >> ${sample_id}_gatk_qc.log


    mkdir QC_reports
    cp trimmed/*report.txt fastqc_clean/*.html fastqc_raw/*.html QC_reports/
    """
}



process IndexReference {

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: reference_with_index
    path "*.dict", emit: reference_dict
    tuple path(reference), path("${reference}.*"), path("*.dict"), emit: reference_with_index2
    path "*log", emit: IndexReference_logs

    script:
    """
    bwa index $reference
    samtools faidx $reference
    java -jar  /usr/local/bin/picard.jar CreateSequenceDictionary \
        R=$reference \
        O=${reference.baseName}.dict

    echo "2. IndexReference \n" >> gatk_reference.log
    echo "============================================================================== \n" >> gatk_reference.log
    echo "number of bases: \$(awk '{sum += \$2} END {print sum}' *.fai)" >> gatk_reference.log
    echo "number of contigs: \$(wc -l < *.fai)" >> gatk_reference.log
    
    """
}

process BWAAlign{
    tag "$sample_id"


    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(index_files)
    val cpus_per_task_ch

    output:
    tuple val(sample_id), path ("${sample_id}.sorted.raw.bam*"), path ("${sample_id}.mappingstats.txt"), emit: aligned_bams_ch
    tuple path ("*sorted.raw.bam*"), path ("*sorted.raw.bam.bai"), emit: reference_index_ch
    path "*.log", emit: BWAAlign_logs


    script:
    """

    read1=${reads[0]}
    read2=${reads[1]}   
    bwa mem -t "$cpus_per_task_ch" "$reference" "\$read1" "\$read2" \
        |samtools view -bS - \
        |samtools sort -@ "$cpus_per_task_ch" -o "${sample_id}.sorted.raw.bam"
    samtools index "${sample_id}.sorted.raw.bam"
    samtools flagstat "${sample_id}.sorted.raw.bam" > "${sample_id}.mappingstats.txt"
    samtools quickcheck "${sample_id}.sorted.raw.bam"


    echo "3. BWAAlign \n" >> ${sample_id}_gatk_align.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_align.log
    cat *.mappingstats.txt >> ${sample_id}_gatk_align.log
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
    path "*log", emit: AddReadGroups_logs

    script:
    """
    #retrieve metadata for the sample
    #order: simple_ID sample_ID read1 read2 instrument flowcell lane barcode sex run_num seq_num
    line=\$(grep $sample_id $metadata)
    simpleID=\$(echo \$line | cut -f1 -d' ')
    instrument=\$(echo \$line | cut -f5 -d' ')
    flowcell=\$(echo \$line | cut -f6 -d' ')
    lane=\$(echo \$line | cut -f7 -d' ')
    barcode=\$(echo \$line | cut -f8 -d' ')
    seqnum=\$(echo \$line | cut -f11 -d' ')

    bam_file=${bam_files[0]}
    echo "launching AddOrReplaceReadGroups for $sample_id with bam file \$bam_file"
    java -Xmx10g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups \
        I=\$bam_file \
        O=${sample_id}.sorted.rg.bam \
        RGSM=\$simpleID \
        RGLB=\${simpleID}.\${seqnum} \
        RGID=\${flowcell}.\${lane} \
        RGPU=\${flowcell}\${lane}.\${barcode} \
        RGPL=\$instrument

    # Index
    java -Xmx10g -jar /usr/local/bin/picard.jar BuildBamIndex I=${sample_id}.sorted.rg.bam


    echo "4. AddReadGroups \n" >> ${sample_id}_gatk_addrg.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_addrg.log
    grep -E "AddOrReplaceReadGroups.*Created"  .command.log >> ${sample_id}_gatk_addrg.log

    """
}

process Dedup {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file)
    path metadata

    output:
    tuple val(sample_id), path("${sample_id}.sorted.dedup.bam"), path("${sample_id}.sorted.dedup.bai"), emit: dedup_bams_ch
    path "*.log", emit: Dedup_logs

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


    echo "5. Dedup \n" >> ${sample_id}_gatk_dedup.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_dedup.log
    grep -A100 "## METRICS CLASS" ${sample_id}.dedup.metrics.txt >> ${sample_id}_gatk_dedup.log
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
    publishDir "${params.outdir}", mode: 'copy'    


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
    """
    gatk --java-options "-Xmx10g" HaplotypeCaller \
        -R "$reference" \
        -I "$bam_file" \
        -O "${sample_id}.g.vcf.gz" \
        -ERC GVCF \
        --sample-ploidy 1


    mkdir -p "${params.outdir}"
    cp "${sample_id}.g.vcf.gz" "${params.outdir}/${sample_id}.g.vcf.gz"
    """

}

process HaplotypeCallerFreeBayes {
    publishDir "${params.outdir}", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    tuple path(reference), path(reference_idx)
    path reference_dict
    val cpus_per_task_ch

    output:
    path("${sample_id}.freebayes.vcf"), emit: freebayes_vcf_ch
    path "*log", emit: HaplotypeCallerFreeBayes_logs

    script:
    """
    freebayes-v1.3.1 \
        -f "$reference" \
        -p 1 \
        --min-alternate-count 3 \
        --min-alternate-fraction 0.2 \
        --pooled-discrete \
        "$bam_file" \
        > "${sample_id}.freebayes.vcf"

    echo "7. freebayes variant number for $sample_id:" >> ${sample_id}_gatk_hcFB.log
    echo "============================================================================== \n" >> ${sample_id}_gatk_hcFB.log
    grep -v "^#" "${sample_id}.freebayes.vcf" | wc -l >> ${sample_id}_gatk_hcFB.log
    mkdir -p "${params.outdir}"
    cp "${sample_id}.freebayes.vcf" "${params.outdir}/${sample_id}.freebayes.vcf"
    """
}


process GenotypeInterval {

    tag "$interval"

    input:
    tuple val(interval),
          path(gvcf),
          path(tbi),
          path(reference_files)

    output:
    tuple path("${interval}.vcf.gz"),
          path("${interval}.vcf.gz.tbi")
    script:

    def reference = reference_files.find { file -> file.name.endsWith(".fasta") }

    """
    gatk GenomicsDBImport \
        --variant ${gvcf} \
        --genomicsdb-workspace-path INTERVAL_${interval}_db \
        --intervals ${interval}

    gatk GenotypeGVCFs \
        -R ${reference} \
        -V gendb://INTERVAL_${interval}_db \
        -O ${interval}.vcf.gz

    rm -rf INTERVAL_${interval}_db
    """
}


process ConcatVCFs {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(vcf_files)

    output:
    path "${params.dataset}_cohort_batch_genotyped.g.vcf.gz", emit: gatk_gvcf
    path "${params.dataset}_cohort_batch_genotyped.g.vcf.gz.tbi"
    path "batch_filter.txt"

    script:
    def vcfs = vcf_files.findAll { file -> file.name.endsWith(".vcf.gz") }.sort { file -> file.name }

    """
    # SLURM-style: create a batch_filter.txt then concat with -f
    printf "%s\\n" ${vcfs.join(' ')} | tr ' ' '\\n' > batch_filter.txt

    bcftools concat \
        -a \
        -Oz \
        -o ${params.dataset}_cohort_batch_genotyped.g.vcf.gz \
        -f batch_filter.txt

    tabix -p vcf ${params.dataset}_cohort_batch_genotyped.g.vcf.gz
    mkdir -p ${params.outdir}
    cp ${params.dataset}_cohort_batch_genotyped.g.vcf.gz ${params.outdir}/${params.dataset}_cohort_batch_genotyped.g.vcf.gz
    """
}




process CompareCallers {

    publishDir "${params.outdir}", mode: 'copy'
    tag "gatk_vs_freebayes"

    input:
    path gatk_gvcf
    path freebayes_vcf

    output:
    path "gatk.raw.vcf"
    path "freebayes.sorted.vcf.gz", emit: freebayes_vcf
    path "gatk.sorted.vcf.gz", emit: gatk_vcf
    path "intersect.position.vcf", emit: intersect_position_vcf
    path "intersect.allele.strict.vcf", emit: intersect_allele_strict_vcf

    script:
    """
    #Extract variants from GATK gVCF
    bcftools view -v snps,indels "$gatk_gvcf" > gatk.raw.vcf

    #Sort
    bcftools sort gatk.raw.vcf -o gatk.sorted.vcf
    bcftools sort "$freebayes_vcf" -o freebayes.sorted.vcf

    #Compress
    bgzip -f gatk.sorted.vcf
    bgzip -f freebayes.sorted.vcf

    #Index
    tabix -f -p vcf gatk.sorted.vcf.gz
    tabix -f -p vcf freebayes.sorted.vcf.gz

    #Position-based intersect (bedtools)
    bedtools intersect -header \
        -a gatk.sorted.vcf.gz \
        -b freebayes.sorted.vcf.gz \
        > intersect.position.vcf

    # Strict allele intersection (bcftools)
    bcftools isec -n=2 -w1 \
        gatk.sorted.vcf.gz \
        freebayes.sorted.vcf.gz \
        -o intersect.allele.strict.vcf


    echo "9. CompareCallers: number of intersected SNPs:" >> ${params.dataset}_gatk_comparison.log
    echo "============================================================================== \n" >> ${params.dataset}_gatk_comparison.log
    echo "Position-based intersection: \n" >> ${params.dataset}_gatk_comparison.log
    grep -v "^#" intersect.allele.strict.vcf | wc -l >> ${params.dataset}_gatk_comparison.log
    echo "Allele-based strict intersection: \n" >> ${params.dataset}_gatk_comparison.log
    grep -v "^#" intersect.position.vcf | wc -l >> ${params.dataset}_gatk_comparison.log
    """
}


process RefineFilter {
    publishDir "${params.outdir}", mode: 'copy'
    tag "refine_filter"

    input:
    path intersect_allele_strict_vcf
    path intersect_position_vcf
    tuple path(reference), path(reference_idx)
    path reference_dict

    output:
    path "*_SNP*.vcf"
    path "*log", emit: RefineFilter_logs
   

    script:
    """
    PREFIX="${params.dataset}"

    SNP_filtered=\${PREFIX}_SNP_filter.vcf
    gatk_filter_flag=\${PREFIX}_SNP_gatk_flagged.vcf
    gatk_filtered=\${PREFIX}_SNP_gatk_filtered
    allele_filtered=\${PREFIX}_SNP.minmax2.mindp5maxdp200.filtered

    ####################################################
    # 1. Select SNPs (flagged)
    ####################################################
    gatk --java-options "-Xmx20g" SelectVariants \
        -R $reference \
        -V $intersect_allele_strict_vcf \
        --select-type-to-include SNP \
        -O \$SNP_filtered

    gatk --java-options "-Xmx20g" VariantFiltration \
        -R $reference \
        -V \$SNP_filtered \
        -O \$gatk_filter_flag \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5" \
        --filter-name "snp_filter"

    ####################################################
    # 2. Remove filtered variants
    ####################################################
    vcftools \
        --vcf \$gatk_filter_flag \
        --recode \
        --remove-filtered-all \
        --out \$gatk_filtered

    ####################################################
    # 3. Keep biallelic SNPs + depth filter
    ####################################################
    vcftools \
        --vcf \$gatk_filtered.recode.vcf \
        --min-alleles 2 \
        --max-alleles 2 \
        --minDP 4 \
        --maxDP 200 \
        --recode \
        --remove-filtered-all \
        --out \$allele_filtered

    echo "Filtering done"   
    echo "10. RefineFilter: number of SNPs after filtering:" >> ${params.dataset}_gatk_filter.log
    echo "============================================================================== \n" >> ${params.dataset}_gatk_filter.log
    grep "After filtering" .command.log >> ${params.dataset}_gatk_filter.log


    mkdir -p ${params.outdir} 
    cp *_SNP_*.vcf ${params.outdir}/ 
    """
}


process Logs {
    publishDir "${params.outdir}/logs", mode: 'copy'
    tag "logs"

    input:
    path QC_logs_ch
    path BWA_logs_ch
    path IndexReference_logs
    path AddReadGroups_logs
    path Dedup_logs
    path HaplotypeCallerFreeBayes_logs
    path RefineFilter_logs

    output:
    path "*_gatk_qc.log"
    path "*_gatk_align.log"
    path "*_gatk_addrg.log"
    path "*_gatk_dedup.log"

    script:
    """
    cat ${QC_logs_ch} > all_samples_gatk_qc.log
    cat ${BWA_logs_ch} > all_samples_gatk_align.log
    cat ${AddReadGroups_logs} > all_samples_gatk_addrg.log
    cat ${Dedup_logs} > all_samples_gatk_dedup.log

    cat all_samples_gatk_qc.log $IndexReference_logs all_samples_gatk_align.log  all_samples_gatk_addrg.log all_samples_gatk_dedup.log $HaplotypeCallerFreeBayes_logs $RefineFilter_logs >> all_samples_gatk_full.log
    """

}