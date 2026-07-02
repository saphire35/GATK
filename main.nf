#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.version='1.0.0'
params.help=null
params.data=null
params.metadata=null
params.outdir='GATK_results'
params.reference=null
params.cpus=4
params.merged_bam=null
params.merged_vcf =null
params.F='2308'
params.f='2'
params.q='20'
//QC
params.galore_quality=30
params.galore_length=100
params.skip_HC=null

params.aligner='BWA'

params.merged_fastq= null

include {QCClean ; MergeFastq ; ReadAlign ; AddReadGroups ; Dedup ; MergeBams ;  HaplotypeCaller ; IndexReference ; HaplotypeCallerFreeBayes ; GenotypeInterval ; ConcatVCFs ; CompareCallers ; IntersectRunVCFs;  RefineFilter ; Logs } from './modules/gatk.nf'

workflow{

        if (params.help) {
                log.info """
                =============================================================================
                G A T K - P I P E L I N E :   V A R I A N T   C A L L I N G ($params.version)
                ============================================================================= 

                Usage:
                nextflow run main.nf --data <DATA_DIR> --metadata <METADATA.tsv> --reference <REF.fasta> [options]

                Required arguments:
                --data              Path to input FASTQ directory
                                Expected format: sample_1.fastq and sample_2.fastq

                --metadata          Path to metadata file used for read groups

                --reference         Path to reference genome FASTA file

                Optional arguments:
                --outdir           Output directory (default: GATK_results)
                --cpus             Total number of CPUs to use (default: 4)
                --merged_bam       If provided, BAM files will be merged before variant calling (default: null).
                --merged_vcf       If provided, VCF files will be intersected before final filtering (default: null).
                --F                BWA filter for proper pairing (default: 2308)
                --f                BWA filter for read group (default: 2)
                --q                BWA filter for mapping quality (default: 20)
                --galore_quality   Trim Galore quality threshold (default: 30)
                --galore_length    Trim Galore minimum length (default: 100)


                Version:
                        ${params.version}

                =============================================================================
                """
                System.exit(0)
        }

        def start_time = new Date()
        log.info("Job started at: ${new Date()}")
        log.info(
                """\

        =============================================================================
        G A T K - P I P E L I N E :   V A R I A N T   C A L L I N 
        =============================================================================
        paths:
        data                                            : ${params.data} (fastq format)
        metadata                                        : ${params.metadata}
        reference genome                                : ${params.reference}
        output directory                                : ${params.outdir}

        merged_bam                                      : ${params.merged_bam}
        merged_vcf                                      : ${params.merged_vcf}
        merged_fastq                                    : ${params.merged_fastq}

        aligner                                         : ${params.aligner}
        skip haplotype caller (--skip_HC)               : ${params.skip_HC}
        galore quality                                  : ${params.galore_quality}
        galore length                                   : ${params.galore_length}
        alignment quality filtering (--q)               : ${params.q}
        alignment check proper pairing filtering (--f)  : ${params.f}
        alignment check read group filtering (--F)      : ${params.F}
        cpus                                            : ${params.cpus}
        =============================================================================
        """
        )
        if (params.data == null) {
                error("Please provide a path to the data directory using --data, file must be in fastq format")
        }
        if (params.metadata == null) {
                error("Please provide a path to the metadata file using --metadata")
        }
        if (params.reference == null) {
                error("Please provide a path to the reference genome using --reference")
        }
        if (params.merged_bam == null && params.merged_vcf == null && params.merged_fastq == null) {
                log.info("No merging option provided, proceeding with individual sample processing")
        }


        //QC with dynamic cpu allocation if //
        def files = file("${params.data}/*_{1,2}.fastq")
        def n_samples = (files.size() / 2) as Integer
        println "[GATK-INFO]: Data channel created with ${n_samples} samples"
        def total = params.cpus as Integer
        def cpus_per_task = Math.max(1, (total / n_samples) as Integer)
        println "[GATK-INFO]: CPUs per task: ${cpus_per_task}"


        data_ch = channel.fromFilePairs("${params.data}/*_{1,2}.fastq", checkIfExists: true)
        data_ch.view { sample_id, reads ->
        "[GATK-INFO]: sample=${sample_id}, reads=${reads}"
        }
        println "[GATK-INFO]: QC analysis with fastqc and trimgalore"

        QCClean(data_ch, cpus_per_task, params.galore_quality, params.galore_length)

        if (params.merged_fastq) {
                MergeFastq(QCClean.out.clean_reads_ch.map { id, reads -> reads }.flatten().collect())
                reads_for_alignment_ch = MergeFastq.out.merged_reads_ch
        } 
        else {
                reads_for_alignment_ch = QCClean.out.clean_reads_ch
        }

        println "[GATK-INFO]: Aligning reads to reference genome"
        IndexReference(file(params.reference))
        ReadAlign(reads_for_alignment_ch, IndexReference.out.reference_with_index, cpus_per_task, params.F,params.f, params.q)


        //adding read groups
        if (params.merged_fastq){
                println "[GATK-INFO]: Adding read groups to merged aligned BAM files using first line of metadata file"
        }
        else{
                println "[GATK-INFO]: Adding read groups to aligned BAM files"
        }
        rg_ch= AddReadGroups(ReadAlign.out.aligned_bams_ch, file(params.metadata))


        //MergeBams
        if (params.merged_bam != null) {
                println "[GATK-INFO]: Merging BAM files with Picard MergeSamFiles"
                merged_ch = MergeBams(rg_ch.sorted_bams_ch.map { _sample_id, bam -> bam }.collect(),
                                        file(params.metadata))

                println "[GATK-INFO]: Marking duplicates with Picard MarkDuplicates"
                Dedup(merged_ch.merged_bam.map { bam -> tuple("merged_sample", bam) }, file(params.metadata))
        }
        else{
                println "[GATK-INFO]: Marking duplicates with Picard MarkDuplicates"
                Dedup(rg_ch.sorted_bams_ch, file(params.metadata))
        }



        //SNP calling
        HaplotypeCallerFreeBayes(Dedup.out.dedup_bams_ch, IndexReference.out.reference_with_index, IndexReference.out.reference_dict, cpus_per_task     )

        if (params.skip_HC==null){
                HaplotypeCaller(Dedup.out.dedup_bams_ch, IndexReference.out.reference_with_index, IndexReference.out.reference_dict, cpus_per_task     )

                if (params.merged_vcf){
                        intervals_ch = IndexReference.out.reference_fai
                                .splitText()
                                .map { line -> line.trim().tokenize('\t')[0] }

                        all_gvcfs_ch = HaplotypeCaller.out.gatk_gvcf_ch
                                .collect(flat: false)
                                .map { records ->
                                        def gvcfs = records.collect { rec -> rec[1] }
                                        def tbis  = records.collect { rec -> rec[2] }

                                        tuple(gvcfs, tbis)
                                }
                        genotype_input_ch = intervals_ch
                                .combine(all_gvcfs_ch)
                                .combine(IndexReference.out.reference_with_index2)
                                .map { interval, gvcfs, tbis, reference, ref_indexes, dict ->
                                        def all_ref_files = [reference] + ref_indexes + [dict]
                                        tuple(interval, "merged_vcf", gvcfs, tbis, all_ref_files)
                                }
                }
                else{
                        println "[GATK-INFO]: Running HaplotypeCaller and GenotypeGVCFs on individual BAM files"
                        intervals_ch =IndexReference.out.reference_fai.splitText()
                                .map { line -> line.trim().tokenize('\t')[0] }

                        genotype_input_ch= intervals_ch
                                .combine(HaplotypeCaller.out.gatk_gvcf_ch)
                                .combine(IndexReference.out.reference_with_index2)
                                .map { interval, sample_id, gvcf, tbi, reference, ref_indexes, dict ->
                                        def all_ref_files = [reference] + ref_indexes + [dict]
                                        tuple(interval, sample_id, gvcf, tbi, all_ref_files)
                                }
                }
                GenotypeInterval(genotype_input_ch)
                concat_input = GenotypeInterval.out.groupTuple()
                ConcatVCFs(concat_input)
                compare_input_ch = ConcatVCFs.out.map { tuple(it[0], it[1]) }
                CompareCallers( compare_input_ch, HaplotypeCallerFreeBayes.out.freebayes_vcf_ch )
                refine_input_ch = CompareCallers.out.intersect_allele_strict_vcf
        }
        else{
                refine_input_ch= HaplotypeCallerFreeBayes.out.freebayes_vcf_ch
        }

        //merged_vcf--> intersect before refinefilter
        if (params.merged_vcf != null ){
                println "[GATK-INFO]: filtering the vcf between the runs"
                run_intersection_input_ch = refine_input_ch
                        .map { sample_id, vcf ->
                                def common_sample_id = sample_id.replaceAll(/_[0-9]+_libLAO[0-9]+$/, '')
                                tuple(common_sample_id, vcf)
                        }
                        .groupTuple()
                println "[GATK-INFO]: merging the vcf between the runs ${run_intersection_input_ch}"
                IntersectRunVCFs(run_intersection_input_ch)

                RefineFilter(IntersectRunVCFs.out.common_runs_allele_vcf_ch,
                IndexReference.out.reference_with_index,IndexReference.out.reference_dict)

                refine_logs = IntersectRunVCFs.out.IntersectRunVCFs_logs
                .mix(RefineFilter.out.RefineFilter_logs)  
        }
        else {
                println "[GATK-INFO]: filtering but not merging the vcf runs"
                RefineFilter( refine_input_ch,  
                        IndexReference.out.reference_with_index, IndexReference.out.reference_dict )
                refine_logs=RefineFilter.out.RefineFilter_logs             
        }



        all_logs_grouped = QCClean.out.QCClean_logs
        .mix(
                ReadAlign.out.ReadAlign_logs,
                AddReadGroups.out.AddReadGroups_logs,
                Dedup.out.Dedup_logs,
                HaplotypeCallerFreeBayes.out.HaplotypeCallerFreeBayes_logs,
                refine_logs,
                CompareCallers.out.CompareCallers_logs
        )
        .map { id, logs ->
                def common_id = id.replaceAll(/_[0-9]+_libLAO[0-9]+$/, '')
                tuple(common_id, logs)
        }
        .groupTuple()
        .map { id, all_logs -> [ id, all_logs.flatten() ] }

        Logs(all_logs_grouped)



        workflow.onComplete {
                def end_time = new Date()
                def duration_ms = end_time.time - start_time.time

                //convertion h/min/s
                def duration_sec = (duration_ms / 1000) as Integer
                def hours = (duration_sec / 3600) as Integer
                def minutes = ((duration_sec % 3600) / 60) as Integer
                def seconds = duration_sec % 60

                log.info("Job finished at: ${end_time}")
                log.info("Total runtime: ${hours}h ${minutes}m ${seconds}s")
        }


}