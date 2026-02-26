#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.version='1.0.0'
params.data=null
params.metadata=null
params.outdir='GATK_results'
params.reference=null
params.cpus=4

//QC
params.galore_quality=30
params.galore_length=100


//
params.merged=null
params.dataset = "GATK"

include {QCClean ; BWAAlign ; AddReadGroups ; Dedup ; MergeBams } from './modules/gatk.nf'
include { HaplotypeCaller ; IndexReference ; HaplotypeCallerFreeBayes ; GenotypeInterval ; ConcatVCFs ; CompareCallers ; RefineFilter ; Logs } from './modules/gatk.nf'





workflow{
        def start_time = new Date()
        log.info("Job started at: ${new Date()}")
        log.info(
                """\

        =============================================================================
        G A T K - P I P E L I N E :   V A R I A N T   C A L L I N G  ${params.version}
        =============================================================================
        paths:
        data                    : ${params.data}
        metadata                : ${params.metadata}
        reference genome        : ${params.reference}
        output directory        : ${params.outdir}
        dataset                 : ${params.dataset}

        merged                  : ${params.merged}

        galore quality          : ${params.galore_quality}
        galore length           : ${params.galore_length}
        cpus                    : ${params.cpus}
        =============================================================================
        """
        )
        if (params.data == null) {
                error("Please provide a path to the data directory using --data")
        }
        if (params.metadata == null) {
                error("Please provide a path to the metadata file using --metadata")
        }
        if (params.reference == null) {
                error("Please provide a path to the reference genome using --reference")
        }


        //QC with dynamic cpu allocation if //
        def files = file("${params.data}/*_{1,2}.fastq")
        def n_samples = (files.size() / 2) as Integer
        println "[GATK-INFO]: Data channel created with ${n_samples} samples"

        def total = params.cpus as Integer
        def cpus_per_task = Math.max(1, (total / n_samples) as Integer)
        println "[GATK-INFO]: CPUs per task: ${cpus_per_task}"

        println "[GATK-INFO]: QC analysis with fastqc and trimgalore"
        data_ch = channel.fromFilePairs("${params.data}/*_{1,2}.fastq", checkIfExists: true)
        QCClean(data_ch, cpus_per_task, params.galore_quality, params.galore_length)

        //BWA Alignment and reference indexing
        println "[GATK-INFO]: Aligning reads to reference genome with BWA-MEM"
        IndexReference(file(params.reference))
        BWAAlign(QCClean.out.clean_reads_ch, IndexReference.out.reference_with_index, cpus_per_task)


        //adding read groups
        println "[GATK-INFO]: Adding read groups to aligned BAM files"
        rg_ch= AddReadGroups(BWAAlign.out.aligned_bams_ch, file(params.metadata))


        //MergeBams
        if (params.merged != null) {
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
        HaplotypeCaller(Dedup.out.dedup_bams_ch, IndexReference.out.reference_with_index, IndexReference.out.reference_dict, cpus_per_task     )
        HaplotypeCallerFreeBayes(Dedup.out.dedup_bams_ch, IndexReference.out.reference_with_index, IndexReference.out.reference_dict, cpus_per_task     )
      //  ConsolidateGenotypes(HaplotypeCaller.out.gvcf_ch.map { _sample_id, gvcf -> gvcf }.collect(),
        //                        IndexReference.out.reference_with_index,IndexReference.out.reference_dict, cpus_per_task)



        def fai_file = file("${params.reference}.fai")
        def n_intervals = fai_file.readLines().size()
        def cpu_total = params.cpus as Integer
        def cpus_per_interval = Math.max(1, (cpu_total / n_intervals) as Integer)
        println "[GATK-INFO]: Number of intervals = ${n_intervals}"
        println "[GATK-INFO]: CPUs per interval = ${cpus_per_interval}"
        intervals_ch = channel.fromList(  file("${params.reference}.fai")
                .readLines()
                .collect { line -> line.tokenize('\t')[0] }
        )

        genotype_input_ch= intervals_ch
                                .combine(HaplotypeCaller.out.gatk_gvcf_ch)
                                .combine(IndexReference.out.reference_with_index2)
                                .map { interval, _sample_id, gvcf, tbi, reference, ref_indexes, dict ->
                                        def all_ref_files = [reference] + ref_indexes + [dict]
                                        tuple(interval, gvcf, tbi, all_ref_files)
                                }
        GenotypeInterval(genotype_input_ch)
        ConcatVCFs(GenotypeInterval.out.collect())
  
       

        CompareCallers( ConcatVCFs.out.gatk_gvcf, HaplotypeCallerFreeBayes.out.freebayes_vcf_ch )
        RefineFilter(CompareCallers.out.intersect_allele_strict_vcf, CompareCallers.out.intersect_position_vcf, IndexReference.out.reference_with_index, IndexReference.out.reference_dict)
        Logs(QCClean.out.QCClean_logs, BWAAlign.out.BWAAlign_logs, IndexReference.out.IndexReference_logs , AddReadGroups.out.AddReadGroups_logs, Dedup.out.Dedup_logs,
                HaplotypeCallerFreeBayes.out.HaplotypeCallerFreeBayes_logs, RefineFilter.out.RefineFilter_logs)





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
