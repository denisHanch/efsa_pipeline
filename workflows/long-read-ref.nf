#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { sv_long; mapping_long }  from '../modules/subworkflow.nf'
include { logUnmapped } from '../modules/logs.nf'
include { calc_unmapped } from '../modules/mapping.nf'

out_folder_name = "long-ref"

workflow long_ref {
    take:
        fastqs
        fasta

    main:
        // qc
        nanoplot(fastqs, out_folder_name)
        
        // mapping
        mapping_long(fastqs, fasta, out_folder_name) | set { indexed_bam }

         // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.long_threshold, out_folder_name)

        // variant calling
        sv_long(fasta, indexed_bam, out_folder_name)

    emit:
        log.info "▶ The long read processing pipeline completed successfully."
    }


workflow {
    // Processing inputs
    log.info  "Processing files in directory: ${params.in_dir}"

    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { fasta }

    Channel.fromPath("${params.in_dir}/tmp2/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { fastqs }
    
    long_ref(fastqs, fasta)
}
