#!/usr/bin/env nextflow

include { nanoplot; multiqc } from '../modules/qc.nf'
include { calc_unmapped; get_unmapped_reads } from '../modules/mapping.nf'
include { sv_long; mapping_long; mapping_long as mapping_long_mod }  from '../modules/subworkflow.nf'
include { logUnmapped } from '../modules/logs.nf'


out_folder_name = "long-mod"

workflow long_mod {
    take:
        fastqs
        ref_fasta
        mod_fasta
        mapping_tag

    main:
        // qc
        nanoplot(fastqs, out_folder_name)
        
        // mapping to reference fasta
        mapping_long(fastqs, ref_fasta, mapping_tag, out_folder_name) | set { indexed_bam }

         // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.long_threshold,  "long-ref-${mapping_tag}")

        // extract unmapped reads
        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }

        // mapping to modified fasta
        mapping_long_mod(unmapped_fastq, mod_fasta, mapping_tag, out_folder_name) | set { indexed_bam }

        // variant calling
        sv_long(mod_fasta, indexed_bam, out_folder_name)

    emit:
        log.info "▶ The long read processing pipeline completed successfully."
    }


workflow {
    // Processing inputs
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*mod.{fa,fna,fasta}") | set { mod_fasta }

    Channel.fromPath("${params.in_dir}/tmp2/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { fastqs }
    
    long_mod(fastqs, ref_fasta, mod_fasta, mapping_tag)
}
