#!/usr/bin/env nextflow

include { freebayes; bcftools_stats } from '../modules/variant_calling.nf'
include { multiqc } from '../modules/qc.nf'
include { logUnmapped } from '../modules/logs.nf'

include { qc; mapping; sv } from '../modules/subworkflow.nf'
include { mapping as mapping_ref } from '../modules/subworkflow.nf'

include { calc_unmapped; bwa_index; get_unmapped_reads } from '../modules/mapping.nf'
include { bwa_index as bwa_index_ref } from '../modules/mapping.nf'


out_folder_name = "short-mod"

workflow mod_ref {
    take:
        trimmed
        ref_fasta
        mod_fasta

    main:
        // map to the reference fasta
        bwa_index_ref(ref_fasta, out_folder_name) | set { ref_index } 
        mapping_ref(ref_fasta, ref_index, trimmed, out_folder_name) | set { indexed_bam }

        // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.short_threshold, out_folder_name)

        // extract unmapped reads
        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }

        // map to modified fasta
        bwa_index(mod_fasta, out_folder_name) | set { mod_index } 
        mapping(mod_fasta, mod_index, unmapped_fastq, out_folder_name) | set { indexed_unmapped_bam }

        // SNPs variant calling
        freebayes(mod_fasta, mod_index, indexed_unmapped_bam, out_folder_name) | set { vcf }
        bcftools_stats(vcf, out_folder_name) | set { bcftools_out }      
        
        multiqc(bcftools_out, out_folder_name, 'varint_calling')

        // SVs variant calling
        sv(mod_fasta, indexed_unmapped_bam, out_folder_name)
        
        emit:
            log.info "▶ The ${out_folder_name} processing pipeline completed successfully."
}



workflow { 
    log.info  "Processing files in directory: ${params.in_dir}"
    
    Channel.fromPath("$params.in_dir/*mod.{fa,fna,fasta}") | set { mod_fasta }
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { ref_fasta }

    
    Channel.fromPath("$params.in_dir/*.fastq.gz")
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    qc(fastqs, out_folder_name) | set { trimmed }

    mod_ref(trimmed, ref_fasta, mod_fasta)
}