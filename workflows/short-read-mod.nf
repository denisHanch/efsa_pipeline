#!/usr/bin/env nextflow

include { freebayes; bcftools_stats } from '../modules/variant_calling.nf'
include { multiqc } from '../modules/qc.nf'
include { qc; mapping; sv } from '../modules/subworkflow.nf'
include { calc_unmapped; bwa_index } from '../modules/mapping.nf'
include { logUnmapped } from '../modules/logs.nf'

out_folder_name = "short-mod"

workflow mod_ref {
    take:
        trimmed
        fasta

    main:
        bwa_index(fasta, out_folder_name) | set { fasta_index } 
        mapping(fasta, fasta_index, trimmed, out_folder_name) | set { indexed_bam }

        // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, out_folder_name)

        // SNPs variant calling
        freebayes(fasta, fasta_index, indexed_bam, out_folder_name) | set { vcf }
        bcftools_stats(vcf, out_folder_name) | set { bcftools_out }      
        
        multiqc(bcftools_out, out_folder_name, 'varint_calling')

        // SVs variant calling
        sv(fasta, indexed_bam, out_folder_name)
        
        emit:
            log.info "▶ The ${out_folder_name} processing pipeline completed successfully."
}



workflow { 
    log.info  "Processing files in directory: ${params.in_dir}"
    
    Channel.fromPath("$params.in_dir/*mod.{fa,fna,fasta}") | set { fasta }
    
    Channel.fromPath("$params.in_dir/*.fastq.gz")
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    qc(fastqs, out_folder_name) | set { trimmed }

    mod_ref(trimmed, fasta)
}