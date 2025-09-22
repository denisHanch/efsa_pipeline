#!/usr/bin/env nextflow

include { freebayes; bcftools_stats } from '../modules/variant_calling.nf'
include { multiqc } from '../modules/qc.nf'
include { qc; mapping; sv; annotate_vcf } from '../modules/subworkflow.nf'
include { calc_unmapped; bwa_index } from '../modules/mapping.nf'
include { logUnmapped } from '../modules/logs.nf'


out_folder_name = "short-ref"

workflow short_ref {
    take:
        trimmed
        fasta

    main:
        bwa_index(fasta, out_folder_name) | set { fasta_index } 
        mapping(fasta, fasta_index, trimmed, out_folder_name) | set { indexed_bam }

        // printout % unmapped reads
        calc_unmapped(indexed_bam) | set { pct }
        logUnmapped(pct, params.short_threshold, out_folder_name)

        // SNPs variant calling
        freebayes(fasta, fasta_index, indexed_bam, out_folder_name) | set { vcf }
        bcftools_stats(vcf, out_folder_name) | set { bcftools_out }
        
        // Annotate vcf
        Channel.fromPath("$params.in_dir/*ref.gtf", checkIfExists: true) | set { gtf }

        annotate_vcf(fasta, gtf, vcf) | set {qc_vcf}
        
        qc_vcf.mix(bcftools_out).collect() | set { qc_out }
        multiqc(qc_out, out_folder_name, 'varint_calling')

        // SVs variant calling
        sv(fasta, indexed_bam, out_folder_name)
        
        emit:
            log.info "▶ The ${out_folder_name} processing pipeline completed successfully."
}


workflow { 
    log.info  "Processing files in directory: ${params.in_dir}"
    
    Channel.fromPath("$params.in_dir/*.fastq.gz") | set { fastqs }
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { fasta }
    
    fastqs
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    qc(fastqs, out_folder_name) | set { trimmed }

    short_ref(trimmed, fasta)
}