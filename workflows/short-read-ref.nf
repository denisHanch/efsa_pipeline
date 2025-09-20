#!/usr/bin/env nextflow

include { convert_bcf_to_vcf; delly; samtools_index; picard_dict } from '../modules/sv_calling.nf'
include { freebayes; snpeff; build_config; bcftools_stats } from '../modules/variant_calling.nf'
include { fastqc; multiqc; trimgalore } from '../modules/qc.nf'
include { bwa_index; bwa_mapping; samtool_index_bam; samtools_sort; samtool_stats; picard; calc_unmapped } from '../modules/mapping.nf'
include { short_mod } from './short-read-mod.nf'

out_folder_name = "short-ref"


workflow short_ref {

    // Processing inputs
    log.info  "Processing files in directory: ${params.in_dir}"
    
    Channel.fromPath("$params.in_dir/*.fastq.gz") | set { fastqs }
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { fasta }
    Channel.fromPath("$params.in_dir/*ref.gtf") | set { gtf }
    
    fastqs
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    trimgalore(fastqs, out_folder_name) | set { trimmed }
    fastqc(trimmed, out_folder_name) | set {fastqc_out}

    // Mapping module
    bwa_index(fasta, out_folder_name) | set { fasta_index } 
    bwa_mapping(fasta, fasta_index, trimmed, out_folder_name) | set { sam } 
    samtools_sort(sam, out_folder_name) | set { bam }
    samtool_stats(bam, out_folder_name) | set { stats_out }
    samtool_index_bam(bam, out_folder_name) | set { indexed_bam }
    picard(fasta, indexed_bam, out_folder_name) | set { picard_out }

    // SNPs variant calling
    freebayes(fasta, fasta_index, indexed_bam, out_folder_name) | set { vcf }
    build_config(fasta, gtf) | set { snpeff_out }
    genome_id = snpeff_out.map {genome_id, snpeff_config -> genome_id}
    snpeff_config = snpeff_out.map {genome_id, snpeff_config -> snpeff_config}
    snpeff(vcf, genome_id, snpeff_config) | set { snpeff_output }
    bcftools_stats(vcf, out_folder_name) | set { bcftools_out }
    annotated_vcfs = snpeff_output.map { id, vcf, html -> tuple(id, vcf) }
    qc_vcf = snpeff_output.map { id, vcf, html -> html }

    // // SVs variant calling
    samtools_index(fasta, out_folder_name) | set { fai }
    picard_dict(fasta, out_folder_name) | set { dict }
    delly(indexed_bam, fasta, fai, dict, out_folder_name) | set { bcf }
    convert_bcf_to_vcf(bcf, out_folder_name) | set { sv_vcf }

    // running multiqc on all files
    fastqc_out.mix(stats_out).mix(picard_out).mix(qc_vcf).mix(bcftools_out).collect() | set { qc_out }
    multiqc(out_folder_name, qc_out)

    log.info "▶ The ${out_folder_name} processing pipeline completed successfully."
}


workflow { 
    short_ref()
}