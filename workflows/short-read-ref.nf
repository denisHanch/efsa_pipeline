#!/usr/bin/env nextflow

include { convert_bcf_to_vcf; delly; samtools_index; picard_dict; svviz } from '../modules/sv_calling.nf'
include { freebayes; snpeff; build_config; bcftools_stats } from '../modules/variant_calling.nf'
include { fastqc; multiqc; trimgalore } from '../modules/qc.nf'
include { bwa_index; bwa_mapping; samtool_index_bam; samtool_stats; picard } from '../modules/mapping.nf'


workflow {
    // Processing inputs
    if (!params.fastqDir) {
            throw new IllegalArgumentException('--fastqDir must be provided')
        } else {

        println("Processing  files in directory: ${params.fastqDir}")
        Channel.fromPath("$params.fastqDir/*.fastq.gz") | set { fastqs }
        }
    Channel.fromPath("$params.fastqDir/ref.fa") | set { fasta }
    Channel.fromPath("$params.fastqDir/ref.gtf") | set { gtf }
    fastqs
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
    | groupTuple(sort: true)
    | set { fastqs }

    // QC and trimming module
    fastqs | trimgalore | set { trimmed }
    fastqc(trimmed) | set {fastqc_out}
    // trimmed.view()

    // Mapping module
    bwa_index(fasta) | set { fasta_index }
    bwa_mapping(fasta, fasta_index, trimmed) | set { bam }
    samtool_stats(bam) | set { stats_out }
    samtool_index_bam(bam) | set { indexed_bam }
    picard(fasta, indexed_bam) | set { picard_out }
    // indexed_bam.view()

    // SNPs variant calling
    freebayes(fasta, fasta_index, indexed_bam) | set { vcf }
    build_config(fasta, gtf) | set { snpeff_config }
    snpeff(vcf, snpeff_config) | set { snpeff_output }
    bcftools_stats(vcf) | set { bcftools_out }
    annotated_vcfs = snpeff_output.map { id, vcf, html -> tuple(id, vcf) }
    qc_vcf = snpeff_output.map { id, vcf, html -> html }
    // qc_vcf.view()

    // SVs variant calling
    samtools_index(fasta) | set { fai }
    picard_dict(fasta) | set { dict }
    delly(indexed_bam, fasta, fai, dict) | set { bcf }
    convert_bcf_to_vcf(bcf) | set { sv_vcf }
    svviz(sv_vcf, indexed_bam, fasta, fai )

    // running multiqc on all files
    fastqc_out.mix(stats_out).mix(picard_out).mix(qc_vcf).mix(bcftools_out).collect() | multiqc
}

// fastqToVcf