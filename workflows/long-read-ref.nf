#!/usr/bin/env nextflow

include { nanoplot } from '../modules/qc.nf'
include { samtool_stats; minimap2; samtools_sort; samtool_index_bam } from '../modules/mapping.nf'
include { sniffles; debreak; cute_sv; survivor } from '../modules/sv_calling.nf'
include { bcftools_stats } from '../modules/variant_calling.nf'


workflow {
    // Processing inputs
    Channel.fromPath("$params.fastqDir/ref.fa") | set { fasta }

    if (!params.fastqDir) {
            throw new IllegalArgumentException('--fastqDir must be provided')
        } else {

        println("Processing  files in directory: ${params.fastqDir}")
        Channel.fromPath("$params.fastqDir/*_subreads.fastq.gz")
        .map { file -> 
        def name = file.baseName.replaceFirst('.fastq', '')
        return [name, file]
    }
    .set { fastqs }

    fastqs.view()
     }
    // QC
    fastqs | nanoplot

    // mapping
    minimap2(fastqs, fasta) | samtools_sort | samtool_index_bam | set { indexed_bam }

    // variant calling
    cute_sv(fasta, indexed_bam) | set { cute_vcf }
    debreak(fasta, indexed_bam) | set { debreak_vcf }
    sniffles(indexed_bam) | set { sniffles_vcf }

    survivor(cute_vcf, debreak_vcf, sniffles_vcf) | set { merged_vcf }
    bcftools_stats(merged_vcf) | multiqc
}