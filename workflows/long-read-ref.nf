#!/usr/bin/env nextflow

// include { convert_bcf_to_vcf } from '../modules/sv_calling.nf'
include { nanoplot } from '../modules/qc.nf'
// include { samtool_stats } from '../modules/mapping.nf'

workflow {
    // Processing inputs
    if (!params.fastqDir) {
            throw new IllegalArgumentException('--fastqDir must be provided')
        } else {

        println("Processing  files in directory: ${params.fastqDir}")
        Channel.fromPath("$params.fastqDir/*_subreads.fastq.gz") | set { fastqs }
        }
        fastqs | nanoplot
}