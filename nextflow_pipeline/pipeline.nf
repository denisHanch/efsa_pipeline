#!/usr/bin/env nextflow

include { qc } from './qc.nf'


workflow {
    if (!params.fastqDir) {
            throw new IllegalArgumentException('--fastqDir must be provided')
        } else {

        println("Processing  files in directory: ${params.fastqDir}")
        Channel.fromPath("$params.fastqDir/*.fastq.gz") | set { fastqs }
        }
    fastqs
    | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq\.gz/)[0][1], it]} 
    | groupTuple(sort: true)
    | set { fastqs }
    
    fastqs.view()
}
