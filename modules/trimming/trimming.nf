#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process trimming {
    container 'biocontainers/trimmomatic:0.39--hdfd78af_2'
    tag { sample.baseName }

    input:
        path sample

    output:
        path "${sample.baseName}.${task.process}.fastq", emit: report
        stdout emit: log

    script:
    """
    # pretend trimming: just copy
    echo "${task.process} ${sample.baseName}"
    cp $sample ${sample.baseName}.${task.process}.fastq
    """
}