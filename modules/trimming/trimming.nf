#!/usr/bin/env nextflow

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
    echo "SAMPLE TRIMMING DATA" >> ${sample.baseName}.${task.process}.fastq
    """
}