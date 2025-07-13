#!/usr/bin/env nextflow

process fastqc {
    container 'biocontainers/fastqc:0.11.9--0'
    tag { sample.baseName }
    
    input:
        path sample
    
    output:
        path "${sample.baseName}.${task.process}.txt", emit: report
        stdout emit: log
    
    script:
    """
    echo "report for ${sample.baseName}"
    echo "report for ${sample.baseName}" > ${sample.baseName}.${task.process}.txt
    """
}