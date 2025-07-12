#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process multiqc {
    container 'ewels/multiqc:1.13--py_0'
    tag { report.baseName }

    input:
        path report
    output:
        path "multiqc_report.${task.process}.txt", emit: report
        stdout emit: log

    script:
    """
    echo "MultiQC aggregated for ${report.baseName}" > multiqc_report.${task.process}.txt
    """
}