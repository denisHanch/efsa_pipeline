#!/usr/bin/env nextflow

process multiqc {
    container 'multiqc/multiqc:v1.30'
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