process multiqc {
    container 'ewels/multiqc:1.13--py_0'
    tag { report.baseName }
    input:
        path report
    output:
        path "multiqc_report.txt"

    script:
    """
    echo "MultiQC aggregated for ${report.baseName}" > multiqc_report.txt
    """
}