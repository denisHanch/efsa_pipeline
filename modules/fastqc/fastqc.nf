process fastqc {
    container 'biocontainers/fastqc:0.11.9--0'
    tag { sample.baseName }
    input:
        path sample
    output:
        path "${sample.baseName}.fastqc.txt"

    script:
    """
    echo "FASTQC report for ${sample.baseName}" > ${sample.baseName}.fastqc.txt
    """
}