process trimming {
    container 'biocontainers/trimmomatic:0.39--hdfd78af_2'
    tag { sample.baseName }
    input:
        path sample
    output:
        path "${sample.baseName}.trimmed.fastq"

    script:
    """
    # pretend trimming: just copy
    cp $sample ${sample.baseName}.trimmed.fastq
    """
}