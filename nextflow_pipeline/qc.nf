#!/usr/bin/env nextflow

/*
 * Step 1. Trim adapters and low quality reads
 */
process trimgalore {

    container 'vibsinglecellnf/trimgalore:trimgalore-0.6.7-cutadapt-4.1'
    tag "$pair_id"
    publishDir "${params.outDir}/trimmed_reads", mode: 'copy'
    cpus 4
    memory '4 GB'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), file('out/*_val_1.fq.gz'), file('out/*_val_2.fq.gz')

    script:
    """
    trim_galore --gzip --cores ${task.cpus} -q 20 --illumina --phred33 --paired -o out ${reads[0]} ${reads[1]}
    """
}



/*
 * Step 2. Check quality of reads
 */
process fastqc {
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "FASTQC on $pair_id"
    publishDir "$params.outDir/fastqc_out", mode: "link"
    cpus 4
    memory '4 GB'

    input:
    tuple val(pair_id), path(reads)

    output:
    file("results_fastqc_${pair_id}")

    script:
    """
    mkdir results_fastqc_${pair_id}
    fastqc $reads -o results_fastqc_${pair_id}
    """
}

/*
 * Step 3. Generate multiqc report
 */
process multiqc {
    container 'staphb/multiqc'
    publishDir "$params.outDir/multiqc", mode: "copy"
    cpus 2
    memory 4.GB

    input:
    file('fastq/*')

    output:
    file('multiqc_report_primary_qc.html')

    script:
    """
    multiqc fastq/
    mv multiqc_report.html multiqc_report_primary_qc.html
    """
}

workflow qc {
    take:
        fastqs
    main:
        fastqs | trimgalore | set { trimmed }

        fastqc(fastqs).collect() | multiqc
    emit:
        fastqs

} 