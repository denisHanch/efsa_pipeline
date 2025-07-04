#!/usr/bin/env nextflow

include { header } from './util'

/*
 * Step 1. Trim adapters and low quality reads
 */
process trimgalore {
    container 'vibsinglecellnf/trimgalore:0.6.7'
    tag "$pair_id"
    publishDir "$params.outDir/trimmed_reads"
    cpus 4
    memory 8.GB

    input:
    tuple val(pair_id), file(r1), file(r2)

    output:
    tuple val(pair_id), file('out/*val_1.fq.gz'), file('out/*val_2.fq.gz')

    script:
    if (params.log) {
        log.info """${header('T R I M - G A L O R E')}
        reads : '${reads}'
        """
    }

    """
    trim_galore --gzip --cores $task.cpus -q 20 --illumina --phred33 --paired -o out $reads

    ls -l *
    """
}


/*
 * Step 2. Check quality of reads
 */
process fastqc {
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "FASTQC on $pair_id"
    publishDir "$params.outDir/fastqc_out"
    cpus 4
    memory 8.GB

    input:
    tuple val(pair_id), file(r1), file(r2)

    output:
    file("results_fastqc_${pair_id}")

    script:
    if (params.log) {
        log.info """${header('F A S T Q C')}
        reads :  '${reads}'
        """
    }

    """
    mkdir results_fastqc_${pair_id}
    fastqc --threads 4 -f fastq -q $reads -o results_fastqc_${pair_id}
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
    if (params.log) {
         log.info """${header(' M U L T I Q C  F A S T Q')}"""
    }
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

        fastqc(trimmed).collect() | multiqc
    emit:
        trimmed

} 