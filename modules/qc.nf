#!/usr/bin/env nextflow

// Processes for  short-read pipeline

/*
 * Trim adapters and low quality reads
 */
process trimgalore {

    container 'vibsinglecellnf/trimgalore:trimgalore-0.6.7-cutadapt-4.1'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/trimmed_reads", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path('out/*.fq.gz')

    script:
    """
    trim_galore --gzip --cores ${task.cpus} -q 20 --illumina --phred33 --paired -o out ${reads}
    """
}



/*
 * Check quality of reads
 */
process fastqc {
    container 'biocontainers/fastqc:v0.11.9_cv8'
    tag "FASTQC on $pair_id"
    publishDir "$params.out_dir/short-ref/fastqc_out", mode: 'copy'


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
 * Generate multiqc report
 */
process multiqc {
    container 'staphb/multiqc'
    publishDir "$params.out_dir/short-ref/multiqc", mode: "copy"

    input:
    file('fastq/*')

    output:
    file('multiqc_report.html')

    script:
    """
    multiqc .
    """
}


// Processes for long-read pipeline

/*
 * Check QC of input samples
*/
process nanoplot {
container 'staphb/nanoplot:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/nanoplot", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    
    output:
    path('nanoplot_report')
    
    script:
    """
    NanoPlot --fastq $reads --outdir nanoplot_report
    """
}