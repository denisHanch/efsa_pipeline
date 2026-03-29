#!/usr/bin/env nextflow

// Processes for  short-read pipeline

/*
 * Trim adapters and low quality reads
 */
process trimgalore {

    publishDir "${params.out_dir}/${out_folder_name}/trimmed_reads", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    val out_folder_name

    output:
    tuple val(pair_id), path("*.fq.gz")

    script:
    """
    if [ \$(echo ${reads} | wc -w) -gt 1 ]; then
        trim_galore --gzip --cores ${task.cpus} -q 20 --illumina --phred33 --paired ${reads}
    else
        trim_galore --gzip --cores ${task.cpus} -q 20 --illumina --phred33 ${reads}
    fi
    """
}



process fastqc {

    publishDir "${params.out_dir}/${out_folder_name}/fastqc_out", mode: "copy"

    input:
    tuple val(pair_id), path(reads)
    val out_folder_name

    output:
    file("results_fastqc_${pair_id}")

    script:
    """
    mkdir results_fastqc_${pair_id}
    fastqc $reads -o results_fastqc_${pair_id}
    """
}


process multiqc {
    
    publishDir "${params.out_dir}/${out_folder_name}/multiqc", mode: "copy"

    input:
    path files
    val out_folder_name
    val filename

    output:
    file("*.html")

    script:
    """
    multiqc --filename ${filename} .
    """
}


// Processes for long-read pipeline

process nanoplot {

    publishDir "${params.out_dir}/${out_folder_name}/nanoplot", mode: "copy"

    input:
    tuple val(pair_id), path(reads)
    val out_folder_name
    
    output:
    path("${pair_id}_report")
    
    script:
    """
    NanoPlot --huge --fastq $reads --outdir ${pair_id}_report --threads ${task.cpus}
    """
}
