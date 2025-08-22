#!/usr/bin/env nextflow

// Processes for short-read pipeline

/*
 * Index the genome
*/
process bwa_index {
    container 'biocontainers/bwa:v0.7.17_cv1'
    tag "$fasta_file"
    publishDir "${params.out_dir}/short-ref/bwa_index", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file}.*" 

    script:
    """
    bwa index ${fasta_file}
    """
}



/*
 * Mapping reads to the genome
*/
process bwa_mapping {
    container 'biocontainers/bwa:v0.7.17_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/bam", mode: 'copy'

    input:
    each path(fasta_file)
    each path(fasta_index)
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.sam")

    script:
    """
    bwa mem ${fasta_file} ${reads} > ${pair_id}.sam
    """
}



/*
 * Indexing BAM file
*/
process samtool_index_bam {
    container 'staphb/samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/bam", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)

    output:
    tuple val(pair_id), path(bam_file), path("${bam_file}.bai")

    script:
    """
    samtools index $bam_file
    """
}


/*
 * Running picard to get mapping statistics
*/
process picard {
    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/picard", mode: 'copy'
    
    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index)


    output:
    path "${pair_id}_alignment_metrics.txt"

    script:
    """
    picard CollectAlignmentSummaryMetrics \
    I=$bam_file \
    O=${pair_id}_alignment_metrics.txt \
    R=$fasta_file
    """
}



/*
 * Running samtools stats to get mapping statistics
*/
process samtool_stats {
    container 'staphb/samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/samtools_stats", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)

    output:
    path "*.stats"

    script:
    """
    samtools stats $bam_file > ${bam_file}.stats
    """
}


// Processes for long-read pipeline

/*
 * Map long reads
*/
process minimap2 {
container 'staphb/minimap2:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/minimap2", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)
    each path(fasta_file)

    output:
    tuple val(pair_id), path("${pair_id}.sam")
    

    script:
    """
    minimap2 -ax map-pb $fasta_file $reads > ${pair_id}.sam
    """
}

/*
 * Sort reads with samtools
*/
process samtools_sort {
container 'staphb/samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/samtools", mode: 'copy'

    input:
    tuple val(pair_id), path(sam)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.bam")

    script:
    """
    samtools sort $sam -o ${pair_id}.bam
    """
}