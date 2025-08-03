#!/usr/bin/env nextflow

/*
 * Index the genome
*/
process bwa_index {
    container 'simonovaekat/bwa-samtools:latest'
    tag "$fasta_file"
    publishDir "${params.out_dir}/bwa_index", mode: 'copy'

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
    container 'simonovaekat/bwa-samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/bam", mode: 'copy'

    input:
    path fasta_file
    path fasta_index
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.bam")

    script:
    """
    bwa mem ${fasta_file} ${reads} | \
        samtools view -bS - | \
        samtools sort -o ${pair_id}.bam
    """
}




/*
 * Indexing BAM file
*/
process samtool_index_bam {
    container 'simonovaekat/bwa-samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/bam", mode: 'copy'

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
    publishDir "${params.out_dir}/picard", mode: 'copy'
    
    input:
    path fasta_file
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
    container 'simonovaekat/bwa-samtools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/samtools_stats", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)

    output:
    path "*.stats"

    script:
    """
    samtools stats $bam_file > ${bam_file}.stats
    """
}



// Long read function pipeline


process minimap2 {
container 'staphb/minimap2:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/minimap2", mode: 'copy'

    input:

    
    output:
    

    script:
    """
    
    """
}