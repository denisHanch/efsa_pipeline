#!/usr/bin/env nextflow

// Processes for short-read pipeline

/*
 * Index the genome
*/
process bwa_index {
    container 'biocontainers/bwa:v0.7.17_cv1'
    tag "$fasta_file"
    publishDir "${params.out_dir}/${out_folder_name}/bwa_index", mode: 'copy'

    input:
    path fasta_file
    val out_folder_name

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

    input:
    each path(fasta_file)
    each path(fasta_index)
    tuple val(pair_id), path(reads)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.sam")

    script:
    """
    bwa mem ${fasta_file} ${reads} -t ${params.max_cpu} > ${pair_id}.sam
    """
}



/*
 * Indexing BAM file
*/
process samtool_index_bam {
    container 'staphb/samtools:1.23'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/bam", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)
    val out_folder_name

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
    publishDir "${params.out_dir}/${out_folder_name}/picard", mode: 'copy'
    
    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

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
    container 'staphb/samtools:1.23'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/samtools_stats", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file)
    val out_folder_name

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
    container 'staphb/minimap2:2.30'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(reads)
    each path(fasta_file)
    val mapping_tag
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.sam")
    

    script:
    """
    minimap2 -t ${params.max_cpu} -ax $mapping_tag $fasta_file $reads > ${pair_id}.sam
    """
}

/*
 * Sort reads with samtools
*/
process samtools_sort {
    container 'staphb/samtools:1.23'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/bam", mode: 'copy'

    input:
    tuple val(pair_id), path(sam)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.bam")

    script:
    """
    samtools sort --threads ${params.max_cpu} $sam -o ${pair_id}.bam
    """
}

/*
 * Calculate total reads
*/
process calc_total_reads {
    container 'staphb/samtools:1.23'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(bam), path(bam_index)

    output:
    env total

    script:
    """
    set -euo pipefail

    total=\$(samtools view -c "$bam" || echo 0)

    export total
    """
}


/*
 * Calculate unmapped reads
*/
process calc_unmapped {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fastq)

    output:
    env reads 

    script:
    """
    set -euo pipefail

    if [[ "$fastq" == *.gz ]]; then
        total_lines=\$(zcat "$fastq" | wc -l)
    else
        total_lines=\$(wc -l < "$fastq")
    fi
    
    reads=\$((total_lines / 4))

    num_files=\$(echo $fastq | wc -w)

    if [[ num_files -eq 2 ]]; then
        reads=\$((reads * 2))
    fi

    export reads
    """
}

/*
 * Export unmapped reads to fastq file
*/
process get_unmapped_reads {
    container 'staphb/samtools:1.23'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/unmapped_fastq", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}_unmapped.fastq")

    script:
    """
    samtools view -@ ${params.max_cpu} -f 4 -b $bam_file | samtools fastq > ${pair_id}_unmapped.fastq
    """
}

/*
 * Calculate unmapped statistics
*/
process compare_unmapped {
    tag "$pair_id1"
    publishDir "${params.out_dir}/unmapped_stats", mode: 'copy'

    input:
    tuple val(pair_id1), path(unmapped_1, stageAs: "unmapped_ref.fastq")
    tuple val(pair_id2), path(unmapped_2, stageAs: "unmapped_mod.fastq")
    val prefix

    output:
    path "${pair_id1}_${prefix}_read_stats.txt"

    script:
    """
    # Extract headers
    awk 'NR%4==1 {gsub(/^@/, "", \$0); print \$0}' unmapped_ref.fastq > ${pair_id1}_ref_headers.txt
    awk 'NR%4==1 {gsub(/^@/, "", \$0); print \$0}' unmapped_mod.fastq > ${pair_id2}_mod_headers.txt

    # Find intersection and unique reads
    comm -12 <(sort ${pair_id1}_ref_headers.txt) <(sort ${pair_id2}_mod_headers.txt) > intersect_ref_mod.txt
    comm -23 <(sort ${pair_id1}_ref_headers.txt) <(sort ${pair_id2}_mod_headers.txt) > unique_${pair_id1}_ref.txt
    comm -13 <(sort ${pair_id1}_ref_headers.txt) <(sort ${pair_id2}_mod_headers.txt) > unique_${pair_id2}_mod.txt

    # Count reads
    intersect_count=\$(wc -l < intersect_ref_mod.txt)
    unique_1_count=\$(wc -l < unique_${pair_id1}_ref.txt)
    unique_2_count=\$(wc -l < unique_${pair_id2}_mod.txt)

    # Output stats
    echo -e "Pair IDs: ${pair_id1} reference vs ${pair_id2} modified" > ${pair_id1}_${prefix}_read_stats.txt
    echo -e "Intersected Reads: \$intersect_count" >> ${pair_id1}_${prefix}_read_stats.txt
    echo -e "Unique Reads ${pair_id1} reference: \$unique_1_count" >> ${pair_id1}_${prefix}_read_stats.txt
    echo -e "Unique Reads ${pair_id2} modified: \$unique_2_count" >> ${pair_id1}_${prefix}_read_stats.txt
    """
}
