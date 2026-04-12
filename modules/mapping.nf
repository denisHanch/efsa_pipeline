#!/usr/bin/env nextflow

// Processes for short-read pipeline

process bwa_index {
    tag "$fasta_file"
    publishDir "${params.out_dir}/${out_folder_name}/bwa_index", mode: "copy"

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


process bwa_mapping {
    tag "$pair_id"

    input:
    path fasta_file
    path fasta_index
    tuple val(pair_id), path(reads)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.sam")

    script:
    """
    bwa mem ${fasta_file} ${reads} -t ${task.cpus} > ${pair_id}.sam
    """
}



process samtools_index_bam {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/bam", mode: "copy"

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


process picard {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/picard", mode: "copy"
    
    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

    output:
    path "${pair_id}_alignment_metrics.txt"

    script:
    """
    java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
    I=$bam_file \
    O=${pair_id}_alignment_metrics.txt \
    R=$fasta_file
    """
}


process samtools_stats {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/samtools_stats", mode: "copy"

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

process build_sv_flank_bed {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(sv_tsv)

    output:
    tuple val(pair_id), path("input_sv.tsv"), path("regions.bed")

    script:
    """
    set -euo pipefail

    cp "${sv_tsv}" input_sv.tsv

    awk '
      BEGIN { FS=OFS="\\t" }
      NR==1 {
        for (i=1; i<=NF; i++) {
          if (\$i=="chrom") c=i
          else if (\$i=="start") s=i
          else if (\$i=="end") e=i
        }
        next
      }
      {
        if (!c || !s || !e) next
        chrom=\$c; start=\$s+0; end=\$e+0
        if (end < start) { tmp=start; start=end; end=tmp }

        b0=start-101; if (b0 < 0) b0=0
        b1=start-1
        if (b1 > b0) print chrom, b0, b1, (NR-1)"|before"

        a0=end; if (a0 < 0) a0=0
        a1=end+100
        if (a1 > a0) print chrom, a0, a1, (NR-1)"|after"
      }
    ' input_sv.tsv > regions.bed
    """
}

process mosdepth {
    tag "$pair_id"
    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index), path(input_sv_tsv), path(regions_bed)
    val source_tag

    output:
    path "${pair_id}_${source_tag}_sv_summary.tsv"

    script:
    def outfile = "${pair_id}_${source_tag}_sv_summary.tsv"
    """
    set -euo pipefail

    coverage_tsv="raw_cov.tsv"
    prefix="${pair_id}_${source_tag}_flank"

    # Default to an empty coverage table when no regions are available.
    : > "\${coverage_tsv}"
    if [ -s "${regions_bed}" ]; then
      mosdepth --threads ${task.cpus} --no-per-base --by "${regions_bed}" "\${prefix}" "${bam_file}"
      gunzip -c "\${prefix}.regions.bed.gz" \
        | awk 'BEGIN{FS=OFS="\\t"} {split(\$4,a,/\\|/); print a[1], a[2], \$5}' > "\${coverage_tsv}"
    fi

    awk '
      BEGIN{FS=OFS="\\t"}
      NR==FNR { cov[\$1 "|" \$2] = \$3; next }
      FNR==1 { print \$0, "coverage_before_100bp", "coverage_after_100bp"; next }
      function get_cov(id, side, key, value) {
        key = id "|" side
        value = cov[key]
        return (value == "" ? "." : value)
      }
      {
        id = FNR - 1
        print \$0, get_cov(id, "before"), get_cov(id, "after")
      }
    ' "\${coverage_tsv}" "${input_sv_tsv}" > "${outfile}"
    """
}


// Processes for long-read pipeline

process minimap2 {
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
    minimap2 -t ${task.cpus} -ax $mapping_tag $fasta_file $reads > ${pair_id}.sam
    """
}


process samtools_sort {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/bam", mode: "copy"

    input:
    tuple val(pair_id), path(sam)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.bam")

    script:
    """
    samtools sort --threads ${task.cpus} $sam -o ${pair_id}.bam
    """
}

process calc_total_reads {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(bam), path(bam_index)

    output:
    env total

    script:
    """
    #!/usr/bin/env bash

    total=\$(samtools view -c "$bam")

    export total
    """
}


process calc_unmapped {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fastq)

    output:
    env reads 

    script:
    """
    #!/usr/bin/env bash
    if [[ "$fastq" == *.gz ]]; then
        total_lines=\$(zcat "$fastq" | wc -l)
    else
        total_lines=\$(wc -l < "$fastq")
    fi
    
    reads=\$((total_lines / 4))

    export reads
    """
}


process get_unmapped_reads {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/unmapped_fastq", mode: "copy"

    input:
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}_unmapped.fastq")

    script:
    """
    samtools view -@ ${task.cpus} -f 4 -b $bam_file | samtools fastq > ${pair_id}_unmapped.fastq
    """
}

process compare_unmapped {
    tag "$pair_id1"
    publishDir "${params.out_dir}/unmapped_stats", mode: "copy"

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
