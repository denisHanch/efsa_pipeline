
// Processes for short-read pipeline


process freebayes {
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/vcf", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}.vcf")

    script:
    """
    freebayes -f $fasta_file --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 --min-alternate-count 3 $bam_file > ${pair_id}.vcf
    """
}

process bcftools_stats {

    publishDir "${params.out_dir}/${out_folder_name}/bcftools_stats", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)
    val out_folder_name

    output:
    path "${pair_id}.bcftools_stats.txt"

    script:
    """
    bcftools stats ${vcf_file} > ${pair_id}.bcftools_stats.txt
    """
}
