/*
 * Calling short variants with FreeBayes
*/
process freebayes {
    container 'biocontainers/freebayes:v1.2.0-2-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/vcf", mode: 'copy'

    input:
    path fasta_file
    path fasta_index
    tuple val(pair_id), path(bam_file), path(bam_index)

    output:
    tuple val(pair_id), path("${pair_id}.vcf")

    script:
    """
    freebayes -f $fasta_file --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 --min-alternate-count 3 $bam_file > ${pair_id}.vcf
    """
}

/*
* SNPEff build
*/

process build_config {
    container 'biocontainers/snpeff:v4.1k_cv3'
    publishDir "${params.out_dir}/vcf", mode: 'copy'


    input:
    path fasta_file
    path gtf_file

    output:
    tuple path("genome_id.txt"), path("snpEff.config")

    script:
    """
    fasta_path='${fasta_file}'
    gtf_path='${gtf_file}'

    # Extract genome ID from first header line
    genome_id=\$(grep '^>' \$fasta_path | head -n 1 | cut -d ' ' -f1 | sed 's/^>//' | tr -cd '[:alnum:]_')
    echo "Using genome ID: \$genome_id"

    mkdir -p data/\$genome_id
    cp \$fasta_path data/\$genome_id/sequences.fa
    cp \$gtf_path data/\$genome_id/genes.gtf

    echo "\$genome_id.genome : Custom Genome" > snpEff.config

    snpEff build -gtf22 -v \$genome_id -c snpEff.config

    echo \$genome_id > genome_id.txt
    """
}


/*
 * Annotating variants with SNPeff
*/
process snpeff {
    container 'biocontainers/snpeff:v4.1k_cv3'
    tag "$pair_id"
    publishDir "${params.out_dir}/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)
    tuple path(genome_id_file), path(snpeff_config)


    output:
    tuple val(pair_id), path("${pair_id}_annotated.vcf"), path("${pair_id}_summary.csv")


    script:
    """
    genome_id=\$(cat ${genome_id_file})

    snpEff -v -c ${snpeff_config} \$genome_id ${vcf_file} > ${pair_id}_annotated.vcf -csvStats "${pair_id}_summary.csv"
    """
}


/*
 * Getting variants statistics with bcftools stats
*/
process bcftools_stats {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/bcftools_stats", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)

    output:
    path "${pair_id}.bcftools_stats.txt"

    script:
    """
    bcftools stats ${vcf_file} > ${pair_id}.bcftools_stats.txt
    """
}
