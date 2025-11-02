
// Processes for short-read pipeline

/*
 * Calling short variants with FreeBayes
*/
process freebayes {
    container 'biocontainers/freebayes:v1.2.0-2-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/vcf", mode: 'copy'

    input:
    each path(fasta_file)
    each path(fasta_index)
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

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
    publishDir "${params.out_dir}/short-ref/vcf", mode: 'copy'


    input:
    each path(fasta_file)
    each path(feature_file)
    val feature_tag
    val build_setting

    output:
    tuple path("genome_id.txt"), path("snpEff.config")

    script:
    """
    fasta_path='${fasta_file}'
    feature_path='${feature_file}'

    # Extract genome ID from first header line
    genome_id=\$(grep '^>' \$fasta_path | head -n 1 | cut -d ' ' -f1 | sed 's/^>//' | tr -cd '[:alnum:]_')
    echo "Using genome ID: \$genome_id"

    mkdir -p data/\$genome_id
    cp \$fasta_path data/\$genome_id/sequences.fa
    cp \$feature_path data/\$genome_id/genes.${feature_tag}

    echo "\$genome_id.genome : Custom Genome" > snpEff.config

    snpEff build -${build_setting} -v \$genome_id -c snpEff.config

    echo \$genome_id > genome_id.txt
    """
}


/*
 * Annotating variants with SNPeff
*/
process snpeff {
    container 'biocontainers/snpeff:v4.1k_cv3'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)
    each path(genome_id_file)
    each path(snpeff_config)


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

// Processes for comparision of VCF files (in main.nf pipeline)

/*
 * Truvari requires sorted vcfs for comparision
*/
process sortVcf {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/truvari", mode: 'copy'

    input:
    tuple val(pair_id),  path(vcf_file)

    output:
    tuple val(pair_id),  path("${pair_id}.vcf.gz")

    script:
    """
    bcftools sort $vcf_file -Oz -o ${pair_id}.vcf.gz
    """
}

/*
 * Truvari require indexed vcfs for comparision
*/
process indexVcf {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/truvari", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)

    output:
    tuple val(pair_id), path(vcf_file), path("${vcf_file}.csi")


    script:
    """
    bcftools index $vcf_file
    """
}

/*
 * Running comparision btw files
*/
process truvari {
    container "${params.registry}/truvari:latest"
    tag "$pair_id1 & $pair_id2"
    publishDir "${params.out_dir}/truvari", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id1), path(vcf1), path(index1), val(pair_id2), path(vcf2), path(index2)

    output:
    path("${pair_id1}_${pair_id2}_truvari")


    script:
    """
    truvari bench -b $vcf1 -c $vcf2 -f $fasta_file -o ${pair_id1}_${pair_id2}_truvari --passonly
    """
}