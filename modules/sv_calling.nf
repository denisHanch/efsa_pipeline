
process samtools_index {
    container 'simonovaekat/bwa-samtools:latest'
    publishDir "${params.out_dir}/samtools_index_dict", mode: 'copy'

    input:
    path fasta_file

    output:
    path("${fasta_file}.fai")

    script:
    """
    samtools faidx $fasta_file
    """   
}

process picard_dict {
    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    publishDir "${params.out_dir}/samtools_index_dict", mode: 'copy'
    
    input:
    path fasta_file

    output:
    path("*.dict")

    script:
    """
    dict_name=\$(basename $fasta_file | sed 's/\\.[^.]*\$/.dict/')
    picard CreateSequenceDictionary R=$fasta_file O=\$dict_name
    """  
}


process delly {
    container 'biocontainers/delly:v0.8.1-2-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index)
    path fasta_file
    path fai
    path dict

    output:
    tuple val(pair_id),  path("${pair_id}_sv.bcf")
    

    script:
    """
    delly call -g $fasta_file -o ${pair_id}_sv.bcf $bam_file
    """
}

process convert_bcf_to_vcf {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/vcf", mode: 'copy'

    input:
    tuple val(pair_id),  path(bcf_file)

    output:
    tuple val(pair_id),  path("*.vcf")

    script:
    """
    bcftools view $bcf_file -Ov -o ${pair_id}_sv.vcf
    """

}


process svviz {
    container 'simonovaekat/svviz2:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/svviz", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    path fasta_file
    path fai


    output:
    path "svviz2_output/index.html"

    script:
    """
    svviz2 --ref $fasta_file --variants $vcf_file $bam_file
    """
}