
// Processes for short-read pipeline

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

// Processes for long-read pipeline


/*
 * variant calling with cuteSV
*/
process cute_sv {
container 'simonovaekat/cutesv:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/cutesv_out", mode: 'copy'

    input:
    path fasta_file
    tuple val(pair_id), path(bam_file), path(bam_index) 

    output:
    tuple val(pair_id), path("${pair_id}_cutesv.vcf")
    

    script:
    """
    mkdir ${pair_id}_out
    cuteSV $bam_file $fasta_file ${pair_id}_cutesv.vcf ${pair_id}_out
    """
}



/*
 * variant calling with debreak
*/
process debreak {
container 'simonovaekat/debreak:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/debreak_out", mode: 'copy'
    cpus 4

    input:
    path fasta_file
    tuple val(pair_id), path(bam_file), path(bam_index) 

    output:
    tuple val(pair_id), path("${pair_id}_debreak.vcf")

    script:
    """
    debreak --bam $bam_file -r $fasta_file -o ${pair_id}_debreak.vcf -t ${task.cpus}
    """
}


/*
 * variant calling with sniffles
*/
process sniffles {
container 'simonovaekat/sniffles:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/sniffles_out", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(bam_file), path(bam_index) 

    output:
    tuple val(pair_id), path("${pair_id}_sniffles.vcf")

    script:
    """
    sniffles --input $bam_file --vcf ${pair_id}_sniffles.vcf --threads ${task.cpus}
    """
}


/*
 * merging SV
*/
process survivor {
container 'simonovaekat/survivor:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/survivor_out", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(sniffles_vcf)
    tuple val(pair_id), path(cute_vcf)
    tuple val(pair_id), path(debreak_vcf)


    output:
    path "${pair_id}_merged.vcf"

    script:
    """
    # Create or overwrite the vcf_list.txt with the VCF file paths
    echo $sniffles_vcf > vcf_list.txt
    echo $cute_vcf >> vcf_list.txt
    echo $debreak_vcf >> vcf_list.txt

    # Run the SURVIVOR merge command
    SURVIVOR merge vcf_list.txt 1000 1 1 1 0 30 ${pair_id}_merged.vcf
    """
}