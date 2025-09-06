
// Processes for short-read pipeline

/*
 * Index fasta file with samtools
*/
process samtools_index {
    container 'staphb/samtools:latest'
    publishDir "${params.out_dir}/short-ref/samtools_index_dict", mode: 'copy'

    input:
    path fasta_file

    output:
    path("${fasta_file}.fai")

    script:
    """
    samtools faidx $fasta_file
    """   
}

/*
 * Create picard dictionary to run delly
*/
process picard_dict {
    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    publishDir "${params.out_dir}/short-ref/samtools_index_dict", mode: 'copy'
    
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

/*
 * Running delly SV caller
*/
process delly {
    container 'biocontainers/delly:v0.8.1-2-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index)
    each path(fasta_file)
    each path(fai)
    each path(dict)

    output:
    tuple val(pair_id),  path("${pair_id}_sv.bcf")
    

    script:
    """
    delly call -g $fasta_file -o ${pair_id}_sv.bcf $bam_file
    """
}

/*
 * Convert bcf (default delly output) to vcf
*/
process convert_bcf_to_vcf {
    container 'biocontainers/bcftools:v1.9-1-deb_cv1'
    tag "$pair_id"
    publishDir "${params.out_dir}/short-ref/vcf", mode: 'copy'
    publishDir "${params.out_dir}/final_vcf", mode: 'copy'

    input:
    tuple val(pair_id),  path(bcf_file)

    output:
    tuple val(pair_id),  path("${pair_id}_sv_short_read.vcf")

    script:
    """
    bcftools view $bcf_file -Ov -o ${pair_id}_sv_short_read.vcf
    """

}

// not working currently - trying to find nice tool for visualization of SVs
process svviz {
    container "${params.registry}/svviz2:latest"
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/svviz", mode: 'copy'

    input:
    tuple val(pair_id), path(vcf_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    each path(fasta_file)
    each path(fai)


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
    container "${params.registry}/cutesv:latest"
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/cutesv_out", mode: 'copy'

    input:
    each path(fasta_file)
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
    container "${params.registry}/debreak:latest"
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/debreak_out", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index) 

    output:
    tuple val(pair_id), path("debreak_out/${pair_id}_debreak.vcf")

    script:
    """
    debreak --bam $bam_file -r $fasta_file -o debreak_out
    mv debreak_out/debreak.vcf debreak_out/${pair_id}_debreak.vcf
    """
}


/*
 * variant calling with sniffles
*/
process sniffles {
    container "${params.registry}/sniffles:latest"
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/sniffles_out", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index) 

    output:
    tuple val(pair_id), path("${pair_id}_sniffles.vcf")

    script:
    """
    sniffles --input $bam_file --vcf ${pair_id}_sniffles.vcf
    """
}


/*
 * merging SV
*/
process survivor {
    container "${params.registry}/survivor:latest"
    tag "$pair_id"
    publishDir "${params.out_dir}/long-ref/survivor_out", mode: 'copy'
    publishDir "${params.out_dir}/final_vcf", mode: 'copy'


    input:
    tuple val(pair_id), path(sniffles_vcf)
    tuple val(pair_id), path(cute_vcf)
    tuple val(pair_id), path(debreak_vcf)


    output:
    tuple val(pair_id), path("${pair_id}_sv_long_read.vcf")

    script:
    """
    # Create or overwrite the vcf_list.txt with the VCF file paths
    cp "${sniffles_vcf}" sniffles.vcf
    cp "${cute_vcf}"     cute.vcf
    cp "${debreak_vcf}"  debreak.vcf

    # Create VCF list file for SURVIVOR
    printf "%s\n" sniffles.vcf cute.vcf debreak.vcf > vcf_list.txt

    # Run the SURVIVOR merge command
    SURVIVOR merge vcf_list.txt 1000 1 1 1 0 30 ${pair_id}_sv_long_read.vcf
    """
}