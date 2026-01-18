
// Processes for short-read pipeline

/*
 * Index fasta file with samtools
*/
process samtools_index {
    container 'staphb/samtools:latest'
    publishDir "${params.out_dir}/${out_folder_name}/samtools_index_dict", mode: 'copy'

    input:
    path fasta_file
    val out_folder_name

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
    publishDir "${params.out_dir}/${out_folder_name}/samtools_index_dict", mode: 'copy'
    
    input:
    path fasta_file
    val out_folder_name

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
    container 'dellytools/delly:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bam_index)
    each path(fasta_file)
    each path(fai)
    each path(dict)
    val out_folder_name

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
    container 'staphb/bcftools:latest'
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(bcf_file)
    val out_folder_name

    output:
    path("${pair_id}_sv_short_read.vcf")

    script:
    """
    bcftools view $bcf_file -Ov -o ${pair_id}_sv_short_read.vcf
    """

}


// Processes for long-read pipeline

/*
 * variant calling with cuteSV
*/
process cute_sv {
    container "${params.registry}/cutesv:v1.0.1"
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/cutesv_out", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index) 
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}_cutesv.vcf")
    

    script:
    """
    mkdir ${pair_id}_out
    cuteSV $bam_file $fasta_file ${pair_id}_cutesv.vcf ${pair_id}_out -t ${params.max_cpu}
    """
}


/*
 * variant calling with debreak
*/
process debreak {
    container "${params.registry}/debreak:v1.0.1"
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/debreak_out", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index)
    val out_folder_name

    output:
    tuple val(pair_id), path("debreak_out/${pair_id}_debreak.vcf")

    script:
    """
    debreak --bam $bam_file -r $fasta_file -o debreak_out -t ${params.max_cpu}
    mv debreak_out/debreak.vcf debreak_out/${pair_id}_debreak.vcf
    """
}


/*
 * variant calling with sniffles
*/
process sniffles {
    container "${params.registry}/sniffles:v1.0.1"
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/sniffles_out", mode: 'copy'

    input:
    each path(fasta_file)
    tuple val(pair_id), path(bam_file), path(bam_index) 
    val out_folder_name

    output:
    tuple val(pair_id), path("${pair_id}_sniffles.vcf")

    script:
    """
    sniffles --input $bam_file --vcf ${pair_id}_sniffles.vcf --reference $fasta_file --threads ${params.max_cpu}
    """
}


/*
 * merging SV
*/
process survivor {
    container "${params.registry}/survivor:v1.0.1"
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/survivor_out", mode: 'copy'

    input:
    tuple val(pair_id), path(sniffles_vcf)
    tuple val(pair_id), path(cute_vcf)
    tuple val(pair_id), path(debreak_vcf)
    val mapping_tag
    val out_folder_name


    output:
    tuple val(pair_id), path("${pair_id}_${mapping_tag}_sv_long_read.vcf")

    script:
    """
    # Create or overwrite the vcf_list.txt with the VCF file paths
    cp "${sniffles_vcf}" sniffles.vcf
    cp "${cute_vcf}"     cute.vcf
    cp "${debreak_vcf}"  debreak.vcf

    # Create VCF list file for SURVIVOR
    printf "%s\n" sniffles.vcf cute.vcf debreak.vcf > vcf_list.txt

    # Run the SURVIVOR merge command
    SURVIVOR merge vcf_list.txt 1000 1 1 1 0 30 ${pair_id}_${mapping_tag}_sv_long_read.vcf
    """
}


/*
 * Convert vcf file to tsv table 
*/
process vcf_to_table {

    container 'staphb/bcftools:1.23'
    publishDir "${params.out_dir}/tables", mode: 'copy'

    input:
    path vcf

    output:
    path "${vcf.simpleName}_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    output="${vcf.simpleName}_sv_summary.tsv"

    if [[ "${vcf.simpleName}" == "ref_x_modsyri" ]]; then
        {
            echo -e "chrom\tstart\tend\tsvtype"
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' "${vcf}"
        } > "\${output}"
    else
        {
            echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tdebreak_type\tsupporting_reads\tscore\tRDCN"
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%ALT\t[%RC]\t%QUAL\t[%RDCN]\n' "${vcf}"
        } >> "\${output}"
    fi
    """
}


process vcf_to_table_long {

    container 'staphb/bcftools:1.23'
    publishDir "${params.out_dir}/tables", mode: 'copy'

    input:
    path vcf

    output:
    path "${vcf.simpleName}_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    output="${vcf.simpleName}_sv_summary.tsv"

    echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsupporting_methods\tscore\tsupporting_reads" > "\${output}"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SUPP\t%QUAL\t[%DR{1}\t]\n' "${vcf}" | awk -F '\t' '{
    sum = 0;
    for(i=8; i<=NF; i++){
        if(\$i != "." && \$i != "") sum += (\$i + 0);
    }
    print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"sum;
    }' >> "\${output}"
    """
}