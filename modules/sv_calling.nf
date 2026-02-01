
// Processes for short-read pipeline

/*
 * Index fasta file with samtools
*/
process samtools_index {
    container params.containers.samtools
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
    container params.containers.picard
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
    container params.containers.delly
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
    container params.containers.bcftools
    tag "$pair_id"
    publishDir "${params.out_dir}/${out_folder_name}/vcf", mode: 'copy'

    input:
    tuple val(pair_id), path(bcf_file)
    val out_folder_name

    output:
    tuple val(pair_id),path("${pair_id}_sv_short_read.vcf")

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
    container params.containers.cutesv
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
    container params.containers.debreak
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
    container params.containers.sniffles
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
    container params.containers.survivor
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

    container params.containers.bcftools
    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'

    input:
    tuple val(name), path(vcf)

    output:
    path "*_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    if [[ "${name}" == "assembly" ]]; then
        {
            echo -e "chrom\tstart\tend\tsvtype"
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' "${vcf}"
        } > "${name}_sv_summary.tsv"
    else
        {
            echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsupporting_reads\tscore\tRDCN"
            bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%ALT\t[%RC]\t%QUAL\t[%RDCN]\n' "${vcf}"
        } >> "${name}_short_sv_summary.tsv"
    fi
    """
}


process vcf_to_table_long {

    container params.containers.bcftools
    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'

    input:
    val tag
    tuple val(pair_id), path(vcf)

    output:
    path "${pair_id}_${tag}_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    output="${pair_id}_${tag}_sv_summary.tsv"

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

process restructure_sv_table {

    publishDir "${params.out_dir}/tables", mode: 'copy'

    input:
    path script
    tuple path(assembly_tsv), path(long_ont_tsv), path(long_pb_tsv), path(short_tsv)

    output:
    path "csv_per_sv_sumary"

    script:
    """
    python ${script} --asm ${assembly_tsv} --short ${short_tsv} --long_ont ${long_ont_tsv} --long_pacbio ${long_pb_tsv} --out csv_per_sv_sumary
    """
}


process create_empty_tbl {
    
    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'
    
    input:
    val prefix

    output:
    path "empty_${prefix}_summary.tsv"

    script:
    """
    echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsupporting_methods\tscore\tsupporting_reads" > empty_${prefix}_summary.tsv
    """
}