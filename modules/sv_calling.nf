
// Processes for short-read pipeline


process samtools_index {
    
    publishDir "${params.out_dir}/${out_folder_name}/samtools_index_dict", mode: 'copy'

    input:
    path fasta_file
    val out_folder_name

    output:
    path("${fasta_file}.fai")

    script:
    """
    samtools faidx "$fasta_file"
    """   
}


process picard_dict {
    
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


process delly {

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


process convert_bcf_to_vcf {

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

process cute_sv {

    publishDir "${params.out_dir}/${out_folder_name}/cutesv_out", mode: "copy"

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



process debreak {

    publishDir "${params.out_dir}/${out_folder_name}/debreak_out", mode: "copy"

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


process sniffles {

    publishDir "${params.out_dir}/${out_folder_name}/sniffles_out", mode: "copy"

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


process survivor {

    publishDir "${params.out_dir}/${out_folder_name}/survivor_out", mode: "copy"

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



process vcf_to_table_asm {

    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'

    input:
    tuple val(name), path(vcf)

    output:
    path "*_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    echo -e "chrom\tstart\tend\tsvtype" > "${name}_sv_summary.tsv"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' "${vcf}" >> "${name}_sv_summary.tsv"

    """
}

process vcf_to_table_short {

    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'

    input:
    tuple val(name), path(vcf)

    output:
    path "*_sv_summary.tsv"

    script:
    """
    set -euxo pipefail

    echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsvlen\tsupporting_reads\tscore\tRDCN" > "${name}_short_sv_summary.tsv"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%ALT\t%INFO/SVLEN\t%QUAL\t[%RDCN]\n' "${vcf}"  >> "${name}_short_sv_summary.tsv"
    """
}

process vcf_to_table_long {

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

    echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsvlen\tsupporting_methods\tscore\tsupporting_reads" > "\${output}"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%QUAL\t[%DR{1}\t]\n' "${vcf}" | awk -F '\t' '{
    sum = 0;
    for(i=9; i<=NF; i++){
        if(\$i != "." && \$i != "") sum += (\$i + 0);
    }
    print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8"\t"sum;
    }' >> "\${output}"
    """
}

process restructure_sv_tbl {

    publishDir "${params.out_dir}/tables", mode: 'copy'

    input:
    path script
    tuple path(assembly_tsv), path(long_ont_tsv), path(long_pb_tsv), path(short_tsv)

    output:
    path "csv_per_sv_summary"

    script:
    """
    python ${script} --asm ${assembly_tsv} --short ${short_tsv} --long_ont ${long_ont_tsv} --long_pacbio ${long_pb_tsv} --out csv_per_sv_summary
    """
}


process create_empty_tbl {
    
    publishDir "${params.out_dir}/tables/tsv", mode: 'copy'
    
    input:
    val prefix

    output:
    path "empty_${prefix}_sv_summary.tsv"

    script:
    """
    echo -e "chrom\tstart\tend\tsvtype\tinfo_svtype\tsupporting_methods\tscore\tsupporting_reads" > empty_${prefix}_sv_summary.tsv
    """
}