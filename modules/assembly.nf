params.workflow_id = "fasta_ref_mod"

process nucmer {
    publishDir "${params.out_dir}/${params.workflow_id}/${prefix}", mode: "copy"

    input:
    tuple val(prefix), path(ref), path(mod)


    output:
    tuple val(prefix), path("${prefix}.delta")



    script:
    """
    nucmer --maxmatch -c 100 -b 500 -l 50 $ref $mod -p ${prefix}
    """

}

process delta_filter {
    publishDir "${params.out_dir}/${params.workflow_id}/${prefix}", mode: "copy"

    input:
    tuple val(prefix), path(delta)

    output:
    tuple val(prefix), path("${prefix}_filtered.delta")

    script:
    """
    delta-filter -m -i 90 -l 100 $delta > ${prefix}_filtered.delta
    """
}


process show_coords {
    publishDir "${params.out_dir}/${params.workflow_id}/${prefix}", mode: "copy"

    input:
    tuple val(prefix), path(filtered_delta)

    output:
    tuple val(prefix), path("${prefix}.filtered.coords")



    script:
    """
    show-coords -THrd $filtered_delta > ${prefix}.filtered.coords
    """

}

process syri {
    publishDir "${params.out_dir}/${params.workflow_id}", mode: "copy"


    input:
    tuple val(prefix), path(ref), path(mod), path(coords), path(filtered_delta)


    output:
    tuple val(prefix), path("${prefix}syri*.vcf")


    script:
    """
    syri -c $coords -d $filtered_delta  -r $ref -q $mod --prefix ${prefix} --nosnp
    """
}

process bgzip_tabix {
    publishDir "${params.out_dir}/${params.workflow_id}/${prefix}", mode: "copy"

    input:
    tuple val(prefix), path(vcf)

    output:
    tuple val(prefix), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi")

    script:
    """
    bgzip -c $vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz
    """
}


process bcftools_concat {
    publishDir "${params.out_dir}/${params.workflow_id}", mode: "copy"

    input:
    path(sv_vcf)
    path(tbis)

    output:
    tuple val('assembly'), path("assembly_concat.vcf")

    script:
    """
    bcftools concat -a -O v -o assembly_concat.vcf *.gz
    """
}