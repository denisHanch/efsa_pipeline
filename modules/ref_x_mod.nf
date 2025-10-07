params.workflow_id = "fasta_ref_mod"

process nucmer {
    container "biocontainers/mummer:v3.23dfsg-4-deb_cv1"
    containerOptions = "--user root"
    publishDir "${params.out_dir}/${params.workflow_id}", mode: 'copy'

    input:
    tuple val(prefix), path(ref), path(mod)


    output:
    path("${prefix}.delta")



    script:
    """
    nucmer --maxmatch -c 100 -b 500 -l 50 $ref $mod -p ${prefix}
    """

}

process deltaFilter {
    container "biocontainers/mummer:v3.23dfsg-4-deb_cv1"
    publishDir "${params.out_dir}/${params.workflow_id}", mode: 'copy'

    input:
    val prefix 
    path delta

    output:
    path("${prefix}_filtered.delta")

    script:
    """
    delta-filter -m -i 90 -l 100 $delta > ${prefix}_filtered.delta
    """
}


process showCoords {
    container "biocontainers/mummer:v3.23dfsg-4-deb_cv1"
    publishDir "${params.out_dir}/${params.workflow_id}", mode: 'copy'

    input:
    val prefix
    path filtered_delta

    output:
    path("${prefix}.filtered.coords")



    script:
    """
    show-coords -THrd $filtered_delta > ${prefix}.filtered.coords
    """

}

process syri {
    container "${params.registry}/syri/:latest"
    publishDir "${params.out_dir}/${params.workflow_id}", mode: 'copy'
    publishDir ( params.map_to_mod_fa ? "${params.out_dir}/${params.workflow_id}" : "${params.out_dir}/final_vcf" ), mode: 'copy', saveAs: { filename -> filename.endsWith(".vcf") ? filename : null }


    input:
    tuple val(prefix), path(ref), path(mod)
    path coords
    path filtered_delta


    output:
    path("${prefix}syri*")


    script:
    """
    syri -c $coords -d $filtered_delta  -r $ref -q $mod --prefix ${prefix} --nosnp
    """
}