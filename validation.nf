#!/usr/bin/env nextflow

process validate {

    tag "validate"
    container "ecomolegmo/validation:v1.0.4"
    publishDir "${projectDir}", mode: 'copy', overwrite: true

    input:
    path config_json

    output:
    path '*'

    script:
    """
    validate --config ${config_json} 
    """
}

workflow {

    main:
    config_ch = Channel.fromPath("${projectDir}/data/inputs/config.json")
    validate(config_ch)
}