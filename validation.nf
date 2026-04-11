#!/usr/bin/env nextflow

include { pipeline } from "./main.nf"


process validate {

    tag "validate"
    container "ecomolegmo/validation:v1.0.4"
    publishDir "${params.out_dir}/data/valid/", mode: 'copy', overwrite: true

    input:
    path config_json

    output:
    path '*'

    script:
    """
    validation.sh --config ${config_json} 
    """
}

workflow {

    main:
    config_ch = Channel.fromPath("${projectDir}/data/inputs/config.json")
    validate(config_ch)
}