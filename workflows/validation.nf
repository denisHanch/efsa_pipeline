#!/usr/bin/env nextflow

process validate {

    tag "validate"
    publishDir "${projectDir}/data/valid/", mode: 'copy', overwrite: true

    input:
    path config_json

    output:
    path 'validated_params.json', emit: params_json
    path 'run_*/**', emit: run_dir

    script:
    """
    validation.sh --config ${config_json} 
    """
}

workflow validation {

    take:
    config_json

    main:
    validate(config_json)

    emit:
    params_json = validate.out.params_json
}