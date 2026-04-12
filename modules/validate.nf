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
    validation.sh --config ${config_json} --threads ${params.max_cpu} 
    """
}
