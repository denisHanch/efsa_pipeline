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
    def val_level_arg  = params.validation_level     ? "--validation-level ${params.validation_level}"  : ""
    def log_level_arg  = params.logging_level        ? "--logging-level ${params.logging_level}"         : ""
    def org_type_arg   = params.organism_type        ? "--type ${params.organism_type}"                  : ""
    def defrag_arg     = params.force_defragment_ref ? "--force-defragment-ref"                          : ""
    """
    validation.sh --config ${config_json} --threads ${params.max_cpu} ${val_level_arg} ${log_level_arg} ${org_type_arg} ${defrag_arg}
    """
}
