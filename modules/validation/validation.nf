#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.mode = params.mode ?: 'prod' // default to 'prod' if not specified

process validation {
    container 'pipeline-validation:dev'
    errorStrategy 'ignore' // Do not fail workflow on non-zero exit

    input:
      path config_file
      path validation_script

    output:
      path "validation.log", emit: log_file
      stdout emit: log

    script:
    """
    echo "Starting validation"
    python $validation_script $config_file > validation.log 2>&1
    echo \$? > exit_code.txt
    """
}

workflow {
    if (params.mode == 'test') {
        // Test workflow: run validation for all test configs
        def configs_ch = Channel.fromPath('./tests/*/config.json')
        def script_ch = Channel.value('/EFSA_workspace/modules/validation/main.py')
        validation(configs_ch, script_ch)
    } else {
        // Prod workflow: run validation for a single config (example)
        def config_ch = Channel.fromPath('/EFSA_workspace/data/inputs/config.json')
        def script_ch = Channel.value('/EFSA_workspace/modules/validation/main.py')
        validation(config_ch, script_ch)
        
    }
}