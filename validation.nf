#!/usr/bin/env nextflow

process validate {

    tag "validate"
    container "ecomolegmo/validation:v1.0.4@sha256:8e7560bcf3e6b831b385f9d8d7734bb2ac61040879a35a99f3c80fb1ab6edc16"
    errorStrategy 'terminate'
    memory '8 GB'
    time '1h'
    publishDir "${projectDir}/data/valid/", mode: 'copy', overwrite: true

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
    new File("data/outputs/tables/csv_per_summary").mkdirs()
    config_ch = Channel.fromPath("${projectDir}/data/inputs/config.json")
    validate(config_ch)
}