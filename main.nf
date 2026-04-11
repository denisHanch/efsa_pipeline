#!/usr/bin/env nextflow

include { validation } from "./workflows/validation.nf"
include { pipeline } from "./workflows/workflow.nf"
include { logWorkflowCompletion } from "./modules/logs.nf"

// Help message
def helpMessage() {
    log.info"""
    Usage:

    nextflow run main.nf
    
    Options:

    -resume          Run pipeline from the point where it was interrupted or failed (Nextflow built-in)
    --out_dir        Output directory                                               (default: ${params.out_dir})
    --in_dir         Input directory                                                (default: ${params.in_dir})
    --max_cpu        Maximum CPUs per process                                       (default: ${params.max_cpu})
    --clean_work     Remove workdir after success                                   (default: ${params.clean_work})
    -with-report     Generate HTML execution report                                 (Nextflow built-in)
    -with-timeline   Produce timeline visualization                                 (Nextflow built-in)
    -with-dag        Produce DAG of workflow                                        (Nextflow built-in)
    --help           Show this help message
    """.stripIndent()
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}

workflow {
    config_ch = Channel.fromPath(params.config_json, checkIfExists: true)
    validation(config_ch)
    pipeline(validation.out.params_json)
}

logWorkflowCompletion("execution of main.nf")

workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    log.error "Check the process execution manifest in ${params.log_dir}/process_manifest.txt for details on which processes failed."
}