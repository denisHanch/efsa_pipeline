#!/usr/bin/env nextflow


// Include workflows
include { BASIC_ANALYSIS } from './workflows/basic.nf'
include { ADVANCED_ANALYSIS } from './workflows/advanced.nf'
include { fastqToVcf } from './workflows/short-read-ref.nf'


// Publishing process to copy final outputs
process publish_results {
    publishDir params.out_dir, mode: 'copy'
    
    input:
    path files
    
    output:
    path files
    
    script:
    """
    echo "Publishing files: ${files}"
    """
}

def debugView = { ch, name ->
   ch
    .map { line -> "[${name}]: $line" }
    .view()
}

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --workflow [basic|advanced] --input <input_file> --reference <reference_file>
    
    Options:
    --workflow     Choose workflow: 'basic' or 'advanced' (default: basic)
    --input        Input FASTQ file
    --reference    Reference FASTA file
    --outdir       Output directory (default: results)
    """.stripIndent()
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}

workflow {

    // 1) Load input FASTQ files & debug
    Channel.fromPath("${params.in_dir}/*.fastq")
           .set { samples_ch }
    debugView(samples_ch, 'RAW_SAMPLES')

    // Select workflow based on parameter
    if (params.workflow == 'basic') {
        log.info "Running BASIC analysis workflow"
        results = BASIC_ANALYSIS(samples_ch)

        // Combine all outputs into a single channel and publish once
        all_outputs = results.fastqc_out.mix(results.trimm_report)
        publish_results(all_outputs)
        
    } else if (params.workflow == 'advanced') {
        log.info "Running ADVANCED analysis workflow"
        results = ADVANCED_ANALYSIS(samples_ch)
        
        
        all_outputs = results.val_out.mix(results.mq_out)
        publish_results(all_outputs)
        
    } else {
        error "Invalid workflow: ${params.workflow}. Choose 'basic' or 'advanced'"
    }

    
}