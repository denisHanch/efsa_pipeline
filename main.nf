#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { fastqc }    from './modules/fastqc/fastqc.nf'
include { trimming }  from './modules/trimming/trimming.nf'
include { multiqc }   from './modules/multiqc/multiqc.nf'
include { validation }  from './modules/validation/validation.nf'

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

workflow {

    // 1) Load input FASTQ files & debug
    Channel.fromPath("${params.in_dir}/*.fastq")
           .set { samples_ch }
    debugView(samples_ch, 'RAW_SAMPLES')

    // 2) FASTQC: access named outputs using .out.emit_name
    fastqc(samples_ch)
    def qc_ch = fastqc.out.report
    def qc_log = fastqc.out.log
    debugView(qc_log, 'FASTQC')

    // 3) Trimming on QCed reads
    trimming(qc_ch)
    def trim_ch = trimming.out.report
    def trim_log = trimming.out.log
    debugView(trim_log, 'TRIMMING')

    // 4) Validation on trimmed reads
    validation(trim_ch, file("${projectDir}/modules/validation/validation.py"))
    def val_ch = validation.out.report
    def val_log = validation.out.log
    debugView(val_log, 'VALIDATION')

    // 5) MultiQC aggregates validation & QC reports
    def combo_ch = qc_ch.mix(val_ch)
    multiqc(combo_ch)
    def mq_ch = multiqc.out.report
    def mq_log = multiqc.out.log
    debugView(mq_log, 'MULTIQC')

    // 6) Publish final aggregated outputs
    publish_results(mq_ch)
}