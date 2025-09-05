#!/usr/bin/env nextflow


// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'
include { long_ref } from './workflows/long-read-ref.nf'
include { short_ref } from './workflows/short-read-ref.nf'


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


def long_ch = Channel.fromPath("$params.in_dir/tmp2/tmp2/*_subreads.fastq.gz")
def short_ch = Channel.fromPath("$params.in_dir/*.fastq.gz")
def fasta_ch = new File(params.in_dir).listFiles().findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') }



workflow {
    def pipelines_running = 0

    if (long_ch) {
        log.info "▶ Running pipeline processing long reads."
        // long_ref()
        pipelines_running++
    }

    if (short_ch) {
        log.info "▶ Running pipeline processing short reads."
        // short_ref()
        pipelines_running++
    }

    if (fasta_ch.size() == 2) {
        log.info "▶ Running pipeline comparing reference and modified fasta."
        // ref_mod(fasta_ch)
        pipelines_running++
    }

    if (pipelines_running == 0) {
        log.warn "⚠ No valid inputs found. Skipping workflows."
        exit 0
    } else {
        log.info "✅ Total workflows started: ${pipelines_running}"
    }

    if (pipelines_running > 2) {

        log.info "✅ Comparing VCF files from pipelines: ${vcfs}"

        vcfs
        // sortVcf
        // indexVcf
        // truvari
    }   

}