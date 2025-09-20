#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'
include { long_ref } from './workflows/long-read-ref.nf'
include { short_ref } from './workflows/short-read-ref.nf'
include { mod_ref } from './workflows/short-read-mod.nf'
include { qc } from './modules/subworkflow.nf'
include { sortVcf; indexVcf; truvari } from './modules/variant_calling.nf'
 

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf -resume
    
    Options:
    -resume          Run pipeline from the point where it was interrupted or previously failed
    --out_dir        Output directory            (default: ${params.out_dir})
    --in_dir         Input directory             (default: ${params.in_dir})
    --registry       Docker/Singularity registry (default: ${params.registry})
    --log            Enable logging              (default: ${params.log})
    --help           Show this help message
    """.stripIndent()
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}

def pipelines_running = 0
def long_ch = Channel.fromPath("$params.in_dir/tmp2/*_subreads.fastq.gz")
def short_ch = Channel.fromPath("$params.in_dir/*.fastq.gz")
def fasta_ch = new File(params.in_dir).listFiles().findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') || it.name.endsWith('.fna') }

out_folder_name = "final_vcf"

workflow {

    if (long_ch) {
        log.info "▶ Running pipeline processing long reads."
        long_ref()

        pipelines_running++
    }

    if (short_ch) {    
        Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}") | set { ref_fasta }
        Channel.fromPath("$params.in_dir/*mod.{fa,fna,fasta}") | set { mod_fasta }
        
        short_ch
        | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
        | groupTuple(sort: true)
        | set { fastqs }

        // QC and trimming module
        qc(fastqs, out_folder_name) | set { trimmed }

        // Running mapping to the reference or modified fasta 
        if (params.map_to_mod_fa) {
            log.info "▶ Running pipeline processing short reads - mapping to the reference fasta."
            short_mod(trimmed, mod_fasta)
        } else {
            log.info "▶ Running pipeline processing short reads - mapping to the modified fasta."
            short_ref(trimmed, ref_fasta)
        }
        pipelines_running++
    }

    if (fasta_ch.size() == 2) {
        log.info "▶ Running pipeline comparing reference and modified fasta."
        ref_mod()

        pipelines_running++
    }

    if (pipelines_running == 0) {
        log.error "⚠ No valid inputs found. Skipping workflows."
        exit 0
    }
}