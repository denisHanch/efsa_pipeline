#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'
include { long_ref } from './workflows/long-read-ref.nf'
include { short_ref } from './workflows/short-read-ref.nf'
include { sortVcf; indexVcf} from './modules/varian'


def debugView = { ch, name ->
   ch
    .map { line -> "[${name}]: $line" }
    .view()
}

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --workflow [basic|advanced]
    
    Options:
    --workflow     Choose workflow: 'basic' or 'advanced' (default: basic)
    --outdir       Output directory (default: results)
    """.stripIndent()
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}


def long_ch = Channel.fromPath("$params.in_dir/tmp2/*_subreads.fastq.gz")
def short_ch = Channel.fromPath("$params.in_dir/*.fastq.gz")
def fasta_ch = new File(params.in_dir).listFiles().findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') || it.name.endsWith('.fna') }


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
        // ref_mod()
        pipelines_running++
    }

    if (pipelines_running == 0) {
        log.warn "⚠ No valid inputs found. Skipping workflows."
        exit 0
    } else {
        log.info "✅ Total workflows started: ${pipelines_running}"
    }

    if (pipelines_running > 2) {

        Channel.fromPath("$params.out_dir/final_vcf/*.vcf").map { file -> 
            def name = file.baseName.replaceFirst('.vcf', '')
            return [name, file]
        }
        .set { vcfs }

        log.info "✅ Comparing VCF files from pipelines: ${vcfs.view()}"

        sortVcf(vcfs) | indexVcf | set { indexed_vcfs }
        // truvari
    }   

}