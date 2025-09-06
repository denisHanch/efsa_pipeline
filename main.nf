#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'
include { long_ref } from './workflows/long-read-ref.nf'
include { short_ref } from './workflows/short-read-ref.nf'
include { sortVcf; indexVcf; truvari } from './modules/variant_calling.nf'


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

def pipelines_running = 0
def long_ch = Channel.fromPath("$params.in_dir/tmp2/*_subreads.fastq.gz")
def short_ch = Channel.fromPath("$params.in_dir/*.fastq.gz")
def fasta_ch = new File(params.in_dir).listFiles().findAll { it.name.endsWith('.fasta') || it.name.endsWith('.fa') || it.name.endsWith('.fna') }


workflow {

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

        // Inputs
        Channel.fromPath("$params.out_dir/final_vcf/*.vcf").map { file -> 
            def name = file.baseName.replaceFirst('.vcf', '')
            return [name, file]
        }
        .set { vcfs }

        Channel.fromPath("$params.in_dir/tmp/*ref.{fa,fna,fasta}") | set { ref_fasta }

        log.info "✅ Comparing VCF files from pipelines"

        // Sorting and indexiing vcfs
        sortVcf(vcfs) | indexVcf | set { indexed_vcfs }

        // preprocessing Channel for truvari input
        split_ch = indexed_vcfs.branch {
            ref_mod: it[0] == "ref_x_modsyri"
            others:  it[0] != "ref_x_modsyri"
        }

        ref_mod_ch = split_ch.ref_mod.collect()
        others_ch  = split_ch.others

        vcf_pairs_ch = ref_mod_ch.combine(others_ch)
        
                
        truvari(ref_fasta, vcf_pairs_ch)
    }   
}