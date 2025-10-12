#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'

include { long_ref as long_ref_pacbio; long_ref as long_ref_ont } from './workflows/long-read-ref.nf'
include { long_mod as long_mod_pacbio; long_mod as long_mod_ont } from './workflows/long-read-mod.nf'

include { short_ref } from './workflows/short-read-ref.nf'
include { short_mod } from './workflows/short-read-mod.nf'

include { truvari_comparison } from './workflows/compare_vcfs.nf'

include { qc } from './modules/subworkflow.nf'
include { describePipeline; logWorkflowCompletion } from './modules/logs.nf'


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

workflow {
    // Inputs
    Channel.fromPath("$params.in_dir/*ref*.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", checkIfExists: true) | set { mod_fasta }

    def pacbio_files = file("${params.in_dir}/pacbio/").listFiles()?.findAll { it.name =~ /_subreads\.fastq\.gz$/ } ?: []

    def ont_files = file("${params.in_dir}/nanopore/").listFiles()?.findAll { it.name =~ /_subreads\.fastq\.gz$/ } ?: []

    def short_read_files = file("$params.in_dir/illumina/").listFiles()?.findAll { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ?: []

    if (pacbio_files) {
        mapping_tag = "map-pb"

        Channel.from(pacbio_files).map { file ->
                def name = file.baseName.replaceFirst(/\.fastq$/, '')
                return [name, file]
            }.set { pacbio_fastqs }

        if (params.map_to_mod_fa) {
            log.info describePipeline("long-pacbio", "modified", mod_fasta)
            long_mod_pacbio(pacbio_fastqs, ref_fasta, mod_fasta, mapping_tag)
        } else {
            log.info describePipeline("long-pacbio", "reference")
            long_ref_pacbio(pacbio_fastqs, ref_fasta, mapping_tag)
        }

        pipelines_running++
    }
    if (ont_files) {
        mapping_tag = "map-ont"
        
        Channel.from(ont_files).map { file ->
                def name = file.baseName.replaceFirst(/\.fastq$/, '')
                return [name, file]
            }
            .set { ont_fastqs }


        if (params.map_to_mod_fa) {
            log.info describePipeline("long-ont", "modified", mod_fasta)
            long_mod_ont(ont_fastqs, ref_fasta, mod_fasta, mapping_tag)
        } else {
            log.info describePipeline("long-ont", "reference")
            long_ref_ont(ont_fastqs, ref_fasta, mapping_tag)
        }

        pipelines_running++
    }

    if (short_read_files) {    
        
        Channel.from(short_read_files)
        .map { [(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it] }
        .groupTuple(sort: true)
        .set { fastqs }

        // QC and trimming module
        qc(fastqs, "short-ref") | set { trimmed }

        // Running mapping to the reference or modified fasta 
        if (params.map_to_mod_fa) {
            log.info describePipeline("short", "modified", mod_fasta)
            short_mod(trimmed, ref_fasta, mod_fasta)
        } else {
            log.info describePipeline("short", "reference")
            short_ref(trimmed, ref_fasta)
        }
        pipelines_running++
    }

    if (ref_fasta && mod_fasta) {
        log.info "▶ Running pipeline comparing reference and modified fasta."
        ref_mod()
        
        pipelines_running++
    }

    if (pipelines_running == 0) {
        log.error "⚠ No valid inputs found. Skipping workflows."
        exit 0
    } else {
        log.info "Performing ${pipelines_running - 1} truvari comparison(s)."
    }

    if (pipelines_running >= 2) {

        // Inputs
        Channel.fromPath("$params.out_dir/final_vcf/*.vcf").map { file -> 
            def name = file.baseName.replaceFirst('.vcf', '')
            return [name, file]
        }
        .set { vcfs }

        truvari_comparison(ref_fasta, vcfs)
    }
}

logWorkflowCompletion("execution of main.nf", true)

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}