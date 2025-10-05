#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'

include { long_ref; long_ref as long_ref_ont } from './workflows/long-read-ref.nf'
include { long_mod; long_mod as long_mod_ont } from './workflows/long-read-mod.nf'

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

def describePipeline = { read_type, fasta_type, mod_fasta = null ->
    def msg = "▶ Running pipeline processing ${read_type} reads - mapping "
    if( mod_fasta )
        msg += "unmapped reads to the ${fasta_type} fasta."
    else
        msg += "to the ${fasta_type} fasta."
    return msg
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}

def pipelines_running = 0

out_folder_name = "final_vcf"



workflow {
    // Inputs
    Channel.fromPath("$params.in_dir/*ref.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*mod.{fa,fna,fasta}", checkIfExists: true) | set { mod_fasta }

    Channel.fromPath("${params.in_dir}/pacbio/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { pacbio_fastqs }

    Channel.fromPath("${params.in_dir}/ont/*_subreads.fastq.gz")
        .map { file -> 
            def name = file.baseName.replaceFirst('.fastq', '')
            return [name, file]
        }
        .set { ont_fastqs }

    Channel.fromPath("$params.in_dir/illumina/*.fastq.gz") | set { short_fastqs }

        if (pacbio_fastqs) {
            mapping_tag = "map-pb"

            if (params.map_to_mod_fa) {
                log.info describePipeline("long-pacbio", "modified", mod_fasta)
                long_mod(pacbio_fastqs, ref_fasta, mod_fasta, mapping_tag)
            } else {
                log.info describePipeline("long-pacbio", "reference")
                long_ref(pacbio_fastqs, ref_fasta, mapping_tag)
            }

            pipelines_running++
        }
        if (ont_fastqs) {
            mapping_tag = "map-ont"
            
            if (params.map_to_mod_fa) {
                log.info describePipeline("long-ont", "modified", mod_fasta)
                long_mod_ont(ont_fastqs, ref_fasta, mod_fasta, mapping_tag)
            } else {
                log.info describePipeline("long-ont", "reference")
                long_ref_ont(ont_fastqs, ref_fasta, mapping_tag)
            }

            pipelines_running++
        }

    if (short_fastqs) {    

        short_fastqs
        | map {[(it.name =~ /^([^_]+)(_((S[0-9]+_L[0-9]+_)?R[12]_001|[12]))?\.fastq.gz/)[0][1], it]}
        | groupTuple(sort: true)
        | set { fastqs }

        // QC and trimming module
        qc(fastqs, out_folder_name) | set { trimmed }

        // Running mapping to the reference or modified fasta 
        if (params.map_to_mod_fa) {
            log.info describePipeline("short", "modified", mod_fasta)
            mod_ref(trimmed, ref_fasta, mod_fasta)
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
    }
}