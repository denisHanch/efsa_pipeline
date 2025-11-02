#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'

include { long_ref as long_ref_pacbio; long_ref as long_ref_ont } from './workflows/long-read-ref.nf'
include { long_mod as long_mod_pacbio; long_mod as long_mod_ont } from './workflows/long-read-mod.nf'

include { short_ref; short_ref as short_mod } from './workflows/short-read-ref.nf'

include { truvari_comparison } from './modules/compare_vcfs.nf'

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
    
    def ref_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def mod_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /mod_plasmid\.(fa|fna|fasta)$/ } ?: []

    def pacbio_files = file("${params.in_dir}/pacbio/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def ont_files = file("${params.in_dir}/ont/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def short_read_files = file("$params.in_dir/illumina/").listFiles()?.findAll { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ?: []


    // Reference to modified fasta comparison
    ref_mod(ref_fasta, mod_fasta)
    
    def vcfs = ref_mod.out.sv_vcf

    pipelines_running++

    // pacbio reads pipeline
    // if (pacbio_files) {
    //     mapping_tag = "map-pb"

    //     Channel.from(pacbio_files).map { file ->
    //             def name = file.baseName.replaceFirst(/\.fastq$/, '')
    //             return [name, file]
    //         }.set { pacbio_fastqs }

    //     if (params.map_to_mod_fa) {
    //         log.info describePipeline("long-pacbio", "modified", mod_fasta)
    //         long_mod_pacbio(pacbio_fastqs, ref_fasta, mod_fasta, mapping_tag)
    //     } else {
    //         log.info describePipeline("long-pacbio", "reference")
    //         long_ref_pacbio(pacbio_fastqs, ref_fasta, mapping_tag) 
    //     }
        
    //     long_ref_pacbio.out.sv_vcf
    //         .map { it[1] }  
    //         .set { sv_pacbio_vcf }

    //     vcfs = vcfs.mix(sv_pacbio_vcf)

    //     pipelines_running++
    // }

    // nanopore reads pipeline
    // if (ont_files) {
    //     mapping_tag = "map-ont"
        
    //     Channel.from(ont_files).map { file ->
    //             def name = file.baseName.replaceFirst(/\.fastq$/, '')
    //             return [name, file]
    //         }
    //         .set { ont_fastqs }


    //     if (params.map_to_mod_fa) {
    //         log.info describePipeline("long-ont", "modified", mod_fasta)
    //         long_mod_ont(ont_fastqs, ref_fasta, mod_fasta, mapping_tag)
    //     } else {
    //         log.info describePipeline("long-ont", "reference")
    //         long_ref_ont(ont_fastqs, ref_fasta, mapping_tag)
    //     }

    //     long_ref_ont.out.sv_vcf
    //         .map { it[1] }  
    //         .set { sv_ont_vcf }

    //     vcfs = vcfs.mix(sv_ont_vcf)

    //     pipelines_running++
    // }

    // short reads pipeline
    if (short_read_files) {    
        
        Channel.from(short_read_files).map { file ->
        def matcher = file.name =~ /^(.+?)(?:[_\.](S[0-9]+_L[0-9]+_)?(R[12]|[12]))?\.f(ast)?q\.gz$/
        if( matcher.matches() ) {
            [ matcher[0][1], file ]
            }
        }
        | filter { it }
        | groupTuple(sort: true)
        | set { fastqs }


        // QC and trimming module
        qc(fastqs, "short-ref") | set { trimmed }

        // Running mapping to the reference and modified fasta
        log.info describePipeline("short", "modified")
        short_mod(trimmed, mod_fasta, "short-mod", mod_plasmid)
    
        log.info describePipeline("short", "reference")
        short_ref(trimmed, ref_fasta, "short-ref", ref_plasmid) 

        vcfs = vcfs.mix(short_ref.out.sv_vcf)

        short_ref.out.unmapped_fastq.view()
        // compare_unmapped(short_ref.out.unmapped_fastq, short_mod.out.unmapped_fastq)

        pipelines_running++
    }

    if (pipelines_running == 0) {
        log.error "⚠ No valid inputs found. Skipping workflows."
        exit 0
    } else {
        log.info "Performing ${pipelines_running - 1} truvari comparison(s)."
    }

    if (pipelines_running >= 2) {

        vcfs = vcfs
        .flatten()
        .map { file ->
            def name = file.getFileName().toString().replaceFirst(/\.vcf$/, '')
            return [name, file]
        }

        truvari_comparison(ref_fasta, vcfs)
    }
}

logWorkflowCompletion("execution of main.nf", true)

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}