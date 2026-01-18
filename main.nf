#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from './workflows/fasta_ref_x_mod.nf'

include { long_ref as long_ref_pacbio; long_ref as long_ref_ont; long_ref as long_mod_pacbio; long_ref as long_mod_ont} from './workflows/long-read-ref.nf'
include { short_ref; short_ref as short_mod } from './workflows/short-read-ref.nf'

include { compare_unmapped; compare_unmapped as compare_unmapped_ont; compare_unmapped as compare_unmapped_pacbio } from './modules/mapping.nf'
include { nanoplot as nanoplo_pacbio } from './modules/qc.nf'
include { nanoplot as nanoplot_pacbio; nanoplot as nanoplot_ont } from './modules/qc.nf'
include { truvari_comparison } from './modules/compare_vcfs.nf'

include { qc } from './modules/subworkflow.nf'
include { describePipeline; logWorkflowCompletion; loadLongFastqFiles; loadShortFastqFiles } from './modules/logs.nf'

include { restructure_sv_table } from './modules/sv_calling.nf'


// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf -resume
    
    Options:
    -resume          Run pipeline from the point where it was interrupted or previously failed
    --out_dir        Output directory             (default: ${params.out_dir})
    --in_dir         Input directory              (default: ${params.in_dir})
    --registry       Docker/Singularity registry  (default: ${params.registry})
    --max_cpu        Maximum CPUs per process     (default: ${params.max_cpu})
    --log            Enable logging               (default: ${params.log})
    --clean_work     Remove workdir after success (default: ${params.clean_work})
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
    Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", checkIfExists: true) | set { mod_fasta }
    
    def ref_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def mod_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /mod_plasmid\.(fa|fna|fasta)$/ } ?: []

    def pacbio_files = file("${params.in_dir}/pacbio/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def ont_files = file("${params.in_dir}/ont/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def short_read_files = file("$params.in_dir/illumina/").listFiles()?.findAll { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ?: []


    // Reference to modified fasta comparison
    ref_mod(ref_fasta, mod_fasta)
    
    def vcfs = ref_mod.out.sv_vcf
    def sv_table = ref_mod.out.sv_table

    pipelines_running++

    // pacbio reads pipeline
    if (pacbio_files) {
        mapping_tag = "map-pb"

        pacbio_fastqs = loadLongFastqFiles("${params.in_dir}/pacbio/*.fastq.gz")
        
        // qc
        nanoplot_pacbio(pacbio_fastqs, "pacbio")
        
        describePipeline("long-pacbio", "reference & modified")
        long_mod_pacbio(pacbio_fastqs, mod_fasta, mapping_tag, mod_plasmid, "pacbio/long-mod")

        long_ref_pacbio(pacbio_fastqs, ref_fasta, mapping_tag, ref_plasmid, "pacbio/long-ref") 

        compare_unmapped_pacbio(long_ref_pacbio.out.unmapped_fastq, long_mod_pacbio.out.unmapped_fastq, "pacbio")

        vcfs = vcfs.mix(long_ref_pacbio.out.sv_vcf.map { it[1] })
        sv_table = sv_table.mix(long_ref_pacbio.out.sv_table.map { it[3] })

        pipelines_running++
    }

    // nanopore reads pipeline
    if (ont_files) {
        mapping_tag = "map-ont"
        
        ont_fastqs = loadLongFastqFiles("${params.in_dir}/ont/*.fastq.gz")

        // qc
        nanoplot_ont(ont_fastqs, "ont")

        describePipeline("long-ont", "reference & modified")
        long_mod_ont(ont_fastqs, mod_fasta, mapping_tag, mod_plasmid, "ont/long-mod")

        long_ref_ont(ont_fastqs, ref_fasta, mapping_tag, ref_plasmid, "ont/long-ref")

        compare_unmapped_ont(long_ref_ont.out.unmapped_fastq, long_mod_ont.out.unmapped_fastq, "ont")

        vcfs = vcfs.mix(long_ref_ont.out.sv_vcf.map { it[1] })
        sv_table = sv_table.mix(long_ref_ont.out.sv_table.map { it[3] })


        pipelines_running++
    }

    // short reads pipeline
    if (short_read_files) {    
        
        fastqs = loadShortFastqFiles(short_read_files)

        // QC and trimming module
        qc(fastqs, "illumina/qc_trimming") | set { trimmed }

        // Running mapping to the reference and modified fasta
        describePipeline("short", "reference & modified")
        short_mod(trimmed, mod_fasta, "illumina/short-mod", mod_plasmid)

        short_ref(trimmed, ref_fasta, "illumina/short-ref", ref_plasmid) 
        
        compare_unmapped(short_ref.out.unmapped_fastq, short_mod.out.unmapped_fastq, "short")

        vcfs = vcfs.mix(short_ref.out.sv_vcf)
        sv_table = sv_table.mix(short_ref.out.sv_table.map { it[3] })

        pipelines_running++
        
        restructure_sv_table(sv_table)

    }

    if (pipelines_running == 0) {
        log.error "❌  No valid inputs found. Skipping workflows.\n"
        exit 0
    } else {
         int comparisons = pipelines_running - 1
        String plural = comparisons == 1 ? "" : "s"
        log.info "ℹ️  Truvari: performing ${comparisons} comparison${plural}.\n"
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

logWorkflowCompletion("execution of main.nf")

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}