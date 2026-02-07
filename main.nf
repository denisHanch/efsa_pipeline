#!/usr/bin/env nextflow

// Include workflows
include { ref_mod } from "./workflows/fasta_ref_x_mod.nf"
include { long_read as long_ref_pacbio; long_read as long_ref_ont; long_read as long_mod_pacbio; long_read as long_mod_ont} from "./workflows/long_read.nf"
include { short_read as short_ref; short_read as short_mod } from "./workflows/short_read.nf"
include { truvari_comparison } from "./workflows/vcf_comparison.nf"
include { qc } from "./workflows/subworkflows.nf"

include { compare_unmapped; compare_unmapped as compare_unmapped_ont; compare_unmapped as compare_unmapped_pacbio } from "./modules/mapping.nf"
include { nanoplot as nanoplot_pacbio; nanoplot as nanoplot_ont } from "./modules/qc.nf"
include { restructure_sv_tbl; create_empty_tbl as create_ont_tbl; create_empty_tbl as create_asm_tbl; create_empty_tbl as create_pacbio_tbl; create_empty_tbl as create_short_tbl } from "./modules/sv_calling.nf"
include { describePipeline; logWorkflowCompletion; loadFastqFiles; loadShortFastqFiles } from "./modules/logs.nf"


// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf -resume
    
    Options:
    -resume          Run pipeline from the point where it was interrupted or previously failed (Nextflow built-in)
    --out_dir        Output directory             (default: ${params.out_dir})
    --in_dir         Input directory              (default: ${params.in_dir})
    --registry       Docker/Singularity registry  (default: ${params.registry})
    --max_cpu        Maximum CPUs per process     (default: ${params.max_cpu})
    --log            Enable logging               (default: ${params.log})
    --clean_work     Remove workdir after success (default: ${params.clean_work})
    -with-report     Generate HTML execution report (Nextflow built-in)
    -with-timeline   Produce timeline visualization (Nextflow built-in)
    -with-dag        Produce DAG of workflow (Nextflow built-in)
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
    def ref_fasta = Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true)
    def mod_fasta = Channel.fromPath("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", checkIfExists: true)

    def ref_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def mod_plasmid = file("$params.in_dir").listFiles()?.findAll { it.name =~ /mod_plasmid\.(fa|fna|fasta)$/ } ?: []

    def pacbio_files = file("${params.in_dir}/pacbio/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def ont_files = file("${params.in_dir}/ont/").listFiles()?.findAll { it.name =~ /\.fastq\.gz$/ } ?: []
    def short_read_files = file("$params.in_dir/illumina/").listFiles()?.findAll { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ?: []


    // Reference to modified fasta comparison - assembly pipeline
    ref_mod(ref_fasta, mod_fasta)
    
    def vcfs = ref_mod.out.sv_vcf
    def sv_tbl = ref_mod.out.sv_tbl

    pipelines_running++

    // pacbio reads pipeline
    if (pacbio_files) {
        mapping_tag = "map-pb"

        pacbio_fastqs = loadFastqFiles("${params.in_dir}/pacbio/*.fastq.gz")
        
        // qc
        nanoplot_pacbio(pacbio_fastqs, "pacbio")
        
        describePipeline("long-pacbio", "reference & modified")
        long_mod_pacbio(pacbio_fastqs, mod_fasta, mapping_tag, mod_plasmid, "pacbio/long-mod")

        long_ref_pacbio(pacbio_fastqs, ref_fasta, mapping_tag, ref_plasmid, "pacbio/long-ref") 

        compare_unmapped_pacbio(long_ref_pacbio.out.unmapped_fastq, long_mod_pacbio.out.unmapped_fastq, "pacbio")

        vcfs = vcfs.mix(long_ref_pacbio.out.sv_vcf)
        sv_tbl = sv_tbl.mix(long_ref_pacbio.out.sv_tbl)

        pipelines_running++
    } else {
        sv_tbl = sv_tbl.mix(create_pacbio_tbl("pb"))
    }

    // nanopore reads pipeline
    if (ont_files) {
        mapping_tag = "map-ont"
        
        ont_fastqs = loadFastqFiles("${params.in_dir}/ont/*.fastq.gz")

        // qc
        nanoplot_ont(ont_fastqs, "ont")

        describePipeline("long-ont", "reference & modified")
        long_mod_ont(ont_fastqs, mod_fasta, mapping_tag, mod_plasmid, "ont/long-mod")

        long_ref_ont(ont_fastqs, ref_fasta, mapping_tag, ref_plasmid, "ont/long-ref")

        compare_unmapped_ont(long_ref_ont.out.unmapped_fastq, long_mod_ont.out.unmapped_fastq, "ont")

        vcfs = vcfs.mix(long_ref_ont.out.sv_vcf)
        sv_tbl = sv_tbl.mix(long_ref_ont.out.sv_tbl)

        pipelines_running++
    } else {
        sv_tbl = sv_tbl.mix(create_ont_tbl("ont"))
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
        sv_tbl = sv_tbl.mix(short_ref.out.sv_tbl)

        pipelines_running++
    
    } else {
        sv_tbl = sv_tbl.mix(create_short_tbl("short"))
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

        truvari_comparison(ref_fasta, vcfs)
    }
    
    script = file("${workflow.projectDir}/modules/utils/create_sv_output_xlsx.py")

    def tbl_channel = sv_tbl.collect().map { list ->
    def asm = list.find { it.name.toLowerCase().contains("assembly") }
    def long_pb = list.find { it.name.toLowerCase().contains("pb") }
    def long_ont = list.find { it.name.toLowerCase().contains("ont") }
    def sht = list.find { it.name.toLowerCase().contains("short") }
        tuple(asm, long_ont, long_pb, sht)
    }

    restructure_sv_tbl(script, tbl_channel)
}

logWorkflowCompletion("execution of main.nf")

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}