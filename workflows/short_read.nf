#!/usr/bin/env nextflow
/*
  workflow: short_read.nf
  Purpose: Process short-read (Illumina) sequencing data: run QC and trimming,
           map reads to reference and optional plasmid sequences, call SNPs and
           structural variants (SVs) where applicable, collect unmapped reads,
           perform basic VCF QC, and convert SV VCFs to TSV summary tables.

  Contract:
  - Inputs:
      - Channel of trimmed short-read FASTQ files (produced by QC/trimming step)
      - Channel of reference FASTA to map against
      - out_folder_name: output prefix / label (controls which analyses run)
      - Optional plasmid FASTA (used to remap unmapped reads)
  - Outputs:
      - Channel of variant VCFs (SNP / SV) or Channel.empty() when skipped
      - Channel of unmapped FASTQ reads (after reference and optional plasmid mapping)
      - Channel of SV summary TSVs (per-run) or Channel.empty() when skipped
*/

include { multiqc } from "../modules/qc.nf"
include { calc_unmapped; calc_unmapped as calc_unmapped_plasmid; calc_total_reads; get_unmapped_reads; bwa_index; bwa_index as bwa_index_plasmid; get_unmapped_reads as get_unmapped_reads_plasmid } from "../modules/mapping.nf"
include { freebayes; bcftools_stats; bcftools_stats as bcftools_stats_plasmid } from "../modules/variant_calling.nf"
include { qc; mapping; sv; annotate_vcf; mapping as mapping_plasmid } from "../workflows/subworkflows.nf"
include { logUnmapped; logUnmapped as logUnmapped_plasmid; logWorkflowCompletion; loadShortFastqFiles } from "../modules/logs.nf"
include { vcf_to_table }  from "../modules/sv_calling.nf"


workflow short_read {
    take:
        trimmed
        fasta
        out_folder_name
        plasmid_fasta

    main:        
        // mapping to the reference
        bwa_index(fasta, out_folder_name) | set { fasta_index } 
        mapping(fasta, fasta_index, trimmed, out_folder_name) | set { indexed_bam }

        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // printout % unmapped reads
        calc_total_reads(indexed_bam) | set { total_reads }
        calc_unmapped(unmapped_fastq) | set { nreads }
        logUnmapped(nreads, total_reads, out_folder_name, "")

        // mapping reads to plasmid
        if (plasmid_fasta) {
            Channel.from(plasmid_fasta) | set { plasmid_fasta }

            bwa_index_plasmid(plasmid_fasta, "${out_folder_name}-plasmid") | set { unmapped_fasta_index } 
            mapping_plasmid(plasmid_fasta, unmapped_fasta_index, unmapped_fastq, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }

            calc_unmapped_plasmid(unmapped_fastq) | set { nreads }
            logUnmapped_plasmid(nreads, total_reads, "${out_folder_name}-plasmid", " against plasmid")
        }

         if (out_folder_name == "illumina/short-ref") { 
            // SNP & variant calling
            freebayes(fasta, indexed_bam, out_folder_name) | set { vcf }
            bcftools_stats(vcf, out_folder_name) | set { bcftools_out }
            
            // Optional annotation of vcf files
            def gff_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref*\.gff3$/ } ?: []

            if (gff_files) {
                Channel.fromPath(gff_files) | set { gff }
                annotate_vcf(fasta, gff, vcf, "gff", "gff3", out_folder_name) | set {qc_vcf}
            
                qc_vcf.mix(bcftools_out).collect() | set { qc_out }
                multiqc(qc_out, out_folder_name, "varint_calling")
            } else {
                multiqc(bcftools_out, out_folder_name, "varint_calling")
            }
        
            // SVs variant calling
            sv(fasta, indexed_bam, out_folder_name) | set { sv_vcf }
            vcf_to_table(sv_vcf) | set { sv_tbl }
        } else {
            sv_vcf = Channel.empty()
            sv_tbl = Channel.empty()
        }
    emit:
        sv_vcf
        unmapped_fastq
        sv_tbl

}

out_folder_name = "illumina/short-ref"

workflow { 

    log.info  "Processing files in directory: ${params.in_dir}"
    
    def plasmid_files = file("$params.in_dir").listFiles()?.findAll { it.name =~ /ref_plasmid\.(fa|fna|fasta)$/ } ?: []
    def short_read_files = file("$params.in_dir/illumina/").listFiles()?.findAll { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ?: []

    Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true) | set { fasta }
        
    fastqs = loadShortFastqFiles(short_read_files)

    qc(fastqs, "illumina/qc_trimming") | set { trimmed }

    short_read(trimmed, fasta, out_folder_name, plasmid_files)
}

workflow.onComplete {
    if (executed) {
        if (workflow.success) {
            log.info "✅ The short-read processing pipeline completed successfully.\n"
        } else {
            log.error "❌ The short-read processing pipeline failed: ${workflow.errorReport}"
        }
    }
}