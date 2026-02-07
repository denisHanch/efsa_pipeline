#!/usr/bin/env nextflow

/*
  workflow: long_read.nf
  Purpose: Process long-read sequencing data (ONT / PacBio): map reads to a reference,
           extract and report unmapped reads, optionally remap unmapped reads to a plasmid
           reference, call structural variants (SVs) for long-read runs, and convert SV VCFs
           to tabular summaries.

  Contract:
  - Inputs:
      - Channel of long-read FASTQ files (e.g., from params.in_dir/ont or /pacbio)
      - Channel of reference FASTA to map against
      - mapping_tag: string used as mapping identifier (e.g., "map-ont" or "map-pb")
      - Optional plasmid FASTA (used to remap unmapped reads)
      - out_folder_name: output path prefix / label (controls conditional SV calling)
  - Outputs:
      - Channel of structural variant VCFs (per-run) or Channel.empty() when skipped
      - Channel of unmapped FASTQ reads (after reference and optional plasmid mapping)
      - Channel of SV summary tables (per-run) or Channel.empty() when skipped
*/

include { nanoplot; multiqc } from "../modules/qc.nf"
include { sv_long; mapping_long; mapping_long as mapping_long_plasmid; sv_long as sv_long_plasmid }  from "../workflows/subworkflows.nf"
include { logUnmapped; logUnmapped as logUnmapped_plasmid; logWorkflowCompletion; loadFastqFiles } from "../modules/logs.nf"
include { calc_unmapped; calc_unmapped as calc_unmapped_plasmid; calc_total_reads; get_unmapped_reads;get_unmapped_reads as get_unmapped_reads_plasmid } from "../modules/mapping.nf"
include { vcf_to_table_long }  from "../modules/sv_calling.nf"


workflow long_read {

    take:
        fastqs
        fasta
        mapping_tag
        plasmid_fasta
        out_folder_name

    main:
        // mapping to the reference
        mapping_long(fastqs, fasta, mapping_tag, out_folder_name) | set { indexed_bam }

        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // printout % unmapped reads
        calc_total_reads(indexed_bam) | set { total_reads }
        calc_unmapped(unmapped_fastq) | set { nreads }
        logUnmapped(nreads, total_reads, out_folder_name, "")

        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            Channel.from(plasmid_fasta) | set { plasmid_fasta }

            mapping_long_plasmid(unmapped_fastq, plasmid_fasta, mapping_tag, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }

            calc_unmapped_plasmid(unmapped_fastq) | set { nreads }
            logUnmapped_plasmid(nreads, total_reads, "${out_folder_name}-plasmid", " against plasmid")
        }

        // SV calling against the reference
        if (out_folder_name == "ont/long-ref" || out_folder_name == "pacbio/long-ref") { 
            sv_long(fasta, indexed_bam, mapping_tag, out_folder_name) | set { sv_vcf }
            vcf_to_table_long(mapping_tag, sv_vcf) | set { sv_tbl }

        } else {
            sv_vcf = Channel.empty()
            sv_tbl = Channel.empty()
        }

    emit:
        sv_vcf
        unmapped_fastq
        sv_tbl
}
