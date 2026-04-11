#!/usr/bin/env nextflow

/*
  workflow: long_read.nf
  Purpose: Process long-read sequencing data (ONT / PacBio): map reads to a reference,
           extract and report unmapped reads, optionally remap unmapped reads to a plasmid
           reference, call structural variants (SVs) for long-read runs, and convert SV VCFs
           to tabular summaries.

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


include { sv_long; mapping_long; mapping_long as mapping_long_plasmid }  from "../workflows/subworkflows.nf"
include { logUnmapped; logUnmapped as logUnmapped_plasmid } from "../modules/logs.nf"
include { calc_unmapped as calc_unmapped_long; calc_unmapped as calc_unmapped_plasmid; calc_total_reads; get_unmapped_reads;get_unmapped_reads as get_unmapped_reads_plasmid; build_sv_flank_bed; mosdepth } from "../modules/mapping.nf"
include { samtools_index; vcf_to_table_long }  from "../modules/sv_calling.nf"

def executed = false

workflow long_read {

    take:
        fastqs
        fasta
        mapping_tag
        plasmid_fasta
        out_folder_name

    main:
        // mapping to the reference

        executed = true

        samtools_index(fasta, out_folder_name) | set { fai }
        mapping_long(fastqs, fasta, mapping_tag, out_folder_name) | set { indexed_bam }

        get_unmapped_reads(indexed_bam, out_folder_name) | set { unmapped_fastq }
        
        // printout % unmapped reads
        calc_total_reads(indexed_bam) | set { total_reads }
        calc_unmapped_long(unmapped_fastq) | set { nreads }
        logUnmapped(nreads, total_reads, out_folder_name, "")

        // mapping reads to plasmid & variant calling
        if (plasmid_fasta) {
            plasmid_fasta.flatten() | set { plasmid_fasta }

            mapping_long_plasmid(unmapped_fastq, plasmid_fasta, mapping_tag, "${out_folder_name}-plasmid") | set { unmapped_bam }
            get_unmapped_reads_plasmid(unmapped_bam, "${out_folder_name}-plasmid") | set { unmapped_fastq }

            calc_unmapped_plasmid(unmapped_fastq) | set { nreads }
            logUnmapped_plasmid(nreads, total_reads, "${out_folder_name}-plasmid", " against plasmid")
        }

        // SV calling against the reference
        if (out_folder_name == "ont/long-ref" || out_folder_name == "pacbio/long-ref") { 
            sv_long(fasta, fai, indexed_bam, mapping_tag, out_folder_name)
            vcf_to_table_long(mapping_tag, sv_long.out.merged_vcf) | set { sv_tbl_raw }

            sv_tbl_keyed = sv_tbl_raw.map { tsv ->
                def suffix = "_${mapping_tag}_sv_summary"
                def pair_id = tsv.baseName.replace(suffix, "")
                tuple(pair_id, tsv)
            }
            build_sv_flank_bed(sv_tbl_keyed) | set { sv_tbl_with_regions }
            mosdepth(indexed_bam.join(sv_tbl_with_regions), mapping_tag) | set { sv_tbl }
            sv_vcf     = sv_long.out.merged_vcf
            supp_reads = sv_long.out.supp_reads

        } else {
            sv_vcf     = Channel.empty()
            sv_tbl     = Channel.empty()
            supp_reads = Channel.empty()
        }


    emit:
        sv_vcf
        unmapped_fastq
        sv_tbl
        supp_reads
}

workflow.onComplete {
    if (executed) {
        if (workflow.success) {
            log.info "✅ The long-read processing pipeline completed successfully.\n"
        } else {
            log.error "❌ The long-read processing pipeline failed: ${workflow.errorReport}"
        }
    }
}
