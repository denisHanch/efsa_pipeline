#!/usr/bin/env nextflow

/*
  workflow: vcf_comparison.nf
  Purpose: Prepare and run pairwise VCF comparisons using Truvari. This workflow
           sorts and indexes VCFs, splits assembly-derived VCFs from others,
           pairs assemblies with other callsets, and executes Truvari comparisons
           to assess structural variant concordance.

  Contract:
  - Inputs:
      - ref_fasta: reference FASTA used by Truvari for sequence context
      - vcfs_ch: channel of tuples (name, vcf_path) or similar VCF channel produced
                by upstream workflows (e.g., assembly and long/short-read callers)
  - Outputs:
      - Truvari outputs are produced in configured output locations per comparison.
      - (No direct channel emit here; Truvari run produces files/dirs consumed by downstream reporting.)
*/

include { sort_vcf; index_vcf; truvari } from '../modules/variant_calling.nf'

workflow truvari_comparison {

    // Inputs
    take:
        ref_fasta
        vcfs_ch

    main:
        // Sorting and indexing VCFs
        sort_vcf(vcfs_ch) | index_vcf | set { indexed_vcfs }

        // Preprocessing channel for Truvari input
        split_ch = indexed_vcfs.branch {
            assembly: it[0] == "assembly"
            others:  it[0] != "assembly"
        }

        assembly_vcfs = split_ch.assembly.collect()
        vcf_pairs_ch  = assembly_vcfs.combine(split_ch.others)

        // Run Truvari comparison
        truvari(ref_fasta, vcf_pairs_ch)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "✅ Truvari: the comparison of vcf files finished successfully.\n"
    } else {
        log.err "❌ Truvari: the comparison of vcf files failed.\n"
    }
}