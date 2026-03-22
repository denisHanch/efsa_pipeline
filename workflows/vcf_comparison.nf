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
def executed = false
include { sort_vcf; index_vcf; truvari } from '../modules/variant_calling.nf'

workflow truvari_comparison {

    take:
        ref_fasta
        vcfs_ch

    main:
        executed = true

        sort_vcf(vcfs_ch) | index_vcf | set { indexed_vcfs }

        split_ch = indexed_vcfs.branch {
            assembly: it[0].endsWith('mod')
            others:  !it[0].endsWith('mod')
        }
    
        assembly_vcfs = split_ch.assembly.collect()
        vcf_pairs_ch  = assembly_vcfs.combine(split_ch.others)

        truvari(ref_fasta, vcf_pairs_ch)
}

workflow.onComplete {
    if (executed) {
 if (workflow.success) {
        log.info "✅ Truvari: the pipeline was executed successfully.\n"
    } else {
        log.error "❌ Truvari: the comparison of vcf files failed.\n"
        }
    }
}