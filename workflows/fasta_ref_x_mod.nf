#!/usr/bin/env nextflow

/*
  workflow: fasta_ref_x_mod.nf
  Purpose: Compare a reference FASTA to a modified/assembled FASTA, run nucmer -> deltaFilter -> show-coords -> syri,
           and convert resulting SV VCFs to TSV summary tables.


  - Inputs:
      - Channel of reference fasta (globbed from params.in_dir)
      - Channel of modified/assembled fasta (globbed from params.in_dir)
  - Outputs:
      - Channel of structural variant VCFs (per-assembly)
      - Channel of SV summary TSVs (per-assembly)
  - Success:
      - For each reference/modified pair, produce a VCF and a corresponding TSV summary.
*/

include { nucmer; delta_filter; show_coords; syri; bcftools_concat; bgzip_tabix } from "../modules/assembly.nf"
include { logWorkflowCompletion; listFiles } from "../modules/logs.nf"
include { vcf_to_table_asm; create_empty_tbl } from "../modules/sv_calling.nf"

def executed = false

workflow ref_mod {
    take:
        ref_fasta
        contigs
    main:
        log.info "▶ Running pipeline comparing reference and modified fasta or reference and contigs."

        contig_files = listFiles("${params.in_dir}/", ".*contig.*\\.fasta")

        ref_mod_fasta = contigs
            .combine(ref_fasta)
            .map { contig, ref -> 
                def prefix = contig.baseName
                tuple(prefix, ref, contig)
            }

        ref_mod_fasta | nucmer | set { delta }

        delta | delta_filter | set { filtered_delta }

        filtered_delta | show_coords | set { coords }

        ref_mod_fasta.join(coords).join(filtered_delta) | set { syri_input_ch }

        syri(syri_input_ch) | set { sv_vcf }
        
        bgzip_tabix(sv_vcf) | set { sv_vcf_bgz }

        vcfs = sv_vcf_bgz.map { prefix, vcf, tbi -> vcf }.collect()

        tbis = sv_vcf_bgz.map { prefix, vcf, tbi -> tbi }.collect()

        bcftools_concat(vcfs, tbis) | vcf_to_table_asm | set { sv_tbl }

    emit: 
        sv_vcf
        sv_tbl
}


workflow.onComplete {
    if (executed) {
        if (workflow.success) {
            log.info "✅ The ref_mod processing pipeline completed successfully.\n"
        } else {
            log.error "❌ The ref_mod processing pipeline failed: ${workflow.errorReport}"
        }
    }
}
