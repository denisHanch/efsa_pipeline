#!/usr/bin/env nextflow

/*
  workflow: fasta_ref_x_mod.nf
  Purpose: Compare a reference FASTA to a modified/assembled FASTA, run nucmer -> deltaFilter -> show-coords -> syri,
           and convert resulting SV VCFs to TSV summary tables.

  Contract:
  - Inputs:
      - Channel of reference fasta (globbed from params.in_dir)
      - Channel of modified/assembled fasta (globbed from params.in_dir)
  - Outputs:
      - Channel of structural variant VCFs (per-assembly)
      - Channel of SV summary TSVs (per-assembly)
  - Success:
      - For each reference/modified pair, produce a VCF and a corresponding TSV summary.
*/

include { nucmer; deltaFilter; showCoords; syri } from "../modules/assembly.nf"
include { logWorkflowCompletion } from "../modules/logs.nf"
include { vcf_to_table } from "../modules/sv_calling.nf"


workflow ref_mod {
    take:
        ref_fasta
        mod_fasta
    main:
        log.info "▶ Running pipeline comparing reference and modified fasta."

        def prefix_name = "assembly"

        ref_mod_fasta = ref_fasta
            .combine(mod_fasta)
            .map { ref, mod -> tuple(prefix_name, ref, mod) }

        ref_mod_fasta | nucmer | set { delta }
        deltaFilter(prefix_name, delta) | set { filtered_delta }
        showCoords(prefix_name, filtered_delta) | set { coords }
        syri(ref_mod_fasta, coords, filtered_delta) | set { sv_vcf }
        vcf_to_table(sv_vcf)  | set { sv_tbl }

    emit: 
        sv_vcf
        sv_tbl
}

workflow {
    Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true) | set { ref_fasta }
    Channel.fromPath("$params.in_dir/*{assembled_genome,mod}.{fa,fna,fasta}", checkIfExists: true) | set { mod_fasta }

    ref_mod(ref_fasta, mod_fasta)
}

logWorkflowCompletion("reference to modified fasta comparision")