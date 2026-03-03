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
include { logWorkflowCompletion; listFiles } from "../modules/logs.nf"
include { vcf_to_table_asm; create_empty_tbl } from "../modules/sv_calling.nf"


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

        delta | deltaFilter | set { filtered_delta }

        filtered_delta | showCoords | set { coords }

        ref_mod_fasta.join(coords).join(filtered_delta) | set { syri_input_ch }

        syri(syri_input_ch) | set { sv_vcf }
        
        if (contig_files.size() == 1) {
            sv_vcf | vcf_to_table_asm | set { sv_tbl }
        } else {
            create_empty_tbl("assembly") | set { sv_tbl }
        }
        

    emit: 
        sv_vcf
        sv_tbl
}

workflow {
    def ref_fasta = Channel.fromPath("$params.in_dir/*{ref,reference_genome}.{fa,fna,fasta}", checkIfExists: true)

    def contigs_ch = Channel.fromPath("$params.in_dir/*_contig_*.fasta")

    ref_mod(ref_fasta, contigs_ch)
}


// logWorkflowCompletion("reference to modified fasta comparision")
