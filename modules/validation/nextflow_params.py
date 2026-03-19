#!/usr/bin/env python3
"""
nextflow_params.py
==================
Utilities to build and write a validated_params.json file after the validation
step completes. The file is consumed by the main Nextflow pipeline via
``-params-file validated_params.json``, overriding the defaults defined in
``nextflow.config`` and ``nextflow_schema.json``.

Produced parameters
-------------------
validated_inputs (hidden, set by validation):
  ref_fasta_validated  – absolute path to the validated reference FASTA
  mod_fasta_validated  – absolute path to the validated modified FASTA
  mod_fasta_avail      – True when a validated modified genome is present

general_options (pipeline execution switches):
  run_ref_x_mod        – True when both ref_genome and mod_genome validation succeeded
  run_syri             – True when 1–n_sequence_limit contigs found (prokaryotic assembly);
                         False for >n_sequence_limit contigs (fragmented/eukaryotic) or when
                         no modified genome is available
  run_truvari          – always False (reserved for future use / manual override)
  run_illumina         – True when validated Illumina reads are present
  run_nanopore         – True when validated Nanopore (ONT) reads are present
  run_pacbio           – True when validated PacBio reads are present
  contig_file_size     – number of contig files from inter-genome characterisation
  run_vcf_annotation   – True when a validated GFF feature file is present

input_output_options (file paths, null when absent):
  pacbio_fastq         – path to the first validated PacBio FASTQ
  nanopore_fastq       – path to the first validated Nanopore FASTQ (or None)
  gff                  – path to the validated reference GFF/GFF3 file
"""

import json
from pathlib import Path


def build_params(validation_results: dict) -> dict:
    """
    Build a Nextflow params dict from validation results.

    Parameters
    ----------
    validation_results : dict
        Expected keys:
          ref_genome    – GenomeOutputMetadata or None
          mod_genome    – GenomeOutputMetadata or None
          genomexgenome – dict with 'contig_files' key (from genomexgenome_validation), or None
          reads         – List[ReadOutputMetadata], each with a .ngs_type attribute
          ref_feature   – FeatureOutputMetadata or None

    Returns
    -------
    dict
        Params dict ready to serialise as JSON for ``-params-file``.
        See module docstring for the full list of produced keys.
    """
    def _path(meta):
        return str(getattr(meta, "output_file", None)) if meta is not None else None

    ref_path = _path(validation_results.get("ref_genome"))
    mod_path = _path(validation_results.get("mod_genome"))
    gxg      = validation_results.get("genomexgenome") or {}
    reads    = validation_results.get("reads") or []
    gff_path = _path(validation_results.get("ref_feature"))

    # reads grouped by ngs_type — keep first validated path per type
    by_type = {}
    for r in reads:
        ngs = getattr(r, "ngs_type", None)
        if ngs and ngs not in by_type:
            by_type[ngs] = _path(r)

    gxg_metadata = gxg.get("metadata") or {}
    contig_file_size = len(gxg_metadata.get("contig_files", []))
    mod_n_sequence_limit = validation_results.get("mod_n_sequence_limit") or 5

    params = {
        # validated_inputs
        "mod_fasta_avail":    mod_path is not None,
        # general_options — pipeline switches
        "run_ref_x_mod":      ref_path is not None and mod_path is not None,
        "run_syri":           1 <= contig_file_size <= mod_n_sequence_limit,
        "run_truvari":        False,
        "run_illumina":       "illumina" in by_type,
        "run_nanopore":       "ont"      in by_type,
        "run_pacbio":         "pacbio"   in by_type,
        "contig_file_size":   contig_file_size,
        "run_vcf_annotation": gff_path is not None,
        # input_output_options — nullable paths
        "nanopore_fastq":     by_type.get("ont"),
    }
    if ref_path:
        params["ref_fasta_validated"] = ref_path
    if mod_path:
        params["mod_fasta_validated"] = mod_path
    if by_type.get("pacbio"):
        params["pacbio_fastq"] = by_type["pacbio"]
    if gff_path:
        params["gff"] = gff_path

    return params


def write_params(params: dict, path: Path) -> None:
    """Serialise params to JSON and write to path."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2, ensure_ascii=False)
    print(f"Validated params written to: {path}")
