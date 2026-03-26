#!/usr/bin/env python3
"""
nextflow_params_handler.py
==================
Utilities to build and write a validated_params.json file after the validation
step completes. The file is consumed by the main Nextflow pipeline via
``-params-file validated_params.json``, overriding the defaults defined in
``nextflow.config`` and ``nextflow_schema.json``.

Produced parameters
-------------------
general_options (pipeline execution switches):
  run_ref_x_mod        – True when both ref_genome and mod_genome and inter-genome validation succeeded;
                         False when any genome is too fragmented (eukaryote type or n_sequences > limit)
  run_truvari          – always False (reserved for future use / manual override)
  run_illumina         – True when validated Illumina reads are present
  run_nanopore         – True when validated Nanopore (ONT) reads are present
  run_pacbio           – True when validated PacBio reads are present
  contig_file_size     – number of contig files from inter-genome characterisation
  run_vcf_annotation   – True when a validated GFF feature file is present

input_output_options (file paths, null when absent):
  ref_fasta_validated  – absolute path to the validated reference FASTA
  mod_fasta_validated  – absolute path to the validated modified FASTA
  pacbio_fastq         – path to the first validated PacBio FASTQ
  nanopore_fastq       – path to the first validated Nanopore FASTQ (or None)
  gff                  – path to the validated reference GFF/GFF3 file
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

_OPTIONAL_PATHS = {"ref_fasta_validated", "mod_fasta_validated", "pacbio_fastq", "gff"}


@dataclass
class NextflowParams:
    # general_options
    run_ref_x_mod: bool = False
    run_truvari: bool = False
    run_illumina: bool = False
    run_nanopore: bool = False
    run_pacbio: bool = False
    contig_file_size: int = 0
    run_vcf_annotation: bool = False
    # input_output_options — always present (may be None)
    nanopore_fastq: Optional[str] = None
    # input_output_options — omitted from JSON when None
    ref_fasta_validated: Optional[str] = None
    mod_fasta_validated: Optional[str] = None
    pacbio_fastq: Optional[str] = None
    gff: Optional[str] = None

    def to_dict(self) -> dict:
        """Produce a dict for Nextflow -params-file JSON.
        Optional path fields are excluded when None; nanopore_fastq is always included."""
        result = {
            "run_ref_x_mod": self.run_ref_x_mod,
            "run_truvari": self.run_truvari,
            "run_illumina": self.run_illumina,
            "run_nanopore": self.run_nanopore,
            "run_pacbio": self.run_pacbio,
            "contig_file_size": self.contig_file_size,
            "run_vcf_annotation": self.run_vcf_annotation,
            "nanopore_fastq": self.nanopore_fastq,
        }
        for k in _OPTIONAL_PATHS:
            v = getattr(self, k)
            if v is not None:
                result[k] = v
        return result


def build_params(
    validation_results: dict,
    force_defragment_ref: bool = False,
) -> NextflowParams:
    """
    Build a Nextflow params dataclass from validation results.

    Parameters
    ----------
    validation_results : dict
        Expected keys:
          ref_genome    – GenomeOutputMetadata or None
          mod_genome    – GenomeOutputMetadata or None
          genomexgenome – dict with 'contig_files' key (from genomexgenome_validation), or None
          reads         – List[ReadOutputMetadata], each with a .ngs_type attribute
          ref_feature   – FeatureOutputMetadata or None
    force_defragment_ref : bool
        When True, GFF validation for the reference is skipped: ``gff`` is set to
        None and ``run_vcf_annotation`` is forced to False, because feature
        coordinates are no longer meaningful on a defragmented reference.

    Returns
    -------
    NextflowParams
        Params dataclass ready to serialise as JSON for ``-params-file``.
        See module docstring for the full list of produced keys.
    """
    def _path(meta):
        return str(getattr(meta, "output_file", None)) if meta is not None else None

    ref_path = _path(validation_results.get("ref_genome"))
    mod_path = _path(validation_results.get("mod_genome"))
    gxg      = validation_results.get("genomexgenome") or {}
    reads    = validation_results.get("reads") or []
    gff_path = None if force_defragment_ref else _path(validation_results.get("ref_feature"))

    # reads grouped by ngs_type — keep first validated path per type
    by_type = {}
    for r in reads:
        ngs = getattr(r, "ngs_type", None)
        if ngs and ngs not in by_type:
            by_type[ngs] = _path(r)

    gxg_metadata = gxg.get("metadata") or {}
    contig_file_size = len(gxg_metadata.get("contig_files", []))

    ref_fragmented = getattr(validation_results.get("ref_genome"), 'fragmented', False)
    mod_fragmented = getattr(validation_results.get("mod_genome"), 'fragmented', False)

    return NextflowParams(
        # general_options — pipeline switches
        run_ref_x_mod=(ref_path is not None and mod_path is not None
                       and not ref_fragmented and not mod_fragmented
                       and gxg.get("passed", False)),
        run_truvari=False,
        run_illumina="illumina" in by_type,
        run_nanopore="ont" in by_type,
        run_pacbio="pacbio" in by_type,
        contig_file_size=contig_file_size,
        run_vcf_annotation=gff_path is not None,
        # input_output_options — nullable paths
        nanopore_fastq=by_type.get("ont"),
        ref_fasta_validated=ref_path,
        mod_fasta_validated=mod_path,
        pacbio_fastq=by_type.get("pacbio"),
        gff=gff_path,
    )


def write_params(params: NextflowParams, path: Path) -> None:
    """Serialise params to JSON and write to path."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(params.to_dict(), fh, indent=2, ensure_ascii=False)
    print(f"Validated params written to: {path}")
