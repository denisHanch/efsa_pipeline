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
  validation_timestamp – timestamp of the validation run (YYYYMMDD_HHMMSS)

input_output_options (file paths, null / empty list when absent):
  ref_fasta_validated  – absolute path to the validated reference FASTA
  mod_fasta_validated  – absolute path to the validated modified FASTA
  ref_plasmid_fasta    – plasmid FASTA for reference (from explicit config or genome validator split)
  mod_plasmid_fasta    – plasmid FASTA for modified (from explicit config or GXG characterisation)
  gff                  – path to the validated reference GFF/GFF3 file
  illumina_fastqs      – list of validated Illumina FASTQ paths
  ont_fastqs           – list of validated Nanopore (ONT) FASTQ paths
  ont_bams             – list of validated Nanopore (ONT) BAM paths
  pacbio_fastqs        – list of validated PacBio FASTQ paths
  pacbio_bams          – list of validated PacBio BAM paths
  contig_files         – list of contig FASTA paths from inter-genome characterisation
"""

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional

_OPTIONAL_PATHS = {"ref_fasta_validated", "mod_fasta_validated", "gff", "ref_plasmid_fasta", "mod_plasmid_fasta"}


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
    validation_timestamp: str = ""
    # input_output_options — omitted from JSON when None
    ref_fasta_validated: Optional[str] = None
    mod_fasta_validated: Optional[str] = None
    ref_plasmid_fasta: Optional[str] = None
    mod_plasmid_fasta: Optional[str] = None
    gff: Optional[str] = None
    # input_output_options — lists, always present (may be empty)
    illumina_fastqs: List[str] = field(default_factory=list)
    ont_fastqs: List[str] = field(default_factory=list)
    ont_bams: List[str] = field(default_factory=list)
    pacbio_fastqs: List[str] = field(default_factory=list)
    pacbio_bams: List[str] = field(default_factory=list)
    contig_files: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Produce a dict for Nextflow -params-file JSON.
        Optional single-path fields are excluded when None; list fields are always included."""
        result = {
            "run_ref_x_mod": self.run_ref_x_mod,
            "run_truvari": self.run_truvari,
            "run_illumina": self.run_illumina,
            "run_nanopore": self.run_nanopore,
            "run_pacbio": self.run_pacbio,
            "contig_file_size": self.contig_file_size,
            "run_vcf_annotation": self.run_vcf_annotation,
            "validation_timestamp": self.validation_timestamp,
            "illumina_fastqs": self.illumina_fastqs,
            "ont_fastqs": self.ont_fastqs,
            "ont_bams": self.ont_bams,
            "pacbio_fastqs": self.pacbio_fastqs,
            "pacbio_bams": self.pacbio_bams,
            "contig_files": self.contig_files,
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
        if meta is None:
            return None
        value = getattr(meta, "output_file", None)
        if value is None:
            return None
        value = str(value).strip()
        return value if value and value != "None" else None

    ref_path         = _path(validation_results.get("ref_genome"))
    mod_path         = _path(validation_results.get("mod_genome"))
    ref_plasmid_path = _path(validation_results.get("ref_plasmid"))
    mod_plasmid_path = _path(validation_results.get("mod_plasmid"))
    gxg      = validation_results.get("genomexgenome") or {}
    gxg_metadata = gxg.get("metadata") or {}

    # Fall back to plasmids detected during genome validation when not explicitly configured
    if ref_plasmid_path is None:
        plasmid_filenames = getattr(validation_results.get("ref_genome"), "plasmid_filenames", None) or []
        if plasmid_filenames:
            ref_plasmid_path = str(plasmid_filenames[0])

    if mod_plasmid_path is None:
        mod_plasmid_path = gxg_metadata.get("plasmid_file") or None

    reads    = [
        r for r in (validation_results.get("reads") or [])
        if r is not None and _path(r) is not None
    ]
    gff_path = None if force_defragment_ref else _path(validation_results.get("ref_feature"))

    # collect all validated paths per ngs_type, split by format (fastq vs bam)
    fastqs_by_type: dict = {}
    bams_by_type: dict = {}
    for r in reads:
        ngs = getattr(r, "ngs_type", None)
        p = _path(r)
        if not ngs or not p:
            continue
        if getattr(r, "input_format", None) == "bam":
            bams_by_type.setdefault(ngs, []).append(p)
        else:
            fastqs_by_type.setdefault(ngs, []).append(p)

    contig_files = gxg_metadata.get("contig_files", [])

    ref_fragmented = getattr(validation_results.get("ref_genome"), 'fragmented', False)
    mod_fragmented = getattr(validation_results.get("mod_genome"), 'fragmented', False)

    return NextflowParams(
        # general_options — pipeline switches
        run_ref_x_mod=(ref_path is not None and mod_path is not None
                       and not ref_fragmented and not mod_fragmented
                       and gxg.get("passed", False)),
        run_truvari=False,
        run_illumina="illumina" in fastqs_by_type,
        run_nanopore=("ont" in fastqs_by_type or "ont" in bams_by_type),
        run_pacbio=("pacbio" in fastqs_by_type or "pacbio" in bams_by_type),
        contig_file_size=len(contig_files),
        run_vcf_annotation=gff_path is not None,
        validation_timestamp=datetime.now().strftime("%Y%m%d_%H%M%S"),
        # input_output_options
        ref_fasta_validated=ref_path,
        mod_fasta_validated=mod_path,
        ref_plasmid_fasta=ref_plasmid_path,
        mod_plasmid_fasta=mod_plasmid_path,
        gff=gff_path,
        illumina_fastqs=fastqs_by_type.get("illumina", []),
        ont_fastqs=fastqs_by_type.get("ont", []),
        ont_bams=bams_by_type.get("ont", []),
        pacbio_fastqs=fastqs_by_type.get("pacbio", []),
        pacbio_bams=bams_by_type.get("pacbio", []),
        contig_files=contig_files,
    )


def write_params(params: NextflowParams, path: Path) -> None:
    """Serialise params to JSON and write to path."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(params.to_dict(), fh, indent=2, ensure_ascii=False)
    print(f"Validated params written to: {path}")
