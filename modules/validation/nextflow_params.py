#!/usr/bin/env python3
"""
nextflow_params.py
==================
Utilities to build and write a validated_params.json file after the validation
step runs. Nextflow picks up this file via `-params-file validated_params.json`.
"""

import json
from pathlib import Path


def build_params(validation_results: dict) -> dict:
    """
    Build a params dict from validation results.

    Parameters
    ----------
    validation_results : dict
        Expected keys:
          ref_genome  – GenomeOutputMetadata or None
          mod_genome  – GenomeOutputMetadata or None

    Returns
    -------
    dict
        Params ready to serialise as JSON for `-params-file`.
    """
    def _path(meta):
        return str(getattr(meta, "output_file", None)) if meta is not None else None

    ref_path = _path(validation_results.get("ref_genome"))
    mod_path = _path(validation_results.get("mod_genome"))

    params = {
        "fasta_ref_x_mod": bool(ref_path and mod_path),
    }
    if ref_path:
        params["ref_fasta_validated"] = ref_path
    if mod_path:
        params["mod_fasta_validated"] = mod_path

    return params


def write_params(params: dict, path: Path) -> None:
    """Serialise params to JSON and write to path."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(params, fh, indent=2, ensure_ascii=False)
    print(f"Validated params written to: {path}")
