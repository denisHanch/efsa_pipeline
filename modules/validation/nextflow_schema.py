#!/usr/bin/env python3
"""
nextflow_schema.py
==================
Static base schema for the EFSA pipeline Nextflow parameters, plus utilities
to enrich it with per-run validation output and write it to disk.

The BASE_SCHEMA dict mirrors the structure required by nf-schema@2.4.2
(JSON Schema Draft 2020-12). After validation runs, build_schema() populates
the validated_inputs section with actual file paths and pipeline flags, then
write_schema() serialises the result to nextflow_schema.json at the project root.
"""

import copy
import json
from pathlib import Path

# ---------------------------------------------------------------------------
# Static base schema
# ---------------------------------------------------------------------------

BASE_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://example.com/nextflow_schema.json",
    "title": "Automated Prokaryotic Genome Assembly Analysis Pipeline",
    "description": "Schema for pipeline parameters",
    "type": "object",
    "$defs": {
        "input_output_options": {
            "title": "Input/output options",
            "type": "object",
            "properties": {
                "in_dir": {
                    "type": "string",
                    "format": "directory-path",
                    "default": "data/valid",
                    "description": "Input directory for Nextflow containing reference and modified genome FASTA files."
                },
                "out_dir": {
                    "type": "string",
                    "format": "directory-path",
                    "default": "data/outputs",
                    "description": "Output directory where pipeline results will be written."
                },
                "log_dir": {
                    "type": "string",
                    "format": "directory-path",
                    "default": "data/outputs/logs",
                    "description": "Directory where log files will be stored."
                }
            }
        },
        "general_options": {
            "title": "General options",
            "type": "object",
            "properties": {
                "help": {
                    "type": "boolean",
                    "default": False,
                    "description": "Display help message and exit."
                },
                "max_cpu": {
                    "type": "integer",
                    "minimum": 1,
                    "default": 1,
                    "description": "Maximum number of CPUs available for pipeline processes."
                },
                "clean_work": {
                    "type": "boolean",
                    "default": False,
                    "description": "Clean the Nextflow work directory after successful completion."
                }
            }
        },
        "validated_inputs": {
            "title": "Validated inputs",
            "type": "object",
            "description": "File paths and pipeline availability flags produced by the validation step.",
            "properties": {
                "ref_fasta_validated": {
                    "type": "string",
                    "format": "file-path",
                    "hidden": True,
                    "description": "Validated reference genome FASTA path (set by validation step)."
                },
                "mod_fasta_validated": {
                    "type": "string",
                    "format": "file-path",
                    "hidden": True,
                    "description": "Validated modified genome FASTA path (set by validation step)."
                },
                "can_run_assembly": {
                    "type": "boolean",
                    "default": False,
                    "hidden": True,
                    "description": "True when both ref and mod genomes validated successfully; enables fasta_ref_x_mod workflow."
                }
            }
        }
    },
    "allOf": [
        {"$ref": "#/$defs/input_output_options"},
        {"$ref": "#/$defs/general_options"},
        {"$ref": "#/$defs/validated_inputs"}
    ]
}

# ---------------------------------------------------------------------------
# Schema builder
# ---------------------------------------------------------------------------

def build_schema(validation_results: dict) -> dict:
    """
    Return a deep copy of BASE_SCHEMA enriched with validation output.

    Parameters
    ----------
    validation_results : dict
        Expected keys:
          ref_genome  – GenomeOutputMetadata or None
          mod_genome  – GenomeOutputMetadata or None

    Returns
    -------
    dict
        Schema with validated_inputs.properties defaults populated where
        validation succeeded.
    """
    schema = copy.deepcopy(BASE_SCHEMA)
    props = schema["$defs"]["validated_inputs"]["properties"]

    def _path(meta):
        return getattr(meta, "output_file", None) if meta is not None else None

    ref_path = _path(validation_results.get("ref_genome"))
    mod_path = _path(validation_results.get("mod_genome"))

    # Only set default when a validated path exists — nf-schema requires omitting
    # the default key entirely when no value should be set.
    if ref_path:
        props["ref_fasta_validated"]["default"] = ref_path
    if mod_path:
        props["mod_fasta_validated"]["default"] = mod_path

    props["can_run_assembly"]["default"] = bool(ref_path and mod_path)

    return schema

# ---------------------------------------------------------------------------
# Schema writer
# ---------------------------------------------------------------------------

def write_schema(schema: dict, path: Path) -> None:
    """Serialise schema to JSON and write to path."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(schema, fh, indent=2, ensure_ascii=False)
    print(f"Schema written to: {path}")
