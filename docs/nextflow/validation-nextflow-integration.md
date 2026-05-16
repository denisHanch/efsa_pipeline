# Validation and Nextflow Integration

This page describes how integration between the validation layer and the Nextflow pipeline is currently implemented.

## Integration Pattern

At the moment, integration is done through a **JSON handoff**:

1. Nextflow runs the `validate` process.
2. `validation.sh` executes the Python validation package.
3. Validation produces `validated_params.json`.
4. The analysis workflow reads and parses that JSON.
5. Parsed values drive downstream workflow selection, file paths, and runtime behavior.

## Where It Happens

- Validation process definition: `modules/validate.nf`
- Validation entrypoint script: `validation.sh` is within the `ecomolegmo/validation:1.0.6` image
- JSON handoff artifact: `validated_params.json`
- Consumer side (Nextflow): analysis workflow that parses the JSON and maps values into pipeline logic


## Why JSON Is Used

The JSON file is the contract between validation.nf and analysis.nf:

- It decouples validation from downstream process definitions.
- It provides a single machine-readable artifact for validated inputs and parameters.
- It allows the analysis workflow to branch dynamically based on validated content.

## Current Implications

Because the integration is JSON parsing based:

- Field names and structure in `validated_params.json` are part of a de facto interface contract.
- Any schema change on the validation side must be reflected in the Nextflow parsing/consumption logic.
- Backward compatibility depends on keeping that JSON contract stable or versioned.

## Genome-size handoff for SV percentage columns

Validation now writes two optional internal parameters when genome sizes are available:

| Parameter | Consumer | Purpose |
|---|---|---|
| `ref_genome_size_bp` | `workflows/analysis.nf` → `restructure_sv_tbl` | Calculates `pct_of_ref_genome` in the final SV CSV tables. |
| `mod_genome_size_bp` | `workflows/analysis.nf` → `restructure_sv_tbl` | Calculates `pct_of_mod_genome` in the final SV CSV tables. |

These values are derived from validation metadata and passed through the workflow automatically. `restructure_sv_tbl` also receives the validated FASTA files and computes genome sizes as a fallback when the optional params are absent. Users do not pass genome sizes to `modules/utils/create_sv_output.py` manually.

## Related Documentation

- [Validation Process](validation-process.md)
- [Validation Overview](../validation/OVERVIEW.md)
- [Configuration Guide](../validation/CONFIG_GUIDE.md)
