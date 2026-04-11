# Validation Process

The pipeline begins with a **validation step** that checks and prepares all input data before any analysis runs. This is implemented as a Nextflow process in `modules/validate.nf`.

## What It Does

The `validate` process wraps `validation.sh`, which invokes the Python validation package. For a detailed description of validation scenarios, preprocessing logic, and output structure, see the [Validation Overview](../validation/OVERVIEW.md).

In summary, the process:

1. Reads the input configuration (`config.json`)
2. Validates and cross-validates all input files
3. Copies validated files into `data/valid/run_YYYYMMDD_HHMMSS/`
4. Produces `validated_params.json` for the downstream analysis workflow

## Process Definition

```groovy
process validate {
    tag "validate"
    publishDir "${projectDir}/data/valid/", mode: 'copy', overwrite: true

    input:
    path config_json

    output:
    path 'validated_params.json', emit: params_json
    path 'run_*/**',              emit: run_dir

    script:
    """
    validation.sh --config ${config_json} --threads ${params.max_cpu}
    """
}
```

## Inputs and Outputs

| Direction | Name | Description |
|-----------|------|-------------|
| Input | `config_json` | Path to the JSON configuration file (default: `data/inputs/config.json`) |
| Output | `params_json` | `validated_params.json` — runtime parameters for the analysis workflow |
| Output | `run_dir` | Timestamped directory with validated/copied input files |

## Container

The validation process runs inside the `ecomolegmo/validation` Docker image, which bundles the Python validation package and all its dependencies. The image version is pinned with a sha256 digest in `nextflow.config`.

## How It Fits in the Pipeline

```text
main.nf
  └─ validate(config_json)
       ├─ params_json ──► analysis workflow (feeds all downstream workflows)
       └─ run_dir     ──► data/valid/run_*/  (validated input files)
```

The `analysis` workflow parses `validated_params.json` to determine which workflows to run and where the validated input files are located.

For full details on the validation configuration format, see [Configuration Guide](../validation/CONFIG_GUIDE.md). For validation levels and settings, see [Settings](../validation/SETTINGS.md).
