# Running the Nextflow Pipeline

## Overview

The main pipeline (`main.nf`) executes **all three workflows** in sequence:

- Short-read processing for Illumina data
- Long-read processing for PacBio/Oxford Nanopore data
- Reference vs modified genome comparison using SyRI


## Running Main Workflow

The pipeline runs validation and processing in a **single command**:

```bash
nextflow run main.nf --max_cpu $(nproc)
```

This first validates input data from `data/inputs/config.json`, then automatically runs the processing workflows (short-read, long-read, ref-vs-mod comparison) based on the validated inputs.

Validated files are written to `data/valid/run_YYYYMMDD_HHMMSS/` and `data/valid/validated_params.json` is produced for runtime consumption. See the [Validation Overview](../validation/OVERVIEW.md) for details on what this file contains.

By default, runtime logs are written to `data/outputs/logs/nextflow.log` and console output is kept quiet.

To use a custom configuration file:

```bash
nextflow run main.nf --config_json /path/to/config.json --max_cpu $(nproc)
```

If you want verbose console output for debugging, override quiet mode at launch time:

```bash
NXF_QUIET=false nextflow run main.nf --max_cpu $(nproc)
```

## Nextflow Options

| Option           | Description                                                                                                                                                                         |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-resume`        | Resume a pipeline run from the point where it previously stopped or failed.                                                                                                         |
| `-with-report`   | Generate a visual HTML report of the workflow execution, including task durations, resource usage, and statuses. The report is saved by default to `data/outputs/logs/report.html`. |
| `-with-timeline` | Generate a timeline visualization showing when each pipeline process started and finished. The timeline is saved by default to `data/outputs/logs/timeline.html`.                   |
| `-with-dag`      | Generate a directed acyclic graph (DAG) illustrating task dependencies in the workflow.                                                                                             |

> **Note:** `-resume` is not supported in this pipeline. The `work/` directory is automatically deleted at the end of every run — including failed runs — so Nextflow has no cached task outputs to resume from.


## Pipeline Options

| Option | Description |
|--------|-------------|
| `--config_json <path>` | Path to the input `config.json` (default: `data/inputs/config.json`) |
| `--out_dir <path>` | Output directory (default: `data/outputs`) |
| `--max_cpu <n>` | Maximum CPUs per process (default: `1`) |
| `--validation_level <level>` | Validation depth: `STRICT`, `TRUST`, or `MINIMAL` — overridden by `config.json` when set there |
| `--logging_level <level>` | Log verbosity: `DEBUG`, `INFO`, `WARNING`, or `ERROR` — overridden by `config.json` when set there |
| `--organism_type <type>` | Organism type: `PROKARYOTE` or `EUKARYOTE` — overridden by `config.json` when set there |
| `--force_defragment_ref` | Force reference defragmentation — unsupported workaround |
| `--help` | Print help message and exit |

For the full parameter list and priority rules, see [Configuration](configuration.md#parameters-params).


## Next Steps

After running the pipeline:

- Review the [Tool Parameter Reference](tool-parameters.md) for details on analysis thresholds and their rationale
- Review [Runtime Messages](runtime-messages.md) to understand pipeline progress
- Explore the [Output Directory Structure](directory-structures.md)
- Check [Pipeline Visualization](pipeline-visualization.md) for workflow diagrams
