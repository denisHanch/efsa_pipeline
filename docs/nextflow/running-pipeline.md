# Running the Nextflow Pipeline

Please run the `docker login` command that is given in the efsa Slack channel - Nextflow notes.

## Overview

The main pipeline (`main.nf`) executes **all three workflows** in sequence:
- Short-read processing
- Long-read processing
- Reference vs modified genome comparison

Each workflow can also be executed individually if required.

## Running Main Workflow

This executes all workflows in sequence:

```bash
nextflow run main.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

## Available Nextflow Options

| Option        | Description |
|---------------|-------------|
| `-resume`     | Resumes a pipeline run from the point where it was interrupted or previously failed. |
| `-with-report` | Generates a visual HTML report of the workflow execution, including task durations, resource usage, and statuses. |
| `-with-timeline` | Produces a timeline visualization showing when each process in the pipeline started and finished. |
| `-with-dag`     | Generates a directed acyclic graph (DAG) showing task dependencies within the workflow. |

## Available Pipeline Options

| Option         | Description                                  | Default                |
| -------------- | -------------------------------------------- | ---------------------- |
| `--in_dir`     | Input directory                              | `data/valid`           |
| `--out_dir`    | Output directory                             | `data/outputs`         |
| `--registry`   | Docker/Singularity container registry        | `ghcr.io/kate-simonova`|
| `--max_cpu`    | Maximum CPUs per process                     | `1`                    |
| `--log`        | Enable logging                               | `true`                 |
| `--clean_work` | Remove work directory after successful run   | `true`                 |
| `--help`       | Display help message                         | –                      |

## Running Individual Pipelines

You can run each of the three sub-pipelines independently.

### Short-read Processing

For Illumina short-read data:

```bash
nextflow run workflows/short-read-ref.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

### Long-read Processing

For Oxford Nanopore / PacBio long reads:

```bash
nextflow run workflows/long-read-ref.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

### Reference vs Modified Genome Comparison

For comparing reference and modified FASTA assemblies:

```bash
nextflow run workflows/fasta_ref_x_mod.nf \
  -process.containerOptions "-u $(id -u):$(id -g)" \
  --max_cpu $(nproc) \
  -resume
```

## Next Steps

After running the pipeline:
- Review [Runtime Messages](runtime-messages.md) to understand pipeline progress
- Explore the [Output Directory Structure](directory-structures.md)
- Check [Pipeline Visualization](pipeline-visualization.md) for workflow diagrams
