# Running the Nextflow Pipeline

## Overview

The main pipeline (`main.nf`) executes **all three workflows** in sequence:

- Short-read processing for Illumina data
- Long-read processing for Pacbio/Oxford Nanopore data
- Reference vs modified genome comparison using SyRI

Each workflow can also be executed individually if required.

## Running Main Workflow

This executes all workflows based on the files located within a subfolder of `data/valid` folder:

```bash
nextflow run main.nf --max_cpu $(nproc) -resume
```

## Nextflow Options

| Option           | Description                                                                                                                                                                         |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-resume`        | Resume a pipeline run from the point where it previously stopped or failed.                                                                                                         |
| `-with-report`   | Generate a visual HTML report of the workflow execution, including task durations, resource usage, and statuses. The report is saved by default to `data/outputs/logs/report.html`. |
| `-with-timeline` | Generate a timeline visualization showing when each pipeline process started and finished. The timeline is saved by default to `data/outputs/logs/timeline.html`.                   |
| `-with-dag`      | Generate a directed acyclic graph (DAG) illustrating task dependencies in the workflow.                                                                                             |


## Pipeline Options

| Option         | Description                                                       | Default                 |
|----------------|-------------------------------------------------------------------|-------------------------|
| `--in_dir`     | Input directory                                                   | `data/valid`            |
| `--out_dir`    | Output directory                                                  | `data/outputs`          |
| `--registry`   | Docker/Singularity container registry                             | `ghcr.io/kate-simonova` |
| `--max_cpu`    | Maximum CPUs per process                                          | `1`                     |
| `--run_truvari`| Enables filtering of VCFs based on truth set                      | `false`                 |
| `--run_syri`   | Enables comparison between the assembly FASTA and reference FASTA | `true`                  |
| `--clean_work` | Remove work directory after successful run                        | `true`                  |
| `--help`       | Display help message                                              | –                       |


## Next Steps

After running the pipeline:

- Review [Runtime Messages](runtime-messages.md) to understand pipeline progress
- Explore the [Output Directory Structure](directory-structures.md)
- Check [Pipeline Visualization](pipeline-visualization.md) for workflow diagrams
