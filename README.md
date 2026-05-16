# EFSA Pipeline

A Nextflow DSL2 pipeline for processing prokaryotic long-read, short-read, and reference-vs-modified genome data, with integrated validation, quality control, mapping, variant calling, and comprehensive structural variant analysis.

## Table of Contents
- [Quick Start](#quick-start)
- [Running the Pipeline](#running-the-pipeline)
- [Documentation](#documentation)
- [Attribution](#attribution)

## Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/denisHanch/efsa_pipeline.git
   ```

2. **Prepare your data**:
   - Place input files in `data/inputs/`
   - Create a configuration file `data/inputs/config.json` (see [Configuration Guide](docs/validation/CONFIG_GUIDE.md) for details)

3. **Start the Docker container**:
   ```bash
   chmod +x run_container.sh
   ./run_container.sh
   ```

4. **Run the pipeline**:
   ```bash
   nextflow run main.nf --max_cpu $(nproc)
   ```

The pipeline automatically runs validation, then processes short-read, long-read, and reference-vs-modified genome workflows based on your configuration.
The modified genome FASTA (`mod_fasta`) is optional after validation. When `mod_fasta` is not provided, the pipeline runs in reference-only mode and skips mod-specific mapping/comparison branches.

## Running the Pipeline

### Full Pipeline Run

```bash
nextflow run main.nf --max_cpu $(nproc)
```

Executes validation first, then runs only the workflow branches enabled by validated inputs.
If a modified genome FASTA is available, both reference and modified branches run; if it is missing, the pipeline runs reference-only processing.

### Individual Workflows

You can also run specific workflows directly:

```bash
# Short-read processing only
nextflow run workflows/short_read.nf --max_cpu $(nproc)

# Long-read processing only
nextflow run workflows/long_read.nf --max_cpu $(nproc)

# Reference vs modified comparison only
nextflow run workflows/fasta_ref_x_mod.nf --max_cpu $(nproc)
```

**Note:** `-resume` is not fully supported in this pipeline. The `work/` directory is automatically deleted after each run. See [Running the Pipeline](docs/nextflow/running-pipeline.md) for details.

### Configuration

All pipeline behavior is controlled via:

- **`data/inputs/config.json`**: Input file specifications and validation settings (see [Configuration Guide](docs/validation/CONFIG_GUIDE.md))
- **`nextflow.config`**: Pipeline parameters, process resources, container images, and profiles (see [Configuration Reference](docs/nextflow/configuration.md))

## Documentation

For detailed information on all aspects of the pipeline, see the [documentation folder](docs/):

### Getting Started
- [Docker Setup](docs/getting-started/docker-setup.md)
- [Quick Start Guide](docs/getting-started/quick-start.md)

### Nextflow Configuration & Execution
- [Configuration Reference](docs/nextflow/configuration.md)
- [Running the Pipeline](docs/nextflow/running-pipeline.md)
- [Validation Process](docs/nextflow/validation-process.md)
- [Tool Parameters Reference](docs/nextflow/tool-parameters.md)
- [Nextflow Subworkflows](docs/nextflow/subworkflows.md)

### Validation
- [Validation Overview](docs/validation/OVERVIEW.md)
- [Configuration Guide](docs/validation/CONFIG_GUIDE.md)
- [File Formats](docs/validation/FILE_FORMATS.md)
- [Settings](docs/validation/SETTINGS.md)
- [Performance Tips](docs/validation/PERFORMANCE.md)

### Outputs
- [Output Directory Guide](docs/outputs/index.md)
- [Short-Read Outputs](docs/outputs/illumina.md)
- [Long-Read Outputs](docs/outputs/long-reads.md)
- [Reference vs Modified Outputs](docs/outputs/fasta-ref-mod.md)
- [SV Table Format & Columns](docs/outputs/sv-tables.md)
- [Logs & Diagnostics](docs/outputs/logs.md)
- [Unmapped Statistics](docs/outputs/unmapped-stats.md)

## Pipeline Overview

The pipeline has three main processing branches:

1. **Short-read (Illumina)** - Quality control, trimming, mapping, variant calling (SNVs, SVs)
2. **Long-read (PacBio/ONT)** - Mapping and structural variant calling via multiple callers
3. **Reference vs Modified** - Genome comparison using MUMmer and SyRI

All outputs are aggregated into unified SV tables in `data/outputs/tables/`.

## Attribution

Portions of this codebase and documentation were written with the assistance of AI coding assistants, including Claude and GitHub Copilot.
