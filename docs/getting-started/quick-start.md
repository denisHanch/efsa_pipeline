# Quick Start

## Prerequisites

- Docker installed on your system
- Git for cloning the repository

## Setup Steps

### 1. Clone the Repository

Download the repository:

```bash
git clone https://github.com/denisHanch/efsa_pipeline.git
```

> **Important!**
>
> Make sure that the data for pipelines are in the folder `data/inputs`.

### 2. Start the Docker Container

```bash
./run_container.sh
```

> **Important!**
>
> Create configuration file `config.json` based on data in `data/inputs`.

### 3. Run Input Validation

Running QC on the input data and processing data for the Nextflow pipeline to `data/valid` folder:

```bash
validate                 - Run validation with default config
validate --config <path> - Run validation with custom config
```

### 4. Execute the Pipeline

Start the pipeline:

```bash
nextflow run main.nf --max_cpu $(nproc) -params-file data/valid/validated_params.json -resume
```

## Next Steps

- Learn more about [Docker Setup](docker-setup.md)
- Configure your [Validation Settings](../validation/CONFIG_GUIDE.md)
- Understand [Pipeline Execution](../nextflow/running-pipeline.md)
