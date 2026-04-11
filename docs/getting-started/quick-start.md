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
> Make sure that the pipeline data is in the `data/inputs` folder.

### 2. Start the Docker Container

```bash
./run_container.sh
```

> **Important!**
>
> Create a configuration file `config.json` based on data in `data/inputs`.

### 3. Run the Pipeline

Run validation and processing in a single command:

```bash
nextflow run main.nf --max_cpu $(nproc)
```

This first validates input data from `data/inputs/config.json`, then automatically runs the processing workflows (short-read, long-read, ref-vs-mod comparison) based on the validated inputs.

> **Note:** Use `--config_json /path/to/config.json` to specify a custom configuration file.

## Next Steps

- Learn more about [Docker Setup](docker-setup.md)
- Configure your [Validation Settings](../validation/CONFIG_GUIDE.md)
- Understand [Pipeline Execution](../nextflow/running-pipeline.md)
