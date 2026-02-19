# Quick Start

## Prerequisites

- Docker installed on your system
- Git for cloning the repository

## Setup Steps

### 1. Clone the Repository

Generate a GitHub token (see instructions in the repository link) and download the repository:

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
> Create configuration file `config.json` on `data/inputs`.

### 3. Run Input Validation

Running QC on the input data and processing data for the Nextflow pipeline to `data/valid` folder:

```bash
python3 ./modules/validation/main.py ./data/inputs/config.json
```

### 4. Execute the Pipeline

Use the GitHub token provided in the Slack channel, then start the pipeline:

```bash
nextflow run main.nf --max_cpu $(nproc) -resume
```

## Next Steps

- Learn more about [Docker Setup](docker-setup.md)
- Configure your [Validation Settings](../validation/CONFIG_GUIDE.md)
- Understand [Pipeline Execution](../nextflow/running-pipeline.md)
