# EFSA Pipeline

## Docker Setup for Users

This guide shows you how to run the EFSA Pipeline in a Docker container with access to input/output folders.

## Prerequisites

- Docker installed on your system
- Git (if cloning the repository)

## Quick Start

### Option 1: Using the Run Script (Recommended)

1. Make sure the script is executable:
   ```bash
   chmod +x run_container.sh
   ```

2. Run the container:
   ```bash
   ./run_container.sh
   ```

3. You'll be dropped into the container shell where you can run CLI commands
4. Type `exit` when done to return to your host system

### Option 2: Using Docker Compose

1. Start the container:
   ```bash
   docker-compose up -d
   ```

2. Enter the container:
   ```bash
   docker-compose exec efsa-pipeline bash
   ```

3. When done, exit the container and stop it:
   ```bash
   exit  # exit the container shell
   docker-compose down
   ```

### Option 3: Manual Docker Commands

1. Build the image:
   ```bash
   docker build -t efsa-pipeline .
   ```

2. Run interactively:
   ```bash
   docker run -it --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name efsa-pipeline-container \
    -v "$(pwd):/EFSA_workspace" \
    -v "$(pwd)/data/inputs:/EFSA_workspace/data/inputs" \
    -v "$(pwd)/data/outputs:/EFSA_workspace/data/outputs" \
    efsa-pipeline
   ```

# Nextflow
From within dev container:

Please run the `docker login` command that is given in the efsa Slack channel - Nextflow notes.

Run pipeline for all the pipelines:
`nextflow run main.nf -resume`

Run pipeline for short-read processing:
`nextflow run workflows/short-read-ref.nf -resume`

Run pipeline for long-read processing:
`nextflow run workflows/long-read-ref.nf -resume`

Run pipeline for ref & mod fasta comparisons:
`nextflow run workflows/fasta_ref_x_mod.nf -resume`