# Docker Container Setup

This guide shows you how to run the EFSA Pipeline in a Docker container with access to input/output folders.

## Prerequisites

- Docker installed on your system
- Git (if cloning the repository)

## Setup Options

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

### Option 2: Manual Docker Commands

1. Build the image:
   ```bash
   docker build -t efsa-pipeline .
   ```

2. Run interactively:
   ```bash
   docker run --privileged -d --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name efsa-pipeline-container \
    -w $(pwd) \
    -v "$(pwd)/data/inputs:/EFSA_workspace/data/inputs" \
    efsa-pipeline

   docker exec -it efsa-pipeline-container /bin/sh
   ```

## What's Next?

After setting up the Docker container, proceed to:
- [Input Validation](../validation/OVERVIEW.md)
- [Pipeline Execution](../nextflow/running-pipeline.md)
