#!/bin/bash
# run_container.sh - Script to run the EFSA Pipeline container interactively

# Download Nextflow on the VM first
echo "Downloading Nextflow on VM..."
curl -fsSL https://get.nextflow.io | bash

# Check if nextflow was downloaded successfully
if [ ! -f "./nextflow" ]; then
    echo "Error: Failed to download Nextflow"
    exit 1
fi

echo "Nextflow downloaded successfully"

# Build the Docker image
echo "Building EFSA Pipeline Docker image..."
docker build -t efsa-pipeline .

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "Error: Docker build failed"
    # Clean up the nextflow binary even if build failed
    rm -f ./nextflow
    exit 1
fi

# Clean up the nextflow binary from VM after successful build
echo "Cleaning up Nextflow binary from VM..."
rm -f ./nextflow

# Check if data directories exist, create if they don't
mkdir -p data/inputs data/outputs

# Run the container interactively with volume mounts
echo "Starting EFSA Pipeline container..."
echo "You will be dropped into the container shell."
echo "Type 'exit' when you're done to return to your host system."
echo ""
docker run -it --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name efsa-pipeline-container \
    -v "$(pwd):/EFSA_workspace" \
    -v "$(pwd)/data/inputs:/EFSA_workspace/data/inputs" \
    -v "$(pwd)/data/outputs:/EFSA_workspace/data/outputs" \
    efsa-pipeline

echo "Container exited. You're back on your host system."