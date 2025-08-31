#!/bin/bash
# run_container.sh - Script to run the EFSA Pipeline container interactively

# Build the Docker image
echo "Building EFSA Pipeline Docker image..."
docker build -t efsa-pipeline .

# Check if data directories exist, create if they don't
mkdir -p data/inputs data/outputs

# Run the container interactively with volume mounts
echo "Starting EFSA Pipeline container..."
echo "You will be dropped into the container shell."
echo "Type 'exit' when you're done to return to your host system."
echo ""

docker run -it --rm \
    --name efsa-pipeline-container \
    -v "$(pwd):/EFSA_workspace" \
    -v "$(pwd)/data/inputs:/EFSA_workspace/data/inputs" \
    -v "$(pwd)/data/outputs:/EFSA_workspace/data/outputs" \
    efsa-pipeline

echo "Container exited. You're back on your host system."