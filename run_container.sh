#!/bin/bash
# run_container.sh - Script to run the EFSA Pipeline container interactively

# Build the Docker image
echo "Building EFSA Pipeline Docker image..."
docker build -t efsa-pipeline .

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "Error: Docker build failed"
    # Clean up the nextflow binary even if build failed
    exit 1
fi


# Check if data directories exist, create if they don't
mkdir -p data/inputs data/outputs

# Use default input directory
INPUT_MOUNT="-v $(pwd)/data/inputs:/EFSA_workspace/data/inputs"
echo "Using default input directory: ./data/inputs"


# Run the container interactively with volume mounts
echo "Starting EFSA Pipeline container..."
echo "You will be dropped into the container shell."
echo "Type 'exit' when you're done to return to your host system."
echo ""

docker run --privileged -d --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name efsa-pipeline-container \
    -v "$(pwd):/EFSA_workspace" \
    $INPUT_MOUNT \
    -v "$(pwd)/data/outputs:/EFSA_workspace/data/outputs" \
    efsa-pipeline

docker exec -it efsa-pipeline-container /bin/sh

echo "Container exited. You're back on your host system."