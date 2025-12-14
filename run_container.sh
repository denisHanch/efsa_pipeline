#!/bin/bash
# run_container.sh - Script to run the EFSA Pipeline container interactively

# Build the Docker image
echo "Building EFSA Pipeline Docker image..."
docker build -t efsa-pipeline .

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "Error: Docker build failed"
    exit 1
fi

# Check if data directories exist, create if they don't
mkdir -p data/inputs data/outputs

# Use default input directory
WORKSPACE_PATH=$(pwd)
INPUT_MOUNT="-v $WORKSPACE_PATH/data/inputs:/EFSA_workspace/data/inputs"
echo "Using default input directory: ./data/inputs"

# Generate unique container name using timestamp and random number
CONTAINER_NAME="efsa-pipeline-$(date +%s)-$$"

# Run the container interactively with volume mounts
echo "Starting EFSA Pipeline container..."
echo "Container name: $CONTAINER_NAME"
echo "Workspace mounted at: $WORKSPACE_PATH"
echo "You will be dropped into the container shell."
echo ""
echo "Available commands:"
echo "  validate              - Run validation with default config"
echo "  validate --config <path> - Run validation with custom config"
echo ""
echo "Type 'exit' when you're done to return to your host system."
echo ""

docker run --privileged -d --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name "$CONTAINER_NAME" \
    -w "$WORKSPACE_PATH" \
    -v "$WORKSPACE_PATH:$WORKSPACE_PATH" \
    $INPUT_MOUNT \
    efsa-pipeline \
    dockerd-entrypoint.sh

# Execute interactive bash shell (not sh)
docker exec -it "$CONTAINER_NAME" /bin/bash

# Stop the container (--rm flag will automatically remove it)
docker stop "$CONTAINER_NAME" 2>/dev/null

echo ""
echo "Container exited and removed. You're back on your host system."