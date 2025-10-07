#!/bin/bash
# run_container.sh - Script to run the EFSA Pipeline container interactively

# Parse command line arguments
INPUT_DIR=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--input INPUT_DIRECTORY]"
            echo "  --input INPUT_DIRECTORY  Use custom input directory instead of default ./data/inputs"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

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
rm -f nxf-tmp*

# Check if data directories exist, create if they don't
mkdir -p data/inputs data/outputs

# Get the absolute path of the current directory
# Mount to the same absolute path inside the container: to prevent mixed up paths, when processes run in separate docker containers
WORKSPACE_PATH=$(pwd)

# Determine input directory to mount
if [ -n "$INPUT_DIR" ]; then
    # User provided custom input directory
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Error: Input directory '$INPUT_DIR' does not exist"
        exit 1
    fi
    INPUT_PATH=$(realpath "$INPUT_DIR")
    INPUT_MOUNT="-v $INPUT_PATH:$INPUT_PATH"
    echo "Using custom input directory: $INPUT_DIR (mounted at $INPUT_PATH)"
else
    # Use default input directory - it will be part of the workspace mount
    INPUT_MOUNT=""
    echo "Using default input directory: ./data/inputs"
fi

# Run the container interactively with volume mounts
echo "Starting EFSA Pipeline container..."
echo "Workspace mounted at: $WORKSPACE_PATH"
echo "You will be dropped into the container shell."
echo "Type 'exit' when you're done to return to your host system."
echo ""

docker run -it --rm \
    --network=host \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    -v /var/run/docker.sock:/var/run/docker.sock \
    --name efsa-pipeline-container \
    -w "$WORKSPACE_PATH" \
    -v "$WORKSPACE_PATH:$WORKSPACE_PATH" \
    $INPUT_MOUNT \
    efsa-pipeline

echo "Container exited. You're back on your host system."
