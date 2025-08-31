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

# Determine input directory to mount
if [ -n "$INPUT_DIR" ]; then
    # User provided custom input directory
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Error: Input directory '$INPUT_DIR' does not exist"
        exit 1
    fi
    INPUT_MOUNT="-v $(realpath "$INPUT_DIR"):/EFSA_workspace/data/inputs"
    echo "Using custom input directory: $INPUT_DIR"
else
    # Use default input directory
    INPUT_MOUNT="-v $(pwd)/data/inputs:/EFSA_workspace/data/inputs"
    echo "Using default input directory: ./data/inputs"
fi

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
    $INPUT_MOUNT \
    -v "$(pwd)/data/outputs:/EFSA_workspace/data/outputs" \
    efsa-pipeline

echo "Container exited. You're back on your host system."