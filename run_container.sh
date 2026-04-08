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
mkdir -p data/inputs data/outputs data/outputs/tables/csv_per_sv_summary

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
echo "  gmo-pipeline                              - Run validation + Nextflow with defaults"
echo "  gmo-pipeline --config <path>              - Custom config file"
echo "  gmo-pipeline --threads <n|auto>           - Validation threads"
echo "  gmo-pipeline --validation-level <level>   - strict | trust | minimal (default: trust)"
echo "  gmo-pipeline --logging-level <level>      - DEBUG | INFO | WARNING | ERROR (default: INFO)"
echo "  gmo-pipeline --type <type>                - prokaryote | eukaryote (default: prokaryote)"
echo "  gmo-pipeline --force-defragment-ref       - UNSUPPORTED: merge fragmented reference contigs"
echo "  gmo-pipeline --profile <profile>          - Nextflow profile (default: standard)"
echo "  gmo-pipeline --resume                     - Resume previous Nextflow run"
echo "  gmo-pipeline --max-cpu <n>                - Max CPUs per Nextflow process"
echo "  gmo-pipeline --out-dir <path>             - Nextflow output directory"
echo "  gmo-pipeline --no-clean-work              - Keep Nextflow work directory"
echo ""
echo "  validate                              - Run validation only with default config"
echo "  validate --config <path>              - Run validation with custom config"
echo "  validate --threads <n|auto>           - Set number of threads"
echo "  validate --validation-level <level>   - strict | trust | minimal (default: trust)"
echo "  validate --logging-level <level>      - DEBUG | INFO | WARNING | ERROR (default: INFO)"
echo "  validate --type <type>                - prokaryote | eukaryote (default: prokaryote)"
echo "  validate --force-defragment-ref       - UNSUPPORTED: merge fragmented reference contigs"
echo ""
echo "Type 'exit' when you're done to return to your host system."
echo ""

docker run --privileged --init -d --rm \
    -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -v /usr/share/ca-certificates:/usr/share/ca-certificates:ro \
    --name "$CONTAINER_NAME" \
    -w "$WORKSPACE_PATH" \
    -v "$WORKSPACE_PATH:$WORKSPACE_PATH" \
    -e "NXF_LOG_FILE=${WORKSPACE_PATH}/data/outputs/logs/nextflow.log" \
    $INPUT_MOUNT \
    efsa-pipeline \
    sh -c "dockerd-entrypoint.sh & tail -f /dev/null"

# Execute interactive bash shell (not sh)
docker exec -it "$CONTAINER_NAME" /bin/bash

# Stop the container (--rm flag will automatically remove it)
docker stop "$CONTAINER_NAME" 2>/dev/null

echo ""
echo "Container exited and removed. You're back on your host system."
