#!/bin/bash
# validate - Wrapper script for validation
set -euo pipefail

# Default config path
DEFAULT_CONFIG="./data/inputs/config.json"
CONFIG_PATH="$DEFAULT_CONFIG"

# Global options (all optional — config.json takes priority if set there)
THREADS=""
VALIDATION_LEVEL=""
LOGGING_LEVEL=""
ORG_TYPE=""
FORCE_DEFRAGMENT=false

VENV_PYTHON="/opt/validation-venv/bin/python"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --validation-level)
            VALIDATION_LEVEL="$2"
            shift 2
            ;;
        --logging-level)
            LOGGING_LEVEL="$2"
            shift 2
            ;;
        --type)
            ORG_TYPE="$2"
            shift 2
            ;;
        --force-defragment-ref)
            FORCE_DEFRAGMENT=true
            shift
            ;;
        -h|--help)
            echo "Usage: validate [--config <path>] [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --config <path>              Path to config file (default: ./data/inputs/config.json)"
            echo ""
            echo "Global options (config.json takes priority if the same option is set there):"
            echo "  --threads <n|auto>           Number of threads, or 'auto' for system default"
            echo "  --validation-level <level>   Validation level: strict, trust, minimal (default: trust)"
            echo "  --logging-level <level>      Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)"
            echo "  --type <type>                Organism type: prokaryote, eukaryote (default: prokaryote)"
            echo "  --force-defragment-ref       UNSUPPORTED WORKAROUND: merge all reference contigs"
            echo "                               into one sequence before validation. Use only when"
            echo "                               the reference is too fragmented for any workflow."
            echo "                               Results are NOT guaranteed to be meaningful."
            echo "  -h, --help                   Show this help message"
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if config file exists
if [ ! -f "$CONFIG_PATH" ]; then
    echo "Error: Config file not found: $CONFIG_PATH"
    exit 1
fi

if [ ! -x "$VENV_PYTHON" ]; then
    echo "Error: Validation venv Python not found: $VENV_PYTHON"
    exit 1
fi

# Create run-stamped staging directory instead of destructively clearing the output dir
RUN_ID="$(date +%Y%m%d_%H%M%S)"
STAGING_DIR="./run_${RUN_ID}"
mkdir -p "$STAGING_DIR"
echo "Validation outputs will be written under: $STAGING_DIR"
export VALIDATION_RUN_DIR="$STAGING_DIR"

# Run the validation script with error handling
echo "Running validation with config: $CONFIG_PATH"
echo "---"

if "$VENV_PYTHON" /opt/validation/main.py "$CONFIG_PATH"; then
    echo "---"
    echo "Validation completed successfully"
    exit 0
else
    EXIT_CODE=$?
    echo "---"
    echo "Error: Validation failed with exit code $EXIT_CODE"
    echo "Container will remain running for debugging"
    exit $EXIT_CODE
fi