#!/bin/bash
# validate - Wrapper script for validation
set -euo pipefail

# Default config path
DEFAULT_CONFIG="./data/inputs/config.json"
CONFIG_PATH="$DEFAULT_CONFIG"

VENV_PYTHON="/opt/validation-venv/bin/python"

# Validation global options (empty = not set, defer to config.json)
THREADS=""
VALIDATION_LEVEL=""
LOGGING_LEVEL=""
ORGANISM_TYPE=""
FORCE_DEFRAGMENT_REF=false

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
            ORGANISM_TYPE="$2"
            shift 2
            ;;
        --force-defragment-ref)
            FORCE_DEFRAGMENT_REF=true
            shift
            ;;
        -h|--help)
            echo "Usage: validate [--config <path>] [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --config <path>              Path to config file (default: ./data/inputs/config.json)"
            echo "  --threads <int>              Number of threads"
            echo "  --validation-level <level>   'STRICTER', 'TRUST', or 'MINIMAL'"
            echo "  --logging-level <level>      'DEBUG', 'INFO', 'WARNING', or 'ERROR'"
            echo "  --type <type>                'PROKARYOTE' or 'EUKARYOTE'"
            echo "  --force-defragment-ref       Merge fragmented reference contigs (unsupported workaround)"
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

# Build extra args array (only include options that were explicitly set)
EXTRA_ARGS=()
[[ -n "$THREADS"           ]] && EXTRA_ARGS+=(--threads "$THREADS")
[[ -n "$VALIDATION_LEVEL"  ]] && EXTRA_ARGS+=(--validation-level "$VALIDATION_LEVEL")
[[ -n "$LOGGING_LEVEL"     ]] && EXTRA_ARGS+=(--logging-level "$LOGGING_LEVEL")
[[ -n "$ORGANISM_TYPE"     ]] && EXTRA_ARGS+=(--type "$ORGANISM_TYPE")
[[ "$FORCE_DEFRAGMENT_REF" == "true" ]] && EXTRA_ARGS+=(--force-defragment-ref)

# Run the validation script with error handling
echo "Running validation with config: $CONFIG_PATH"
echo "---"

if "$VENV_PYTHON" /opt/validation/main.py "$CONFIG_PATH" "${EXTRA_ARGS[@]}"; then
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