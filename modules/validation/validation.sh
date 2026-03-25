#!/bin/bash
# validate - Wrapper script for EFSA validation

# Default config path
DEFAULT_CONFIG="./data/inputs/config.json"
CONFIG_PATH="$DEFAULT_CONFIG"
FORCE_DEFRAGMENT=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --force-defragment-ref)
            FORCE_DEFRAGMENT=true
            shift
            ;;
        -h|--help)
            echo "Usage: validate [--config <path>] [--force-defragment-ref]"
            echo ""
            echo "Options:"
            echo "  --config <path>         Path to config file (default: ./data/inputs/config.json)"
            echo "  --force-defragment-ref  UNSUPPORTED WORKAROUND: merge all reference contigs"
            echo "                          into one sequence before validation. Use only when"
            echo "                          the reference is too fragmented for any workflow."
            echo "                          Results are NOT guaranteed to be meaningful."
            echo "  -h, --help              Show this help message"
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



# Clear valid directory before running validation
VALID_DIR="./data/valid"
echo "Clearing valid directory: $VALID_DIR"
rm -rf "${VALID_DIR:?}"/*

# Run the validation script with error handling
echo "Running EFSA validation with config: $CONFIG_PATH"
echo "---"

EXTRA_ARGS=""
[ "$FORCE_DEFRAGMENT" = true ] && EXTRA_ARGS="--force-defragment-ref"

if python3 ./modules/validation/main.py "$CONFIG_PATH" $EXTRA_ARGS; then
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