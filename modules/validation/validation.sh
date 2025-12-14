#!/bin/bash
# validate-efsa - Wrapper script for EFSA validation

# Default config path
DEFAULT_CONFIG="./data/inputs/config.json"
CONFIG_PATH="$DEFAULT_CONFIG"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: validate-efsa [--config <path>]"
            echo ""
            echo "Options:"
            echo "  --config <path>    Path to config file (default: ./data/inputs/config.json)"
            echo "  -h, --help         Show this help message"
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

# Run the validation script with error handling
echo "Running EFSA validation with config: $CONFIG_PATH"
echo "---"

if python3 ./modules/validation/main.py "$CONFIG_PATH"; then
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