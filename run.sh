#!/usr/bin/env bash
# run.sh - Run validation followed by Nextflow pipeline
set -euo pipefail

VENV_PYTHON="/opt/validation-venv/bin/python"
DEFAULT_CONFIG="./data/inputs/config.json"

# Validation flags
CONFIG_PATH="$DEFAULT_CONFIG"
THREADS=""
VALIDATION_LEVEL=""
LOGGING_LEVEL=""
ORG_TYPE=""
FORCE_DEFRAGMENT=false

# Nextflow flags
NF_PROFILE="standard"
NF_RESUME=false
NF_MAX_CPU=""
NF_OUT_DIR=""
NF_IN_DIR=""
NF_CLEAN_WORK=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Validation options:
  --config <path>              Path to config file (default: $DEFAULT_CONFIG)
  --threads <n|auto>           Number of threads, or 'auto' for system default
  --validation-level <level>   strict | trust | minimal (default: trust)
  --logging-level <level>      DEBUG | INFO | WARNING | ERROR (default: INFO)
  --type <type>                prokaryote | eukaryote (default: prokaryote)
  --force-defragment-ref       UNSUPPORTED: merge fragmented reference contigs

Nextflow options:
  --profile <profile>          Nextflow profile (default: standard)
  --resume                     Resume previous run
  --max-cpu <n>                Maximum CPUs per process (default: 1)
  --out-dir <path>             Output directory (default: data/outputs)
  --in-dir <path>              Input directory (default: data/valid)
  --no-clean-work              Keep work directory after success

  -h, --help                   Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --config)             CONFIG_PATH="$2";       shift 2 ;;
        --threads)            THREADS="$2";           shift 2 ;;
        --validation-level)   VALIDATION_LEVEL="$2";  shift 2 ;;
        --logging-level)      LOGGING_LEVEL="$2";     shift 2 ;;
        --type)               ORG_TYPE="$2";          shift 2 ;;
        --force-defragment-ref) FORCE_DEFRAGMENT=true; shift ;;
        --profile)            NF_PROFILE="$2";        shift 2 ;;
        --resume)             NF_RESUME=true;          shift ;;
        --max-cpu)            NF_MAX_CPU="$2";        shift 2 ;;
        --out-dir)            NF_OUT_DIR="$2";        shift 2 ;;
        --in-dir)             NF_IN_DIR="$2";         shift 2 ;;
        --no-clean-work)      NF_CLEAN_WORK=false;    shift ;;
        -h|--help)            usage; exit 0 ;;
        *)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# Step 1: Validation
# ============================================================================
if [ ! -f "$CONFIG_PATH" ]; then
    echo "Error: Config file not found: $CONFIG_PATH"
    exit 1
fi

if [ ! -x "$VENV_PYTHON" ]; then
    echo "Error: Validation venv Python not found: $VENV_PYTHON"
    exit 1
fi

VALID_DIR="./data/valid"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
STAGING_DIR="${VALID_DIR}/run_${RUN_ID}"
mkdir -p "$STAGING_DIR"
echo "Validation outputs will be written under: $STAGING_DIR"
export VALIDATION_RUN_DIR="$STAGING_DIR"

VALIDATION_ARGS=()
[ "$FORCE_DEFRAGMENT" = true ] && VALIDATION_ARGS+=(--force-defragment-ref)
[ -n "$THREADS" ]            && VALIDATION_ARGS+=(--threads "$THREADS")
[ -n "$VALIDATION_LEVEL" ]   && VALIDATION_ARGS+=(--validation-level "$VALIDATION_LEVEL")
[ -n "$LOGGING_LEVEL" ]      && VALIDATION_ARGS+=(--logging-level "$LOGGING_LEVEL")
[ -n "$ORG_TYPE" ]           && VALIDATION_ARGS+=(--type "$ORG_TYPE")

echo "Running validation with config: $CONFIG_PATH"
echo "---"
if ! "$VENV_PYTHON" ./modules/validation/main.py "$CONFIG_PATH" "${VALIDATION_ARGS[@]}"; then
    EXIT_CODE=$?
    echo "---"
    echo "Error: Validation failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi
echo "---"
echo "Validation completed successfully"

# ============================================================================
# Step 2: Nextflow
# ============================================================================
PARAMS_FILE="${VALID_DIR}/validated_params.json"
if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Nextflow params file not found: $PARAMS_FILE"
    exit 1
fi

NF_ARGS=(-params-file "$PARAMS_FILE" -profile "$NF_PROFILE")
[ "$NF_RESUME" = true ] && NF_ARGS+=(-resume)
[ -n "$NF_MAX_CPU" ]    && NF_ARGS+=(--max_cpu "$NF_MAX_CPU")
[ -n "$NF_OUT_DIR" ]    && NF_ARGS+=(--out_dir "$NF_OUT_DIR")
[ -n "$NF_IN_DIR" ]     && NF_ARGS+=(--in_dir "$NF_IN_DIR")
[ "${NF_CLEAN_WORK:-}" = false ] && NF_ARGS+=(--clean_work false)

echo ""
echo "Running Nextflow pipeline..."
echo "---"
nextflow run main.nf "${NF_ARGS[@]}"
