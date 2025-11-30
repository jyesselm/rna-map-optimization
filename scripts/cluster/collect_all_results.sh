#!/bin/bash
# Convenience script to collect top results from all optimization runs

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_DIR="${PROJECT_ROOT}/optimization_results"
ENV_NAME="rna-map-optimization"

echo "=========================================="
echo "Collecting Top Results from All Cases"
echo "=========================================="
echo ""

# Check if results directory exists
if [ ! -d "${RESULTS_DIR}" ]; then
    echo "ERROR: Results directory not found: ${RESULTS_DIR}"
    echo "Please run optimization jobs first"
    exit 1
fi

# Check if conda environment exists
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "ERROR: Conda environment '${ENV_NAME}' not found"
    echo "Please run: bash ${SCRIPT_DIR}/setup_optimization_env.sh"
    exit 1
fi

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${ENV_NAME}

# Change to project directory
cd ${PROJECT_ROOT}

# Run collection script
python scripts/optimization/collect_top_results.py \
    --results-dir "${RESULTS_DIR}" \
    --top-n 100 \
    --output-dir "${RESULTS_DIR}/aggregated"

echo ""
echo "✅ Results collection complete!"
echo ""
echo "Check aggregated results in:"
echo "  ${RESULTS_DIR}/aggregated/"

