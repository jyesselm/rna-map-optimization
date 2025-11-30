#!/bin/bash
# Setup script for Optuna testing environment

set -e  # Exit on error

ENV_NAME="rna-map-optuna"
ENV_FILE="environment_optuna.yml"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Setting up Optuna testing environment"
echo "=========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    exit 1
fi

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' already exists."
    read -p "Do you want to remove and recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n "${ENV_NAME}" -y
    else
        echo "Using existing environment. Activate with: conda activate ${ENV_NAME}"
        exit 0
    fi
fi

# Create environment
echo "Creating conda environment from ${ENV_FILE}..."
cd "${SCRIPT_DIR}"
conda env create -f "${ENV_FILE}"

# Install rna_map package in development mode
echo ""
echo "Installing rna_map package in development mode..."
conda run -n "${ENV_NAME}" pip install -e .

echo ""
echo "=========================================="
echo "Environment created successfully!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate ${ENV_NAME}"
echo ""
echo "To test the optimization script, run:"
echo "  conda activate ${ENV_NAME}"
echo "  python3 scripts/optimize_bowtie2_params_optuna.py \\"
echo "    --fasta test/resources/case_1/test.fasta \\"
echo "    --fastq1 test/resources/case_1/test_mate1.fastq \\"
echo "    --fastq2 test/resources/case_1/test_mate2.fastq \\"
echo "    --n-trials 5 \\"
echo "    --output-dir optimization_optuna_test"
echo ""

