#!/bin/bash
# Setup conda environment for Bowtie2 parameter optimization on cluster

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
ENV_NAME="rna-map-optimization"
ENV_FILE="${PROJECT_ROOT}/environment_optuna.yml"

echo "=========================================="
echo "Setting up Optimization Environment"
echo "=========================================="
echo ""
echo "Project root: ${PROJECT_ROOT}"
echo "Environment name: ${ENV_NAME}"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "ERROR: Conda or Mamba not found"
    echo "Please install Miniconda or Mambaforge:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Use mamba if available (faster), otherwise conda
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba (faster)"
else
    CONDA_CMD="conda"
    echo "Using conda"
fi

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' already exists."
    read -p "Do you want to remove and recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        ${CONDA_CMD} env remove -n "${ENV_NAME}" -y
    else
        echo "Using existing environment."
        echo ""
        echo "To activate: conda activate ${ENV_NAME}"
        echo "To update rna_map package:"
        echo "  conda activate ${ENV_NAME}"
        echo "  cd ${PROJECT_ROOT}/src/rna_map"
        echo "  python -m pip install -e ."
        exit 0
    fi
fi

echo ""
echo "Creating conda environment from ${ENV_FILE}..."
echo "This will install:"
echo "  - Python 3.10+"
echo "  - Bowtie2"
echo "  - pandas, numpy, tabulate, pyyaml"
echo "  - optuna (for Bayesian optimization)"
echo ""

# Create environment
cd "${PROJECT_ROOT}"
${CONDA_CMD} env create -f "${ENV_FILE}" -n "${ENV_NAME}"

echo ""
echo "Installing rna_map package in editable mode..."
echo ""

# Activate environment and install rna_map from src/rna_map
${CONDA_CMD} run -n "${ENV_NAME}" bash -c "cd ${PROJECT_ROOT}/src/rna_map && python -m pip install -e ."

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✅ Setup complete!"
    echo "=========================================="
    echo ""
    echo "Environment: ${ENV_NAME}"
    echo ""
    echo "To activate:"
    echo "  conda activate ${ENV_NAME}"
    echo ""
    echo "To verify installation:"
    echo "  conda activate ${ENV_NAME}"
    echo "  python -c 'import rna_map; print(\"rna_map installed successfully\")'"
    echo ""
else
    echo ""
    echo "⚠️  rna_map installation had issues"
    echo "You can install it manually:"
    echo "  conda activate ${ENV_NAME}"
    echo "  cd ${PROJECT_ROOT}/src/rna_map"
    echo "  python -m pip install -e ."
    echo ""
fi

