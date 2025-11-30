# Installation Guide

## Prerequisites

- Python 3.10 or higher
- Conda or Miniconda
- Git

## Step-by-Step Installation

### 1. Install Main RNA-MAP Package

The optimization toolkit requires the main `rna-map` Python package. Install it first:

```bash
# Option 1: Install from GitHub (recommended)
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python

# Option 2: If published to PyPI
pip install rna-map
```

### 2. Create Conda Environment

```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate environment
conda activate rna-map-optimization
```

### 3. Install Main Package in Environment

After activating the environment, install the main `rna-map` package:

```bash
conda activate rna-map-optimization
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python
```

### 4. Verify Installation

```bash
# Check that rna-map is installed
python -c "import rna_map; print('✅ rna-map installed')"

# Check that optimization scripts work
python scripts/optimize_bowtie2_params_optuna.py --help
```

## Alternative: Install as Package

If this repository is packaged (has `pyproject.toml`):

```bash
# Install main package first
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python

# Install optimization package
cd rna-map-optimization
pip install -e .
```

## Troubleshooting

### Import Error: No module named 'rna_map'

**Solution**: Install the main RNA-MAP package:
```bash
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python
```

### Optuna Import Error

**Solution**: Install Optuna:
```bash
pip install optuna>=3.0
```

### Bowtie2 Not Found

**Solution**: Install Bowtie2 via conda:
```bash
conda install -c bioconda bowtie2
```

## Development Setup

For development, install with dev dependencies:

```bash
# Install main package
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python

# Install optimization package with dev dependencies
pip install -e ".[dev]"
```

