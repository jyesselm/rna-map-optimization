# Installation Guide

## Prerequisites

- Python 3.10 or higher
- Conda or Miniconda
- Git

## Step-by-Step Installation

### 1. Create Conda Environment

```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate environment
conda activate rna-map-optimization
```

### 2. Install Optimization Package

The `rna-map-mini` dependency will be automatically installed from GitHub:

```bash
# Install the package (rna-map-mini will be installed automatically)
cd rna-map-optimization
pip install -e .
```

### 3. Verify Installation

```bash
# Check that rna-map-mini is installed
python -c "import rna_map_mini; print('✅ rna-map-mini installed')"

# Check that optimization command works
rna-map-optimize --help
```

## Quick Install

Simply install the package - all dependencies including `rna-map-mini` will be installed automatically:

```bash
cd rna-map-optimization
pip install -e .
```

The `rna-map-mini` package will be automatically installed from GitHub as specified in `pyproject.toml`.

## Troubleshooting

### Import Error: No module named 'rna_map_mini'

**Solution**: Install the rna-map-mini package:
```bash
pip install git+https://github.com/jyesselm/rna-map-mini.git
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
# Install optimization package with dev dependencies
# rna-map-mini will be installed automatically
pip install -e ".[dev]"
```

