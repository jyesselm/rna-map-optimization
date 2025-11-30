# RNA-MAP Optimization Toolkit

Parameter optimization toolkit for [RNA-MAP Nextflow workflows](https://github.com/jyesselm/rna_map_nextflow).

## Overview

This repository provides tools for optimizing Bowtie2 alignment parameters for RNA mutational profiling (MaP) experiments. It uses Bayesian optimization (Optuna) and grid search to find optimal parameter combinations that maximize alignment quality and mutation detection sensitivity.

**Features:**
- **Automated parameter optimization** using Optuna (Bayesian optimization)
- **Grid search** for systematic parameter exploration
- **Cluster-based optimization** for large-scale runs
- **Analysis tools** for interpreting results
- **Recommended parameters** based on extensive analysis

## Installation

### Prerequisites

1. **Install the main RNA-MAP package** (required dependency):

   ```bash
   # Option 1: Install from GitHub (recommended)
   pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python
   
   # Option 2: If published to PyPI
   pip install rna-map
   ```

2. **Install optimization toolkit**:

   ```bash
   # Clone this repository
   git clone https://github.com/jyesselm/rna-map-optimization.git
   cd rna-map-optimization
   
   # Install dependencies
   conda env create -f environment.yml
   conda activate rna-map-optimization
   
   # Or install via pip (if packaged)
   pip install -e .
   ```

### Quick Setup

```bash
# 1. Install main RNA-MAP package
pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python

# 2. Setup optimization environment
conda env create -f environment.yml
conda activate rna-map-optimization

# 3. Verify installation
python scripts/optimize_bowtie2_params_optuna.py --help
```

## Quick Start

### Use Recommended Parameters

Recommended parameters from optimization are available in:
- **[docs/recommended_params/best_parameters.txt](./docs/recommended_params/best_parameters.txt)** - Full parameter string
- **[docs/BEST_PARAMETERS.md](./docs/BEST_PARAMETERS.md)** - Detailed breakdown

### Run Your Own Optimization

#### Local Optimization (Optuna - Recommended)

```bash
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --n-trials 100 \
    --output-dir optimization_results
```

#### Grid Search Optimization

```bash
python scripts/optimize_bowtie2_params.py \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --output-dir optimization_results
```

#### Baseline Comparison

```bash
python scripts/run_baseline_test.py \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --output-dir baseline_results
```

#### Cluster-Based Optimization

For large-scale optimization on a cluster, see [scripts/README.md](./scripts/README.md) for detailed instructions.

## Directory Structure

```
rna-map-optimization/
├── README.md                      # This file
├── LICENSE                        # Apache 2.0 License
├── pyproject.toml                 # Python package configuration
├── environment.yml                # Conda environment
├── .gitignore                     # Git ignore rules
│
├── scripts/                       # Executable scripts
│   ├── optimize_bowtie2_params_optuna.py  # Bayesian optimization (recommended)
│   ├── optimize_bowtie2_params.py         # Grid search
│   ├── run_baseline_test.py               # Baseline comparison
│   ├── collect_top_results.py            # Result aggregation
│   ├── cluster/                           # Cluster scripts
│   └── container/                         # Container build scripts
│
├── docs/                          # Documentation
│   ├── README.md                 # Documentation index
│   ├── BEST_PARAMETERS.md        # Recommended parameters
│   ├── TOP_100_PARAMETER_ANALYSIS.md  # Analysis results
│   ├── recommended_params/      # Parameter files
│   ├── examples/                 # Example scripts
│   └── archive/                  # Historical documentation
│
├── config/                        # Configuration files
│   └── slurm_optimized.config    # SLURM config for cluster
│
└── test/                          # Test cases and results
    ├── test_cases/                # Small test data (committed)
    └── results/                   # Optimization results (gitignored)
```

## Recommended Parameters

From analysis of top 100 parameter combinations:

**Constants (use these always):**
- Seed length: 18
- Mismatch penalty: 6,2
- Gap penalties: read=8,4, ref=5,3
- Sensitivity: fast-local
- Min insert size: 50

**Recommended defaults:**
- Max insert size: 200
- Seed mismatches: 1
- Score minimum: L,10,0.2
- Repetitive effort: 4

See **[docs/TOP_100_PARAMETER_ANALYSIS.md](./docs/TOP_100_PARAMETER_ANALYSIS.md)** for complete analysis.

## Documentation

### Essential Guides
- **[docs/BEST_PARAMETERS.md](./docs/BEST_PARAMETERS.md)** - Recommended parameters and usage
- **[docs/TOP_100_PARAMETER_ANALYSIS.md](./docs/TOP_100_PARAMETER_ANALYSIS.md)** - Detailed analysis results
- **[scripts/README.md](./scripts/README.md)** - Cluster-based optimization guide

### API Reference
- Optimization scripts are command-line tools
- Can be imported as Python modules for programmatic use
- See script docstrings for API details

## Usage Examples

### Quick Optimization (Local)

```bash
conda activate rna-map-optimization
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test.fasta \
    --fastq1 test_R1.fastq \
    --fastq2 test_R2.fastq \
    --n-trials 50 \
    --output-dir results
```

### Cluster Optimization

See [scripts/README.md](./scripts/README.md) for detailed cluster workflow.

## Results and Analysis

Optimization results include:
- Best parameter combinations
- Quality scores and statistics
- Comparison with baseline
- Visualization of optimization history (Optuna)

## Relationship to Main Repository

This optimization toolkit depends on the main [rna-map-nextflow](https://github.com/jyesselm/rna_map_nextflow) repository:

- **Main repo**: Production workflow and pipeline code
- **This repo**: Parameter optimization and research tools

The optimization code uses the `rna-map` Python package from the main repository to execute pipelines and evaluate results.

## Contributing

Contributions are welcome! Please see the main repository for contribution guidelines.

## License

Apache License 2.0 - See [LICENSE](./LICENSE) file for details.

## Citation

If you use this optimization toolkit, please cite the main RNA-MAP repository.

## Support

For issues and questions:
- Optimization-specific: Open an issue in this repository
- Main workflow: Open an issue in [rna-map-nextflow](https://github.com/jyesselm/rna_map_nextflow)
