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

1. **Install optimization toolkit**:

   ```bash
   # Clone this repository
   git clone https://github.com/jyesselm/rna-map-optimization.git
   cd rna-map-optimization
   
   # Install dependencies (rna-map-mini will be installed automatically)
   conda env create -f environment.yml
   conda activate rna-map-optimization
   
   # Install the package
   pip install -e .
   ```

### Quick Setup

```bash
# 1. Setup optimization environment
conda env create -f environment.yml
conda activate rna-map-optimization

# 2. Install package (rna-map-mini installed automatically)
pip install -e .

# 3. Verify installation
rna-map-optimize --help
```

## Quick Start

### Use Recommended Parameters

Recommended parameters from optimization are available in:
- **[docs/recommended_params/best_parameters.txt](./docs/recommended_params/best_parameters.txt)** - Full parameter string
- **[docs/BEST_PARAMETERS.md](./docs/BEST_PARAMETERS.md)** - Detailed breakdown

### Run Your Own Optimization

#### Example: Run with Test Case

Try the optimization with the included test case:

```bash
# Activate environment
conda activate rna-map-optimization

# Run optimization with test_cases/case_1
rna-map-optimize \
    --case-dir test_cases/case_1 \
    --n-trials 50 \
    --output-dir test_results

# The tool automatically detects:
# - test.fasta (reference)
# - test_mate1.fastq (R1 reads)
# - test_mate2.fastq (R2 reads)
```

Results will be saved to `test_results/case_1/` including:
- Best parameter combinations
- Optimization history visualizations
- Trial results and statistics

#### Local Optimization (Optuna - Recommended)

**Option 1: Single case directory (auto-detects files) - Recommended**
```bash
rna-map-optimize \
    --case-dir /path/to/case \
    --n-trials 100 \
    --output-dir optimization_results
```
The tool automatically finds `*.fasta`/`*.fa` and `*_R1*.fastq`/`*_R2*.fastq` files in the directory.

**Option 2: Explicit file paths**
```bash
rna-map-optimize \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --n-trials 100 \
    --output-dir optimization_results
```

**Option 3: Batch mode (multiple cases)**
```bash
rna-map-optimize \
    --cases-dir /path/to/cases \
    --n-trials 100 \
    --output-dir optimization_results
```
Processes all subdirectories as separate test cases.

**Note**: You can also use `python -m rna_map_optimization.cli` instead of `rna-map-optimize` (they're equivalent).

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

#### Alignment Verification (Secondary Check)

Verify alignments using alternative methods for small subsamples:

```bash
# Install optional verification dependencies
pip install biopython  # or: pip install parasail

# Verify alignments
python scripts/verify_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --fastq reads.fastq \
    --method biopython \
    --max-reads 100
```

See [docs/ALIGNMENT_VERIFICATION.md](./docs/ALIGNMENT_VERIFICATION.md) for details.

#### View Alignments (Human-Readable)

View alignments in a beautiful, human-readable format:

```bash
# Install optional rich library for colored output
pip install rich

# View alignments
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 10 \
    --min-mapq 20 \
    --output alignments.html
```

See [docs/VIEW_ALIGNMENTS.md](./docs/VIEW_ALIGNMENTS.md) for details.

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
├── src/                           # Source code
│   └── rna_map_optimization/      # Main package
│       ├── cli.py                 # Command-line interface
│       ├── optimizer.py           # Optuna-based optimization
│       ├── bowtie2.py             # Bowtie2 utilities
│       └── utils.py               # Utility functions
│
├── scripts/                       # Helper scripts and docs
│   ├── cluster/                   # Cluster job scripts
│   └── container/                 # Container build scripts
│
├── docs/                          # Documentation
│   ├── README.md                 # Documentation index
│   ├── BEST_PARAMETERS.md        # Recommended parameters
│   ├── ALIGNMENT_VERIFICATION.md # Alignment verification guide
│   ├── VIEW_ALIGNMENTS.md        # Alignment viewer guide
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
rna-map-optimize \
    --case-dir test_cases/case_1 \
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

## Relationship to RNA-MAP-Mini

This optimization toolkit depends on [rna-map-mini](https://github.com/jyesselm/rna-map-mini):

- **rna-map-mini**: Core Python package for RNA-MAP analysis
- **This repo**: Parameter optimization and research tools

The optimization code uses the `rna-map-mini` Python package to execute pipelines and evaluate results.

## Contributing

Contributions are welcome! Please see the main repository for contribution guidelines.

## License

Apache License 2.0 - See [LICENSE](./LICENSE) file for details.

## Citation

If you use this optimization toolkit, please cite the rna-map-mini repository.

## Support

For issues and questions:
- Optimization-specific: Open an issue in this repository
- Core analysis: Open an issue in [rna-map-mini](https://github.com/jyesselm/rna-map-mini)
