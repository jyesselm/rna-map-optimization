# Optuna Environment Setup Guide

This guide explains how to set up and use the conda environment for Optuna-based Bowtie2 parameter optimization.

## Quick Setup

### Option 1: Using the Setup Script (Recommended)

```bash
# Make the setup script executable (if not already)
chmod +x setup_optuna_env.sh

# Run the setup script
./setup_optuna_env.sh
```

This script will:
1. Create a conda environment named `rna-map-optuna`
2. Install all required dependencies (Optuna, Plotly, Bowtie2, etc.)
3. Install the `rna_map` package in development mode

### Option 2: Manual Setup

```bash
# Create the environment from the YAML file
conda env create -f environment_optuna.yml

# Activate the environment
conda activate rna-map-optuna

# Install the rna_map package in development mode
cd src/rna_map
python -m pip install -e .
cd ../..
```

## Verifying the Installation

Run the test script to verify everything is set up correctly:

```bash
conda activate rna-map-optuna
python test_optuna_env.py
```

You should see all tests pass (✓).

## Running Optimization

Once the environment is set up, you can run the optimization script:

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 100 \
    --output-dir optimization_optuna
```

### Command-Line Options

- `--fasta`: Path to reference FASTA file (required)
- `--fastq1`: Path to first FASTQ file (required)
- `--fastq2`: Path to second FASTQ file (optional, for paired-end)
- `--output-dir`: Output directory for results (default: `bowtie2_optimization_optuna`)
- `--n-trials`: Number of Optuna trials to run (default: 100)
- `--threads`: Number of threads for alignment (default: 4)
- `--mapq-cutoff`: MAPQ cutoff for high-quality alignments (default: 20)
- `--timeout`: Timeout in seconds (optional)
- `--study-name`: Name for Optuna study (default: `bowtie2_optimization`)
- `--storage`: Optuna storage URL (e.g., `sqlite:///study.db`) for resuming studies

## Environment Contents

The `rna-map-optuna` environment includes:

- **Python**: >=3.10
- **Core dependencies**:
  - pandas >=1.5
  - numpy >=1.21
  - tabulate >=0.9
  - pyyaml >=6.0
- **Optuna dependencies**:
  - optuna >=3.0
  - plotly >=5.0 (for visualizations)
- **Bioinformatics tools**:
  - bowtie2 (for read alignment)
- **Local package**:
  - rna_map (installed in development mode)

## Output Files

After running optimization, the output directory will contain:

- `optuna_study.json`: Summary of the best trial and parameters
- `optuna_summary.csv`: CSV file with all trial results
- `optuna_study.pkl`: Pickled Optuna study object (for resuming)
- `visualizations/`: HTML visualization files
  - `optimization_history.html`: Optimization progress over time
  - `param_importances.html`: Parameter importance analysis
  - `parallel_coordinate.html`: Parallel coordinate plot
- `results/trial_N/`: Individual trial results (SAM files, etc.)
- `index/`: Bowtie2 index files

## Troubleshooting

### Module Not Found Errors

If you get import errors, make sure:
1. The environment is activated: `conda activate rna-map-optuna`
2. The package is installed: `cd src/rna_map && python -m pip install -e . && cd ../..`
3. You're using the correct Python: `which python` should point to the conda environment

### Bowtie2 Not Found

If bowtie2 is not found:
```bash
conda activate rna-map-optuna
conda install -c bioconda bowtie2
```

### Updating Dependencies

To update the environment:
```bash
conda activate rna-map-optuna
conda env update -f environment_optuna.yml --prune
```

## Removing the Environment

To remove the environment:

```bash
conda env remove -n rna-map-optuna
```

## Example Workflow

```bash
# 1. Setup environment (one-time)
./setup_optuna_env.sh

# 2. Activate environment
conda activate rna-map-optuna

# 3. Verify installation
python test_optuna_env.py

# 4. Run optimization
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 50 \
    --output-dir my_optimization

# 5. View results
cat my_optimization/optuna_study.json
open my_optimization/visualizations/optimization_history.html
```

