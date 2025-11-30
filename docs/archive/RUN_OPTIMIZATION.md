# How to Run Bowtie2 Parameter Optimization

This guide explains how to run the improved Optuna-based optimization script.

## Prerequisites

### 1. Set up the conda environment (one-time setup)

```bash
# Make sure you're in the project directory
cd /Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/nextflow/rna_map_netflow

# Run the setup script (creates environment and installs dependencies)
./setup_optuna_env.sh

# Or manually:
conda env create -f environment_optuna.yml
conda activate rna-map-optuna
cd src/rna_map
python -m pip install -e .
cd ../..
```

### 2. Verify the environment

```bash
conda activate rna-map-optuna
python test_optuna_env.py
```

All tests should pass (✓).

## Running Optimization

### Basic Usage (Recommended for testing)

Start with a small number of trials to test:

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 50 \
    --threads 4 \
    --output-dir my_optimization
```

### Full Optimization (Production run)

For a comprehensive search, use more trials:

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 200 \
    --threads 8 \
    --output-dir optimization_results
```

### Single-end Reads

If you have single-end reads (no mate2 file):

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta reference.fasta \
    --fastq1 reads.fastq \
    --n-trials 100 \
    --threads 4 \
    --output-dir optimization_results
```

### With Thread Optimization

If you want to optimize thread count as well (usually not necessary):

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 100 \
    --threads 8 \
    --optimize-threads \
    --output-dir optimization_results
```

### With Custom Settings

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --n-trials 150 \
    --threads 6 \
    --mapq-cutoff 25 \
    --read-length 150 \
    --output-dir my_optimization \
    --study-name "my_study"
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--fasta` | Path to reference FASTA file | **Required** |
| `--fastq1` | Path to first FASTQ file | **Required** |
| `--fastq2` | Path to second FASTQ file (paired-end) | Optional |
| `--output-dir` | Output directory for results | `bowtie2_optimization_optuna` |
| `--n-trials` | Number of Optuna trials | `100` |
| `--threads` | Number of threads for alignment | `4` |
| `--optimize-threads` | Optimize thread count | `False` |
| `--mapq-cutoff` | MAPQ cutoff for filtering | `20` |
| `--read-length` | Expected read length (bp) | `150` |
| `--timeout` | Timeout in seconds | `None` |
| `--study-name` | Name for Optuna study | `bowtie2_optimization` |
| `--storage` | Optuna storage URL (for resuming) | `None` |

## What Gets Optimized

The script optimizes these parameters:

1. **seed_length**: 10-22 (step 2)
2. **seed_mismatches**: 0-1
3. **maxins**: 200-1200 (step 100)
4. **score_min**: 20+ options (L-type, G-type, S-type)
5. **mismatch_penalty**: 4 options
6. **gap_penalty_read**: 4 options
7. **gap_penalty_ref**: 4 options
8. **sensitivity_mode**: 5 options
9. **threads**: 2-16 (if `--optimize-threads` is used)

## Output Files

After running, you'll find these files in the output directory:

```
optimization_results/
├── optuna_study.json          # Summary with best parameters
├── optuna_summary.csv         # All trial results (spreadsheet-friendly)
├── optuna_study.pkl          # Pickled study (for resuming/analysis)
├── index/                     # Bowtie2 index files
├── results/                   # Individual trial results
│   ├── trial_0/
│   ├── trial_1/
│   └── ...
└── visualizations/            # HTML visualization files
    ├── optimization_history.html
    ├── param_importances.html
    └── parallel_coordinate.html
```

## Viewing Results

### Quick summary:

```bash
# View best parameters
cat optimization_results/optuna_study.json | python -m json.tool

# View all trials in CSV format
cat optimization_results/optuna_summary.csv
```

### Open visualizations:

```bash
# Open in browser (macOS)
open optimization_results/visualizations/optimization_history.html
open optimization_results/visualizations/param_importances.html
```

### View best Bowtie2 arguments:

The best Bowtie2 arguments are printed at the end of the run and saved in `optuna_study.json`:

```bash
grep "best_bowtie2_args" optimization_results/optuna_study.json
```

## Using Best Parameters

After optimization, use the best parameters in your Nextflow pipeline:

```bash
nextflow run main.nf \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --bt2_alignment_args "--local --no-unal --no-discordant --no-mixed -L 16 -X 1000 --score-min S,20,0.15 --mp 6,2 --rdg 8,4 --rfg 5,3"
```

## Tips

1. **Start small**: Run 20-50 trials first to test
2. **Review early results**: Check if signal-to-noise is improving
3. **Use more trials**: 100-200 trials for comprehensive search
4. **Parallel runs**: Run multiple optimizations in parallel on different datasets
5. **Resume studies**: Use `--storage sqlite:///study.db` to save and resume
6. **Review CSV**: Use Excel or pandas to analyze trends in `optuna_summary.csv`

## Troubleshooting

### Environment issues:
```bash
# Check environment is activated
conda activate rna-map-optuna
which python  # Should point to conda environment

# Reinstall if needed
conda env remove -n rna-map-optuna
./setup_optuna_env.sh
```

### Bowtie2 not found:
```bash
conda activate rna-map-optuna
conda install -c bioconda bowtie2
```

### Out of memory:
- Reduce `--threads`
- Reduce `--n-trials`
- Process smaller datasets

### Slow optimization:
- Use fewer trials (`--n-trials 50`)
- Increase threads (`--threads 8`)
- Check disk space (temporary SAM files created)

## Example Workflow

```bash
# 1. Setup (one-time)
./setup_optuna_env.sh
conda activate rna-map-optuna

# 2. Quick test (5-10 minutes)
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 20 \
    --output-dir test_optimization

# 3. Review results
cat test_optimization/optuna_study.json | python -m json.tool
open test_optimization/visualizations/optimization_history.html

# 4. Full optimization (1-3 hours)
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 200 \
    --threads 8 \
    --output-dir full_optimization

# 5. Use best parameters in pipeline
# Copy the best_bowtie2_args from optuna_study.json
```

## Quick Reference

**Minimum command to get started:**

```bash
conda activate rna-map-optuna && \
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 50 \
    --output-dir my_optimization
```

That's it! The script will:
1. Build the Bowtie2 index
2. Run optimization trials
3. Generate results and visualizations
4. Print the best parameters

