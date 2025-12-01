# Cluster-Based Bowtie2 Optimization

This directory contains scripts for running large-scale Bowtie2 parameter optimization on a cluster.

## Setup

### Option 1: Using Container (Recommended)

#### Option 1a: Build Locally with Docker (Easiest - No sudo needed)

Build on your local machine with Docker, then transfer to cluster:

```bash
# On your local machine (with Docker installed)
bash scripts/optimization/build_with_docker.sh rna-map-optimization.sif

# Transfer to cluster
scp rna-map-optimization.sif user@cluster:/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/
```

This will:
- Build Docker image locally (no sudo needed for Docker)
- Automatically convert to Apptainer format
- Create a `.sif` file ready to transfer

#### Option 1b: Build on Cluster

Build an Apptainer/Singularity container with all dependencies:

```bash
bash scripts/optimization/build_optimization_container.sh rna-map-optimization.sif
```

This will:
- Create a container image with all dependencies
- Include Bowtie2, Python, Optuna, and rna-map-mini package
- No need to manage conda environments on compute nodes
- Results are written to mounted directories

**No sudo required**: The build script automatically uses fakeroot when available. If you don't have sudo access:
1. **Request fakeroot capability** from your cluster admin: `sudo apptainer capability add --user $USER fakeroot`
2. **Build on another system** with sudo access, then transfer the `.sif` file to the cluster
3. **Use Docker** (if available) to build, then convert to Apptainer format

**Usage with container:**
```bash
export CONTAINER_PATH=/path/to/rna-map-optimization.sif
bash scripts/optimization/submit_optimization_jobs.sh
```

### Option 2: Using Conda Environment

Alternatively, create a conda environment:

```bash
bash scripts/optimization/setup_optimization_env.sh
```

This will:
- Create a conda environment named `rna-map-optimization`
- Install Bowtie2, Python dependencies, and Optuna
- Installs `rna-map-mini` package from GitHub via pip (included in environment.yml)

### 2. Prepare Test Cases

Create a `data/` directory in the project root and organize test cases:

```
data/
├── case_1/
│   ├── reference.fasta
│   ├── reads_R1.fastq.gz
│   └── reads_R2.fastq.gz
├── case_2/
│   ├── reference.fasta
│   ├── reads_R1.fastq.gz
│   └── reads_R2.fastq.gz
└── ...
```

**Required files per case:**
- `*.fasta` or `*.fa` - Reference sequence(s)
- `*_R1*.fastq*` or `*_1.fastq*` - Read 1 (required)
- `*_R2*.fastq*` or `*_2.fastq*` - Read 2 (optional, for paired-end)

### 3. Configure Optimization Settings

Edit `scripts/optimization/cluster_optimization_config.yml` to set:

- **Optimization parameters**: Number of combinations, threads, cutoffs
- **Cluster settings**: Time limits, memory, CPUs
- **Output settings**: Base directory, whether to keep intermediates

Example configuration:

```yaml
optimization:
  max_combinations: 200
  read_length: 150
  threads: 8
  mapq_cutoff: 20

cluster:
  time: "24:00:00"
  memory: "16G"
  cpus: 8
```

### 4. Submit Jobs

Run the submission script:

```bash
bash scripts/optimization/submit_optimization_jobs.sh
```

This will:
- Scan `data/` for all test cases
- Create SLURM job scripts for each case
- Submit jobs to the cluster
- Create output directories for results

## Output Structure

Results will be organized as:

```
optimization_results/
├── case_1/
│   ├── <job_id>.out          # Job stdout
│   ├── <job_id>.err          # Job stderr
│   ├── index/                 # Bowtie2 index
│   ├── results/               # Alignment results
│   │   ├── baseline/          # Baseline results
│   │   └── combo_XXX/         # Parameter combination results
│   ├── best_parameters.json   # Best parameters found
│   ├── optimization_results.json  # All results
│   └── optimization_summary.csv   # Summary CSV
├── case_2/
│   └── ...
└── job_scripts/               # Submitted job scripts
    ├── case_1_optimization.sh
    └── case_2_optimization.sh
```

## Monitoring Jobs

Check job status:

```bash
squeue -u $USER
```

View job output:

```bash
tail -f optimization_results/case_1/<job_id>.out
```

Cancel a job:

```bash
scancel <job_id>
```

## Collecting Results

After all jobs complete, collect and aggregate the top results:

### Quick Method (Recommended)

```bash
bash scripts/optimization/collect_all_results.sh
```

This automatically:
- Activates the conda environment
- Collects top 100 from each case
- Saves aggregated results to `optimization_results/aggregated/`

### Manual Method

```bash
# Activate environment
conda activate rna-map-optimization

# Collect top 100 from each case
rna-map-optimize collect-results \
    --results-dir optimization_results \
    --top-n 100 \
    --output-dir optimization_results/aggregated
```

### Output Files

The collection script creates:
- `aggregated_top_results.json` - Full aggregated results with all cases
- `combined_top_results.csv` - All top results from all cases combined
- `top_100_overall.csv` - Top 100 combinations across all cases (sorted by quality_score)
- `summary_statistics.json` - Summary statistics and best per case

### Options

- `--top-n`: Number of top results per case (default: 100)
- `--min-quality`: Minimum quality score threshold (optional, filters results)
- `--output-dir`: Where to save aggregated results (default: results-dir/aggregated)
- `--results-dir`: Base directory with optimization results (default: optimization_results)

### Example: Filter by Quality

```bash
# Only include results with quality_score >= 0.8
rna-map-optimize collect-results \
    --results-dir optimization_results \
    --top-n 100 \
    --min-quality 0.8
```

## Per-Case Overrides

You can override settings for specific cases in the config file:

```yaml
case_overrides:
  case_1:
    max_combinations: 500
    time: "48:00:00"
    memory: "32G"
  case_2:
    max_combinations: 100
    time: "12:00:00"
```

## Single vs Multiple Sequences

The optimization script automatically detects the number of sequences:

- **Single sequence**: Uses permissive parameters for mutation detection
- **Multiple sequences**: Uses stricter parameters for better discrimination

The script will:
- Use `-k 1` to only report best alignments
- Adjust seed length, scoring thresholds, and effort parameters
- Track multi-mapping statistics

## Troubleshooting

### Environment not found
```bash
# Recreate environment
bash scripts/optimization/setup_optimization_env.sh
```

### Jobs failing
- Check `.err` files for error messages
- Verify all input files exist and are readable
- Check that Bowtie2 is available in the environment
- Ensure sufficient time/memory allocated

### No test cases found
- Verify `data/` directory exists
- Check that subdirectories contain `.fasta` and `.fastq` files
- Ensure file naming matches expected patterns

## Example Workflow

```bash
# 1. Setup environment (once)
bash scripts/optimization/setup_optimization_env.sh

# 2. Prepare data
mkdir -p data/case_1
cp reference.fasta data/case_1/
cp reads_R1.fastq.gz data/case_1/
cp reads_R2.fastq.gz data/case_1/

# 3. Configure (edit config file)
vim scripts/optimization/cluster_optimization_config.yml

# 4. Submit jobs
bash scripts/optimization/submit_optimization_jobs.sh

# 5. Monitor
squeue -u $USER

# 6. Check results
ls optimization_results/case_1/
cat optimization_results/case_1/best_parameters.json
```

