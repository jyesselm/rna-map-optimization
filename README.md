# RNA-MAP Optimization Toolkit

Parameter optimization toolkit for RNA-MAP experiments using Bowtie2 alignment.

## Installation

```bash
# Clone repository
git clone https://github.com/jyesselm/rna-map-optimization.git
cd rna-map-optimization

# Create conda environment
conda env create -f environment.yml
conda activate rna-map-optimization

# Install package
pip install -e .
```

## Quick Start

### Run Optimization

```bash
rna-map-optimize optimize \
    --case-dir test_cases/case_1 \
    --n-trials 100 \
    --output-dir results
```

The tool automatically detects FASTA and FASTQ files in the case directory.

### Analyze Results

```bash
# Collect top results
rna-map-optimize collect-results \
    --results-dir results \
    --top-n 100

# Analyze parameter importance
rna-map-optimize analyze-parameters \
    --results-dir results/case_1 \
    --top-n 100
```

## Building Container for Cluster

### Build SIF File on Cluster

1. **Load Apptainer/Singularity module:**
   ```bash
   module load apptainer  # or: module load singularity
   ```

2. **Build the container:**
   ```bash
   # From the project root directory
   apptainer build --fakeroot rna-map-optimization.sif containers/rna-map-optimization.def
   ```
   
   This will create `rna-map-optimization.sif` in the current directory. The build takes 10-20 minutes.
   
   **Note:** If your cluster doesn't support `--fakeroot`, you may need to use `sudo` or build on a system with root access. The definition file is located at `containers/rna-map-optimization.def`.

3. **Move to accessible location:**
   ```bash
   # Move to a shared location accessible to compute nodes
   mv rna-map-optimization.sif /path/to/shared/storage/
   ```

### Test Container

1. **Test basic functionality:**
   ```bash
   apptainer exec rna-map-optimization.sif rna-map-optimize --help
   ```

2. **Test with a small case:**
   ```bash
   # Create test directories
   mkdir -p test_container/{input,output}
   cp test_cases/case_1/* test_container/input/
   
   # Run test (5 trials for quick test)
   apptainer exec \
       -B test_container/input:/data:ro \
       -B test_container/output:/results \
       rna-map-optimization.sif \
       rna-map-optimize optimize \
           --case-dir /data \
           --output-dir /results \
           --n-trials 5
   
   # Check results
   ls test_container/output/
   ```

3. **Verify output files exist:**
   ```bash
   # Should see these files:
   ls -lh test_container/output/
   # - optuna_study.json          (best parameters)
   # - optuna_summary.csv         (all trial results)
   # - final_bit_vector_metrics.json  (mutation stats)
   # - visualizations/            (HTML plots)
   ```

4. **Verify container has all dependencies:**
   ```bash
   apptainer exec rna-map-optimization.sif python -c "import optuna; print('✓ Optuna OK')"
   apptainer exec rna-map-optimization.sif python -c "import rna_map_mini; print('✓ rna-map-mini OK')"
   apptainer exec rna-map-optimization.sif bowtie2 --version
   ```

## Running on Cluster with SLURM

### Example SLURM Script

Create `run_optimization.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=bt2_optimize
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=optimization_%j.out
#SBATCH --error=optimization_%j.err

# Update these paths
CONTAINER="/path/to/rna-map-optimization.sif"
CASE_DIR="/path/to/test_cases/case_1"
OUTPUT_DIR="/path/to/results/case_1"

# Use apptainer if available, otherwise singularity
CONTAINER_CMD=$(command -v apptainer || command -v singularity)

# Run optimization
${CONTAINER_CMD} exec \
    -B ${CASE_DIR}:/data:ro \
    -B ${OUTPUT_DIR}:/results \
    ${CONTAINER} \
    rna-map-optimize optimize \
        --case-dir /data \
        --output-dir /results \
        --n-trials 200 \
        --threads 8
```

Submit job:
```bash
sbatch run_optimization.sh
```

See `scripts/example_sbatch_container.sh` for a complete example.

## Configuration

Parameter search ranges are configured in `config/optimization_config.yml`. Edit this file to adjust:
- Parameter ranges (min/max values)
- Categorical options
- Search space

## Output Files

After optimization completes, you'll find in the output directory:

- `optuna_study.json` - Best parameters and study metadata
- `optuna_summary.csv` - All trial results
- `final_bit_vector_metrics.json` - Detailed mutation statistics
- `visualizations/` - HTML plots of optimization history
- `index/` - Bowtie2 index files (reusable)

The `results/` directory (containing trial SAM files) is automatically removed after optimization to save space.

## Command Reference

### Optimize
```bash
rna-map-optimize optimize \
    --case-dir DIR          # Case directory (auto-detects files)
    --cases-dir DIR         # Batch mode (multiple cases)
    --n-trials N            # Number of trials (default: 100)
    --output-dir DIR        # Output directory
    --threads N             # Threads for alignment
    --mapq-cutoff N         # MAPQ cutoff (default: 20)
    --keep-intermediates    # Keep trial files (default: cleanup)
```

### Analyze Sequences
```bash
rna-map-optimize analyze-sequences \
    --fasta FILE           # Input FASTA
    --output FILE          # Output FASTA (optional)
    --min-common-length N  # Min common sequence length
```

### Collect Results
```bash
rna-map-optimize collect-results \
    --results-dir DIR      # Results directory
    --top-n N              # Top N results per case
    --output-dir DIR       # Aggregated output
```

### Analyze Parameters
```bash
rna-map-optimize analyze-parameters \
    --results-dir DIR      # Results directory
    --top-n N              # Top N trials to analyze
    --output FILE           # JSON output (optional)
```

## Directory Structure

```
rna-map-optimization/
├── src/rna_map_optimization/  # Main package
├── containers/                # Container build files
│   ├── rna-map-optimization.def  # Apptainer definition
│   └── Dockerfile            # Docker build file
├── scripts/                   # Helper scripts
├── config/                    # Configuration files
└── test_cases/               # Test data
```

## License

Apache License 2.0
