# Using Container for Optimization

This guide explains how to use the Apptainer/Singularity container for running Bowtie2 parameter optimization.

## Why Use a Container?

1. **No Environment Setup**: All dependencies are pre-installed
2. **Consistent Results**: Same environment on all compute nodes
3. **Faster Execution**: No package installation during jobs
4. **Portability**: Works across different clusters
5. **Reproducibility**: Exact same environment every time

## Building the Container

### Prerequisites

- **For local Docker build**: Docker installed on your local machine
- **For cluster build**: Apptainer or Singularity installed on cluster
- Sufficient disk space (~3-4 GB)

### Build Process

#### Method 1: Build Locally with Docker (Recommended)

Build on your local machine, then transfer to cluster:

```bash
# On your local machine (with Docker)
# Build Docker image
docker build -f scripts/container/Dockerfile -t rna-map-optimization:latest ..

# Build Apptainer/Singularity from Docker
apptainer build rna-map-optimization.sif docker-daemon://rna-map-optimization:latest

# Or build directly with Apptainer
apptainer build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def

# Transfer to cluster
scp rna-map-optimization.sif user@cluster:/path/to/destination/
```

**Advantages**:
- No sudo required (Docker works without root on most systems)
- No need for fakeroot capability
- Can build anywhere Docker is available

#### Method 2: Build on Cluster

Build directly on the cluster:

```bash
# From project root on cluster
bash scripts/optimization/build_optimization_container.sh rna-map-optimization.sif
```

**Requires**: Fakeroot capability or sudo access

The build process:
1. Creates conda environment from `environment.yml`
2. Installs Bowtie2, Python dependencies, and Optuna
3. Installs `rna-map-mini` package from GitHub (via pip in environment.yml)
4. Packages everything into a single `.sif` file

**Build time**: ~10-20 minutes (depending on network speed)

**No sudo required**: The build script automatically uses fakeroot if available. If you don't have sudo access:

1. **Request fakeroot capability** from your cluster admin:
   ```bash
   sudo apptainer capability add --user $USER fakeroot
   ```
   Then the build script will work without sudo.

2. **Build on another system** with sudo access, then transfer the `.sif` file:
   ```bash
   # On a system with sudo (e.g., your local machine)
   bash scripts/optimization/build_optimization_container.sh rna-map-optimization.sif
   # Copy to cluster
   scp rna-map-optimization.sif user@cluster:/path/to/destination/
   ```

3. **Use Docker** (if available) to build, then convert:
   ```bash
   # Build Docker image
   docker build -f docker/Dockerfile -t rna-map-optimization .
   # Convert to Apptainer (still needs fakeroot or sudo)
   apptainer build rna-map-optimization.sif docker-daemon://rna-map-optimization:latest
   ```

### Verify Container

```bash
# Test that rna-map-mini is installed
apptainer exec rna-map-optimization.sif python -c "import rna_map_mini; print('✓ rna-map-mini installed')"

# Test that Bowtie2 is available
apptainer exec rna-map-optimization.sif bowtie2 --version

# Test that Optuna is available
apptainer exec rna-map-optimization.sif python -c "import optuna; print('✓ optuna installed')"
```

## Using the Container

### Method 1: Automatic (Recommended)

Set the container path and submit jobs normally:

```bash
export CONTAINER_PATH=/path/to/rna-map-optimization.sif
bash scripts/optimization/submit_optimization_jobs.sh
```

The submission script will automatically:
- Detect the container
- Generate SLURM job scripts that use the container
- Mount data and results directories
- Run optimization inside the container

### Method 2: Manual Execution

Run optimization directly with the container:

```bash
apptainer exec \
    -B /path/to/data:/data \
    -B /path/to/results:/results \
    -B /path/to/project:/work \
    rna-map-optimization.sif \
    python /work/scripts/optimize_bowtie2_params.py \
        --fasta /data/case_1/reference.fasta \
        --fastq1 /data/case_1/reads_R1.fastq.gz \
        --fastq2 /data/case_1/reads_R2.fastq.gz \
        --output-dir /results/case_1 \
        --max-combinations 200
```

### Directory Mounting

The container uses these mount points:
- `/data` - Input data directory (FASTA, FASTQ files)
- `/results` - Output directory for results
- `/work` - Project root (for scripts)

**Important**: Paths inside the container are different from host paths!

## Container Contents

The container includes:
- **Python 3.10+** with conda
- **Bowtie2** aligner
- **Optuna** for Bayesian optimization
- **rna-map-mini** package (installed from GitHub)
- **All Python dependencies** (pandas, numpy, pyyaml, tabulate, plotly)

## Job Script Example

When using the container, job scripts look like this:

```bash
#!/bin/bash
#SBATCH --job-name=bt2_opt_case1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

CONTAINER_PATH="/path/to/rna-map-optimization.sif"
BIND_ARGS="-B /data:/data -B /results:/results -B /work:/work"

apptainer exec ${BIND_ARGS} ${CONTAINER_PATH} \
    python /work/scripts/optimize_bowtie2_params.py \
        --fasta /data/case_1/reference.fasta \
        --fastq1 /data/case_1/reads_R1.fastq.gz \
        --output-dir /results/case_1
```

## Advantages Over Conda

| Feature | Container | Conda |
|---------|-----------|-------|
| Setup time | Once (build) | Per node |
| Consistency | Guaranteed | Depends on env |
| Disk usage | Single file | Multiple files |
| Portability | High | Medium |
| Reproducibility | Excellent | Good |

## Troubleshooting

### Container not found
```bash
# Check container exists
ls -lh rna-map-optimization.sif

# Verify path in CONTAINER_PATH
echo $CONTAINER_PATH
```

### Permission errors
```bash
# Make sure container is readable
chmod 644 rna-map-optimization.sif

# Check Apptainer/Singularity permissions
apptainer exec rna-map-optimization.sif whoami
```

### Mount point issues
```bash
# Verify directories exist and are accessible
ls -ld /path/to/data /path/to/results

# Test mounting
apptainer exec -B /path/to/data:/data rna-map-optimization.sif ls /data
```

### Import errors
```bash
# Test imports inside container
apptainer exec rna-map-optimization.sif python -c "import rna_map_mini; import optuna; print('OK')"
```

## Best Practices

1. **Store container on shared storage** - Accessible to all compute nodes
2. **Use absolute paths** - Avoid path resolution issues
3. **Mount scratch space** - Use `$SCRATCH` for temporary files
4. **Version containers** - Tag containers with version numbers
5. **Test before submitting** - Verify container works with a small test case

## Example Workflow

```bash
# 1. Build container (once)
bash scripts/optimization/build_optimization_container.sh \
    /shared/containers/rna-map-optimization-v1.0.sif

# 2. Set container path
export CONTAINER_PATH=/shared/containers/rna-map-optimization-v1.0.sif

# 3. Submit jobs
bash scripts/optimization/submit_optimization_jobs.sh

# 4. Monitor jobs
squeue -u $USER

# 5. Collect results (same as before)
bash scripts/optimization/collect_all_results.sh
```

## Updating the Container

To update the container with new code:

```bash
# Rebuild container
bash scripts/optimization/build_optimization_container.sh \
    /shared/containers/rna-map-optimization-v1.1.sif

# Update CONTAINER_PATH
export CONTAINER_PATH=/shared/containers/rna-map-optimization-v1.1.sif

# Resubmit jobs
bash scripts/optimization/submit_optimization_jobs.sh
```

