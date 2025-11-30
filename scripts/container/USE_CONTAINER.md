# Using the Optimization Container

Once you have built the `rna-map-optimization.sif` container on the cluster, follow these steps to use it.

## Step 1: Set Container Path

Set the environment variable to point to your container:

```bash
export CONTAINER_PATH=/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/rna-map-optimization.sif
```

**Or add to your shell profile** (`.bashrc` or `.zshrc`) to make it persistent:

```bash
echo 'export CONTAINER_PATH=/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/rna-map-optimization.sif' >> ~/.bashrc
source ~/.bashrc
```

## Step 2: Verify Container Works

Test that the container is accessible and imports work:

```bash
# Test rna_map import
apptainer exec $CONTAINER_PATH python -c "import rna_map; print('✓ rna_map installed')"

# Test Optuna import
apptainer exec $CONTAINER_PATH python -c "import optuna; print('✓ optuna installed')"

# Test Bowtie2 is available
apptainer exec $CONTAINER_PATH bowtie2 --version
```

## Step 3: Prepare Your Data

Organize your test cases in a `data/` directory:

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

## Step 4: Configure Optimization Settings (Optional)

Edit `scripts/optimization/cluster_optimization_config.yml` if needed:

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

## Step 5: Submit Optimization Jobs

```bash
# Make sure CONTAINER_PATH is set
export CONTAINER_PATH=/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/rna-map-optimization.sif

# Submit jobs
bash scripts/optimization/submit_optimization_jobs.sh
```

The script will:
- Detect the container automatically
- Generate SLURM job scripts that use the container
- Mount data and results directories
- Submit jobs to the cluster

## Step 6: Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View job output
tail -f optimization_results/case_1/<job_id>.out

# View job errors
tail -f optimization_results/case_1/<job_id>.err
```

## Step 7: Collect Results

After all jobs complete:

```bash
bash scripts/optimization/collect_all_results.sh
```

This will:
- Collect top results from each case
- Aggregate results
- Create summary files in `optimization_results/aggregated/`

## Troubleshooting

### Container not found
```bash
# Verify container exists
ls -lh $CONTAINER_PATH

# Check CONTAINER_PATH is set
echo $CONTAINER_PATH
```

### Import errors in container
```bash
# Test imports
apptainer exec $CONTAINER_PATH python -c "import rna_map; import optuna; print('OK')"
```

### Jobs failing
- Check `.err` files for error messages
- Verify data files exist and are readable
- Check container path in job scripts

### Container path not set in jobs
Make sure `CONTAINER_PATH` is exported before running `submit_optimization_jobs.sh`:
```bash
export CONTAINER_PATH=/path/to/rna-map-optimization.sif
bash scripts/optimization/submit_optimization_jobs.sh
```

