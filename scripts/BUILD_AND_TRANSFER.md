# Build Container Locally and Transfer to Cluster

## Quick Start

### Step 1: Start Docker Desktop

Make sure Docker Desktop is running on your local machine:
- **macOS**: Open Docker Desktop from Applications
- **Linux**: `sudo systemctl start docker` (or your distro's equivalent)
- **Windows**: Open Docker Desktop

Verify Docker is running:
```bash
docker ps
```

### Step 2: Build Container

```bash
cd /path/to/rna_map_nextflow
bash scripts/optimization/build_with_docker.sh rna-map-optimization.sif
```

This will:
1. Build Docker image (takes ~10-20 minutes)
2. Convert to Apptainer format (creates `.sif` file)

### Step 3: Transfer to Cluster

```bash
# Replace with your cluster details
scp rna-map-optimization.sif username@cluster.hostname:/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/
```

### Step 4: Use on Cluster

```bash
# SSH to cluster
ssh username@cluster.hostname

# Set container path
export CONTAINER_PATH=/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/rna-map-optimization.sif

# Run optimization
bash scripts/optimization/submit_optimization_jobs.sh
```

## Troubleshooting

### Docker daemon not running
- Start Docker Desktop application
- Wait for it to fully start (whale icon in menu bar)
- Run `docker ps` to verify

### Build fails
- Check you have enough disk space (~3-4 GB)
- Ensure network connection is stable
- Check Docker has enough resources allocated

### Transfer fails
- Verify cluster path exists
- Check you have write permissions on cluster
- Ensure SSH key is set up correctly

