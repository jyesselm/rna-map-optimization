# Building Container on Cluster (No Docker Required)

Since you don't have Docker on the cluster, the best approach is to build the Apptainer container directly on the cluster using the definition file.

## Method 1: Build from Definition File on Cluster (Recommended)

### Step 1: Transfer Required Files to Cluster

```bash
# From your local machine (from project root)
scp scripts/container/rna-map-optimization.def user@cluster:/path/on/cluster/
scp environment.yml user@cluster:/path/on/cluster/
scp -r scripts/ user@cluster:/path/on/cluster/
```

### Step 2: Build on Cluster

```bash
# SSH to cluster
ssh user@cluster

# Navigate to directory with files
cd /path/on/cluster

# Build container (requires fakeroot or sudo)
apptainer build rna-map-optimization.sif rna-map-optimization.def
```

**If you have fakeroot capability:**
```bash
apptainer build --fakeroot rna-map-optimization.sif scripts/optimization/rna-map-optimization.def
```

**If you need sudo:**
```bash
sudo apptainer build rna-map-optimization.sif scripts/optimization/rna-map-optimization.def
```

### Step 3: Use Container

```bash
export CONTAINER_PATH=/path/on/cluster/rna-map-optimization.sif
bash scripts/cluster/submit_optimization_jobs.sh
```

## Method 2: Use Conda Environment (No Container)

If building the container is problematic, use conda environment instead:

```bash
# On cluster - create environment
conda env create -f environment.yml
conda activate rna-map-optimization

# Then run optimization (no container needed)
bash scripts/cluster/submit_optimization_jobs.sh
```

## Method 3: Transfer Docker Image and Convert (If Docker Available)

If Docker becomes available on cluster later:

```bash
# On local machine - save Docker image
docker save rna-map-optimization:latest | gzip > rna-map-optimization.tar.gz

# Transfer to cluster
scp rna-map-optimization.tar.gz user@cluster:/path/

# On cluster - load and convert
docker load < rna-map-optimization.tar.gz
apptainer build rna-map-optimization.sif docker-daemon://rna-map-optimization:latest
```

## Quick Transfer Script

Create a script to transfer all needed files:

```bash
#!/bin/bash
# transfer_build_files.sh

CLUSTER_USER="your_username"
CLUSTER_HOST="cluster.hostname"
CLUSTER_PATH="/work/yesselmanlab/jyesselm/installs/rna_map_nextflow"

echo "Transferring files to cluster..."
scp scripts/container/rna-map-optimization.def ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_PATH}/
scp environment.yml ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_PATH}/
scp -r scripts/ ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_PATH}/

echo "✅ Files transferred!"
echo ""
echo "On cluster, run:"
echo "  cd ${CLUSTER_PATH}"
echo "  apptainer build --fakeroot rna-map-optimization.sif rna-map-optimization.def"
```

## Troubleshooting

### Fakeroot not available
Request from cluster admin:
```bash
sudo apptainer capability add --user $USER fakeroot
```

### Build fails
- Check that all files transferred correctly
- Verify paths in definition file match cluster structure
- Check disk space (container needs ~3-4 GB)

### Can't build on cluster
Use conda environment instead (Method 2) - no container needed!

