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
