#!/bin/bash
# Minimal example: Run optimization with Apptainer container
#
# Usage: Update paths and run: bash example_container_minimal.sh

CONTAINER="/path/to/rna-map-optimization.sif"
CASE_DIR="test_cases/case_1"
OUTPUT_DIR="results"

# Use apptainer if available, otherwise singularity
CONTAINER_CMD=$(command -v apptainer || command -v singularity)

# Run with container
${CONTAINER_CMD} exec \
    -B ${CASE_DIR}:/data:ro \
    -B ${OUTPUT_DIR}:/results \
    ${CONTAINER} \
    rna-map-optimize optimize \
        --case-dir /data \
        --output-dir /results \
        --n-trials 100

