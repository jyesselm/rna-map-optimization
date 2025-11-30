#!/bin/bash
# Transfer files needed to build container on cluster
#
# Usage:
#   ./scripts/optimization/transfer_build_files.sh user@cluster /path/on/cluster
#
# Example:
#   ./scripts/optimization/transfer_build_files.sh jyesselm@cluster.hostname /work/yesselmanlab/jyesselm/installs/rna_map_nextflow

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 user@cluster /path/on/cluster"
    echo ""
    echo "Example:"
    echo "  $0 jyesselm@cluster.hostname /work/yesselmanlab/jyesselm/installs/rna_map_nextflow"
    exit 1
fi

CLUSTER="$1"
CLUSTER_PATH="$2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

echo "=========================================="
echo "Transferring Build Files to Cluster"
echo "=========================================="
echo ""
echo "Cluster: $CLUSTER"
echo "Path: $CLUSTER_PATH"
echo ""

# Create remote directory structure
echo "Creating remote directory structure..."
ssh "$CLUSTER" "mkdir -p ${CLUSTER_PATH}/scripts/optimization ${CLUSTER_PATH}/src"

# Transfer files
echo ""
echo "Transferring files..."
echo "  - Definition file..."
scp "${SCRIPT_DIR}/rna-map-optimization.def" "${CLUSTER}:${CLUSTER_PATH}/scripts/optimization/"

echo "  - Environment file..."
scp "${PROJECT_ROOT}/environment_optuna.yml" "${CLUSTER}:${CLUSTER_PATH}/"

echo "  - rna_map package source..."
scp -r "${PROJECT_ROOT}/src/rna_map" "${CLUSTER}:${CLUSTER_PATH}/src/"

echo ""
echo "=========================================="
echo "✅ Files transferred successfully!"
echo "=========================================="
echo ""
echo "On the cluster, run:"
echo "  cd ${CLUSTER_PATH}"
echo "  apptainer build rna-map-optimization.sif scripts/optimization/rna-map-optimization.def"
echo ""
echo "Or with fakeroot (if available):"
echo "  apptainer build --fakeroot rna-map-optimization.sif scripts/optimization/rna-map-optimization.def"
echo ""

