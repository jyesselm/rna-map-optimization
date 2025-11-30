#!/bin/bash
# Convert Docker image to Apptainer format locally
#
# This script converts a Docker image to Apptainer .sif format
# so you can transfer it to a cluster that doesn't have Docker.
#
# Usage:
#   ./scripts/optimization/convert_docker_to_apptainer.sh [docker_image] [output_path]
#
# Example:
#   ./scripts/optimization/convert_docker_to_apptainer.sh rna-map-optimization:latest rna-map-optimization.sif

set -e

DOCKER_IMAGE="${1:-rna-map-optimization:latest}"
OUTPUT_PATH="${2:-rna-map-optimization.sif}"

echo "=========================================="
echo "Converting Docker Image to Apptainer"
echo "=========================================="
echo ""
echo "Docker image: $DOCKER_IMAGE"
echo "Output: $OUTPUT_PATH"
echo ""

# Check if Docker image exists
if ! docker images | grep -q "$(echo $DOCKER_IMAGE | cut -d: -f1)"; then
    echo "ERROR: Docker image not found: $DOCKER_IMAGE"
    echo "Available images:"
    docker images | head -10
    exit 1
fi

# Check if Apptainer/Singularity is available
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "ERROR: Neither Apptainer nor Singularity found"
    echo ""
    echo "Options:"
    echo "1. Install Apptainer locally:"
    echo "   brew install apptainer"
    echo "   (May require sudo for fakeroot capability)"
    echo ""
    echo "2. Save Docker image and convert on cluster:"
    echo "   docker save $DOCKER_IMAGE | gzip > rna-map-optimization.tar.gz"
    echo "   scp rna-map-optimization.tar.gz user@cluster:/path/"
    echo "   # On cluster (if Docker available):"
    echo "   docker load < rna-map-optimization.tar.gz"
    echo "   apptainer build rna-map-optimization.sif docker-daemon://$DOCKER_IMAGE"
    echo ""
    echo "3. Use the Apptainer definition file directly on cluster:"
    echo "   scp scripts/optimization/rna-map-optimization.def user@cluster:/path/"
    echo "   # On cluster:"
    echo "   apptainer build rna-map-optimization.sif rna-map-optimization.def"
    exit 1
fi

echo "Using: $CONTAINER_CMD"
echo "Converting Docker image to Apptainer format..."
echo ""

# Try fakeroot first, then sudo
if [ "$CONTAINER_CMD" = "apptainer" ] && $CONTAINER_CMD capability list 2>/dev/null | grep -q "fakeroot"; then
    echo "Using fakeroot for conversion..."
    if $CONTAINER_CMD build --fakeroot "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
        echo ""
        echo "✅ Apptainer image created: $OUTPUT_PATH"
    else
        echo "Fakeroot conversion failed, trying sudo..."
        if sudo $CONTAINER_CMD build "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
            echo ""
            echo "✅ Apptainer image created: $OUTPUT_PATH"
        else
            echo "❌ Conversion failed"
            exit 1
        fi
    fi
else
    echo "Using sudo for conversion (you may be prompted for password)..."
    if sudo $CONTAINER_CMD build "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
        echo ""
        echo "✅ Apptainer image created: $OUTPUT_PATH"
    else
        echo "❌ Conversion failed"
        exit 1
    fi
fi

echo ""
echo "=========================================="
echo "✅ Conversion complete!"
echo "=========================================="
echo ""
echo "Container: $OUTPUT_PATH"
echo "Size: $(ls -lh $OUTPUT_PATH | awk '{print $5}')"
echo ""
echo "Transfer to cluster:"
echo "  scp $OUTPUT_PATH user@cluster:/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/"
echo ""

