#!/bin/bash
# Build optimization container using Docker locally, then convert to Apptainer
#
# This script builds the container using Docker (which works without sudo on most systems)
# and then converts it to Apptainer format for use on the cluster.
#
# Usage:
#   ./scripts/optimization/build_with_docker.sh [output_path]
#
# Example:
#   ./scripts/optimization/build_with_docker.sh rna-map-optimization.sif
#   ./scripts/optimization/build_with_docker.sh /path/to/rna-map-optimization.sif

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Default output path
OUTPUT_PATH="${1:-${PROJECT_ROOT}/rna-map-optimization.sif}"

echo "=========================================="
echo "Building Optimization Container with Docker"
echo "=========================================="
echo ""
echo "Output: $OUTPUT_PATH"
echo "Project root: $PROJECT_ROOT"
echo ""

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker not found in PATH"
    echo "Please install Docker to build containers locally"
    exit 1
fi

# Check if Apptainer/Singularity is available for conversion
CONTAINER_CMD=""
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "WARNING: Apptainer/Singularity not found"
    echo "Will build Docker image only. You'll need to convert it manually on the cluster."
    CONTAINER_CMD=""
fi

# Check if Dockerfile exists
DOCKERFILE="${SCRIPT_DIR}/Dockerfile"
if [ ! -f "$DOCKERFILE" ]; then
    echo "ERROR: Dockerfile not found: $DOCKERFILE"
    exit 1
fi

# Check if environment file exists
ENV_FILE="${PROJECT_ROOT}/environment.yml"
if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: Environment file not found: $ENV_FILE"
    exit 1
fi

# Verify environment.yml exists (rna-map-mini is installed via pip from GitHub)
ENV_FILE="${PROJECT_ROOT}/environment.yml"
if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: environment.yml not found: $ENV_FILE"
    exit 1
fi

echo "Step 1: Building Docker image..."
echo "This may take 10-20 minutes depending on network speed..."
echo ""

# Build Docker image
cd "$PROJECT_ROOT"
DOCKER_IMAGE="rna-map-optimization:latest"

if docker build -f "$DOCKERFILE" -t "$DOCKER_IMAGE" . 2>&1; then
    echo ""
    echo "✅ Docker image built successfully: $DOCKER_IMAGE"
else
    echo ""
    echo "❌ Docker build failed"
    exit 1
fi

# Convert to Apptainer/Singularity format if available
if [ -n "$CONTAINER_CMD" ]; then
    echo ""
    echo "Step 2: Converting Docker image to Apptainer format..."
    echo "Output: $OUTPUT_PATH"
    echo ""
    
    # Try fakeroot first, then sudo
    if [ "$CONTAINER_CMD" = "apptainer" ] && $CONTAINER_CMD capability list 2>/dev/null | grep -q "fakeroot"; then
        echo "Using fakeroot for conversion..."
        if $CONTAINER_CMD build --fakeroot "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
            echo ""
            echo "✅ Apptainer image created successfully: $OUTPUT_PATH"
        else
            echo "Fakeroot conversion failed, trying sudo..."
            if sudo $CONTAINER_CMD build "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
                echo ""
                echo "✅ Apptainer image created successfully: $OUTPUT_PATH"
            else
                echo ""
                echo "❌ Conversion failed"
                echo "Docker image is available as: $DOCKER_IMAGE"
                echo "You can convert it manually on the cluster with:"
                echo "  apptainer build $OUTPUT_PATH docker-daemon://${DOCKER_IMAGE}"
                exit 1
            fi
        fi
    else
        echo "Using sudo for conversion (you may be prompted for password)..."
        if sudo $CONTAINER_CMD build "$OUTPUT_PATH" "docker-daemon://${DOCKER_IMAGE}" 2>&1; then
            echo ""
            echo "✅ Apptainer image created successfully: $OUTPUT_PATH"
        else
            echo ""
            echo "⚠️  Conversion failed (may need sudo or fakeroot)"
            echo "Docker image is available as: $DOCKER_IMAGE"
            echo "You can convert it manually on the cluster with:"
            echo "  apptainer build $OUTPUT_PATH docker-daemon://${DOCKER_IMAGE}"
            echo ""
            echo "Or transfer the Docker image and convert on cluster:"
            echo "  docker save $DOCKER_IMAGE | gzip > rna-map-optimization.tar.gz"
            echo "  scp rna-map-optimization.tar.gz user@cluster:/path/"
            echo "  # On cluster: docker load < rna-map-optimization.tar.gz"
            echo "  # Then: apptainer build rna-map-optimization.sif docker-daemon://${DOCKER_IMAGE}"
            exit 1
        fi
    fi
else
    echo ""
    echo "⚠️  Apptainer/Singularity not available for conversion"
    echo "Docker image built: $DOCKER_IMAGE"
    echo ""
    echo "Since you don't have Docker on the cluster, build the container on the cluster instead:"
    echo ""
    echo "Option 1: Build from definition file on cluster (RECOMMENDED):"
    echo "  # Transfer files to cluster:"
    echo "  scp scripts/optimization/rna-map-optimization.def user@cluster:/path/"
                echo "  scp environment.yml user@cluster:/path/"
                echo "  # No need to transfer src/rna_map - rna-map-mini installed from GitHub"
    echo "  # On cluster:"
    echo "  apptainer build rna-map-optimization.sif rna-map-optimization.def"
    echo ""
    echo "Option 2: Use conda environment instead (no container needed):"
    echo "  # On cluster:"
    echo "  bash scripts/optimization/setup_optimization_env.sh"
    echo ""
    echo "See scripts/optimization/BUILD_ON_CLUSTER.md for detailed instructions"
    exit 0
fi

echo ""
echo "=========================================="
echo "✅ Build complete!"
echo "=========================================="
echo ""
echo "Container location: $OUTPUT_PATH"
echo ""
echo "To transfer to cluster:"
echo "  scp $OUTPUT_PATH user@cluster:/work/yesselmanlab/jyesselm/installs/rna_map_nextflow/"
echo ""
echo "To test the container:"
echo "  $CONTAINER_CMD exec $OUTPUT_PATH python -c 'import rna_map; print(\"rna_map installed successfully\")'"
echo ""
echo "To use with optimization scripts:"
echo "  export CONTAINER_PATH=$OUTPUT_PATH"
echo "  bash scripts/optimization/submit_optimization_jobs.sh"
echo ""

