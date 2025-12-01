# Container Build Guide

This guide explains how to build Docker and Apptainer/Singularity containers for RNA-MAP optimization.

## Prerequisites

### For Docker:
- Docker installed and running
- Access to pull base images from Docker Hub

### For Apptainer/Singularity:
- Apptainer or Singularity installed
- For `--fakeroot` builds: configured fakeroot support
- For root builds: root/sudo access

## Quick Start

### Docker

```bash
# From the project root directory
docker build -f scripts/container/Dockerfile -t rna-map-optimization:latest .
```

### Apptainer/Singularity

```bash
# From the project root directory (with fakeroot)
apptainer build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def

# Or without fakeroot (requires root/sudo)
sudo singularity build rna-map-optimization.sif scripts/container/rna-map-optimization.def
```

## Detailed Build Instructions

### Option 1: Build Docker Image

#### Basic Build

```bash
cd /path/to/rna-map-optimization
docker build -f scripts/container/Dockerfile -t rna-map-optimization:latest .
```

#### Build with Custom Tag

```bash
docker build -f scripts/container/Dockerfile \
    -t rna-map-optimization:v1.0.0 \
    -t rna-map-optimization:latest \
    .
```

#### Build Arguments (if needed in future)

The Dockerfile currently uses default paths. If you need to customize:

```bash
docker build -f scripts/container/Dockerfile \
    --build-arg ENV_FILE=environment.yml \
    -t rna-map-optimization:latest \
    .
```

#### Verify Docker Build

```bash
# Test that the container runs
docker run --rm rna-map-optimization:latest python --version

# Test imports
docker run --rm rna-map-optimization:latest python -c "import rna_map_mini; print('✓ rna-map-mini')"
docker run --rm rna-map-optimization:latest python -c "import optuna; print('✓ optuna')"

# Test bowtie2
docker run --rm rna-map-optimization:latest bowtie2 --version
```

### Option 2: Build Apptainer/Singularity Image

#### Build with Fakeroot (Recommended, if available)

Fakeroot allows building without root privileges:

```bash
cd /path/to/rna-map-optimization

# Build with fakeroot
apptainer build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def

# Or with Singularity
singularity build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def
```

**Requirements for fakeroot:**
- Apptainer 3.0+ or Singularity 3.5+
- `/etc/subuid` and `/etc/subgid` configured (see below)

**Setup fakeroot (one-time, on Linux):**
```bash
# Check if configured
cat /etc/subuid
cat /etc/subgid

# If not configured, edit these files to add your user:
# /etc/subuid: your_username:100000:65536
# /etc/subgid: your_username:100000:65536
```

#### Build with Root/Sudo

If fakeroot is not available:

```bash
cd /path/to/rna-map-optimization

# Build with sudo
sudo apptainer build rna-map-optimization.sif scripts/container/rna-map-optimization.def

# Or with Singularity
sudo singularity build rna-map-optimization.sif scripts/container/rna-map-optimization.def
```

#### Build from Docker Image (Alternative)

If you've already built a Docker image:

```bash
# First build Docker image (see Option 1)
docker build -f scripts/container/Dockerfile -t rna-map-optimization:latest .

# Then convert to Apptainer/Singularity
apptainer build rna-map-optimization.sif docker-daemon://rna-map-optimization:latest

# Or with Singularity
singularity build rna-map-optimization.sif docker-daemon://rna-map-optimization:latest
```

#### Verify Apptainer/Singularity Build

```bash
# Test that the container runs
apptainer exec rna-map-optimization.sif python --version

# Test imports
apptainer exec rna-map-optimization.sif python -c "import rna_map_mini; print('✓ rna-map-mini')"
apptainer exec rna-map-optimization.sif python -c "import optuna; print('✓ optuna')"

# Test bowtie2
apptainer exec rna-map-optimization.sif bowtie2 --version
```

## What Gets Installed

Both containers include:

1. **Base System**:
   - Ubuntu-based system with Miniconda3
   - Build tools (gcc, make, etc.)
   - Git, wget, curl

2. **Conda Environment** (`rna-map-optimization`):
   - Python >= 3.10
   - Bowtie2 (from bioconda)
   - pandas, numpy, pyyaml, tabulate

3. **Python Packages** (via pip):
   - rna-map-mini (from GitHub)
   - optuna >= 3.0
   - plotly >= 5.0

4. **Optimization Scripts**:
   - All scripts from `scripts/` directory
   - Made executable automatically

## File Structure in Container

```
/opt/
├── environment.yml          # Conda environment file
├── scripts/                 # Optimization scripts
└── conda/
    └── envs/
        └── rna-map-optimization/  # Conda environment

/data                        # Mount point for input data
/results                     # Mount point for output results
/work                        # Working directory
```

## Building for Different Architectures

### Docker Multi-Platform Build

```bash
# Build for multiple platforms (requires buildx)
docker buildx build --platform linux/amd64,linux/arm64 \
    -f scripts/container/Dockerfile \
    -t rna-map-optimization:latest \
    .
```

### Apptainer/Singularity

Apptainer/Singularity builds for the host architecture automatically. For different architectures, build on a machine with that architecture.

## Troubleshooting

### Docker Build Issues

**Issue: "Cannot connect to Docker daemon"**
```bash
# Start Docker service
sudo systemctl start docker

# Or on macOS/Windows, start Docker Desktop
```

**Issue: "Out of space during build"**
```bash
# Clean up Docker
docker system prune -a

# Check disk space
df -h
```

**Issue: "Failed to fetch" (network issues)**
```bash
# Use build with --network=host
docker build --network=host -f scripts/container/Dockerfile -t rna-map-optimization:latest .
```

### Apptainer/Singularity Build Issues

**Issue: "fakeroot: command not found"**
- Install fakeroot: `sudo apt-get install fakeroot` (Debian/Ubuntu)
- Or use sudo build instead

**Issue: "Failed to get section"**
- Check that the .def file is valid YAML/Singularity format
- Ensure all file paths in %files section exist

**Issue: "Permission denied"**
- Use `sudo` for root builds
- Or configure fakeroot properly

**Issue: "Failed to pull Docker image"**
- Check internet connection
- Verify Docker Hub is accessible
- Try building from a local Docker image instead

## Container Size Optimization

### Docker

The current Dockerfile already optimizes size by:
- Using miniconda instead of full conda
- Cleaning conda cache (`conda clean -afy`)
- Removing apt cache (`rm -rf /var/lib/apt/lists/*`)

For further optimization, use multi-stage builds (advanced).

### Apptainer/Singularity

Apptainer/Singularity images are already compressed. To reduce size further:
- Remove unnecessary files before building
- Use `--sandbox` for development, convert to `.sif` for production

## Transferring Containers

### Docker

```bash
# Save to tar file
docker save rna-map-optimization:latest | gzip > rna-map-optimization.tar.gz

# Transfer to another machine
scp rna-map-optimization.tar.gz user@remote:/path/

# Load on remote machine
gunzip -c rna-map-optimization.tar.gz | docker load
```

### Apptainer/Singularity

```bash
# Transfer .sif file (already compressed)
scp rna-map-optimization.sif user@remote:/path/

# Use directly - no need to load
apptainer exec /path/to/rna-map-optimization.sif python script.py
```

## Using Containers

See [CONTAINER_USAGE.md](./CONTAINER_USAGE.md) for detailed usage instructions.

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Build Container

on:
  push:
    tags:
      - 'v*'

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Build Docker image
        run: |
          docker build -f scripts/container/Dockerfile \
            -t rna-map-optimization:${{ github.ref_name }} \
            -t rna-map-optimization:latest \
            .
      - name: Push to registry
        run: |
          docker push rna-map-optimization:${{ github.ref_name }}
          docker push rna-map-optimization:latest
```

## Version Tags

Recommended versioning:
- `rna-map-optimization:latest` - Latest build
- `rna-map-optimization:v1.0.0` - Specific version
- `rna-map-optimization:dev` - Development builds

## Support

For issues:
1. Check that all prerequisites are installed
2. Verify file paths are correct
3. Check disk space and permissions
4. Review build logs for specific errors

