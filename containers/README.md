# Container Documentation

This directory contains files and documentation for building Docker and Apptainer/Singularity containers for RNA-MAP optimization.

## Files

- **`Dockerfile`** - Docker container definition
- **`rna-map-optimization.def`** - Apptainer/Singularity container definition
- **`BUILD.md`** - Complete guide for building containers (Docker and Apptainer)
- **`CONTAINER_USAGE.md`** - Guide for using containers on clusters
- **`build_with_docker.sh`** - Script to build Docker image and convert to Apptainer
- **`build_optimization_container.sh`** - Script to build Apptainer container directly

## Quick Start

### Build with Docker

```bash
# From project root
docker build -f scripts/container/Dockerfile -t rna-map-optimization:latest .
```

### Build with Apptainer

```bash
# From project root (with fakeroot)
apptainer build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def

# Or with sudo
sudo apptainer build rna-map-optimization.sif scripts/container/rna-map-optimization.def
```

### Verify Build

```bash
apptainer exec rna-map-optimization.sif python -c "import rna_map_mini; print('✓ rna-map-mini')"
apptainer exec rna-map-optimization.sif python -c "import optuna; print('✓ optuna')"
apptainer exec rna-map-optimization.sif bowtie2 --version
```

## Documentation

### For Building Containers
- **[BUILD.md](./BUILD.md)** - Complete build instructions for Docker and Apptainer
  - Prerequisites
  - Step-by-step build instructions
  - Troubleshooting
  - Multi-platform builds
  - Container size optimization

### For Using Containers
- **[CONTAINER_USAGE.md](./CONTAINER_USAGE.md)** - How to use containers for optimization
  - Setting up containers on clusters
  - Running optimization jobs
  - Mount points and paths
  - Job script examples
  - Best practices

## Container Contents

Both containers include:
- **Base**: Ubuntu-based system with Miniconda3
- **Python**: Python 3.10+ with all dependencies
- **Bowtie2**: Latest version from bioconda
- **rna-map-mini**: Installed from GitHub via pip
- **Optuna**: >= 3.0 for Bayesian optimization
- **Plotly**: For visualization
- **All scripts**: From `scripts/` directory

## What's Included

### System Packages
- Build tools (gcc, make, etc.)
- Git, wget, curl
- All required libraries

### Python Packages
- rna-map-mini (from GitHub)
- optuna >= 3.0
- plotly >= 5.0
- pandas >= 1.5
- numpy >= 1.21
- pyyaml >= 6.0
- tabulate >= 0.9

### Bioinformatic Tools
- Bowtie2 (latest from bioconda)

## Build Scripts

### `build_with_docker.sh`
Builds Docker image locally, then converts to Apptainer format.

**Usage:**
```bash
./scripts/container/build_with_docker.sh [output_path]
```

**Requirements:**
- Docker installed
- Apptainer/Singularity (optional, for conversion)

### `build_optimization_container.sh`
Builds Apptainer/Singularity container directly.

**Usage:**
```bash
./scripts/container/build_optimization_container.sh [output_path]
```

**Requirements:**
- Apptainer/Singularity installed
- Fakeroot capability or sudo access

## Container Paths

Inside the container:
- `/opt/environment.yml` - Conda environment file
- `/opt/scripts/` - All optimization scripts
- `/data` - Mount point for input data
- `/results` - Mount point for output results
- `/work` - Working directory

## Version History

- **v1.0.0**: Initial release
  - Uses `rna-map-mini` from GitHub
  - Includes all optimization scripts
  - Python 3.10+, Optuna 3.0+, Bowtie2

## Support

For issues:
1. Check [BUILD.md](./BUILD.md) for build problems
2. Check [CONTAINER_USAGE.md](./CONTAINER_USAGE.md) for usage issues
3. Verify all prerequisites are installed
4. Review build logs for specific errors

## Related Documentation

- [Main README](../../README.md) - Project overview
- [Installation Guide](../../INSTALL.md) - Local installation
- [Configuration Guide](../../config/README.md) - Optimization configuration

