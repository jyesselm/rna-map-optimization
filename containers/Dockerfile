FROM continuumio/miniconda3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
        build-essential \
        wget \
        curl \
        git \
        && rm -rf /var/lib/apt/lists/*

# Configure git to allow cloning into /tmp (fixes "dubious ownership" error)
# This is needed for pip to clone rna-map-mini from GitHub
RUN git config --global --add safe.directory '*'

# Copy environment file and optimization scripts
COPY environment.yml /opt/environment.yml
COPY scripts/ /opt/scripts/

# Create conda environment from environment.yml
RUN /opt/conda/bin/conda env create -f /opt/environment.yml -n rna-map-optimization

# Activate environment (rna-map-mini will be installed via pip from environment.yml)
# The environment.yml already includes rna-map-mini in pip dependencies
RUN . /opt/conda/etc/profile.d/conda.sh && \
        conda activate rna-map-optimization && \
        echo "rna-map-mini will be installed from GitHub via environment.yml"

# Clean up conda cache
RUN /opt/conda/bin/conda clean -afy

# Create directories for data and results
RUN mkdir -p /data /results /work

# Make optimization scripts executable
RUN chmod +x /opt/scripts/*.sh /opt/scripts/*.py 2>/dev/null || true

# Set PATH to include conda environment bin directory
ENV PATH=/opt/conda/envs/rna-map-optimization/bin:/opt/conda/bin:$PATH
ENV CONDA_DEFAULT_ENV=rna-map-optimization

# Set working directory
WORKDIR /work

# Default command
CMD ["/opt/conda/envs/rna-map-optimization/bin/python"]
