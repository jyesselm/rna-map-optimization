# Configuration File

The `optimization_config.yml` file allows you to configure all optimization parameters, cluster settings, and output preferences without modifying code.

## Usage

### Loading the Configuration

The configuration file can be loaded using PyYAML:

```python
import yaml
from pathlib import Path

config_path = Path("config/optimization_config.yml")
with open(config_path) as f:
    config = yaml.safe_load(f)

# Access settings
n_trials = config["optimization"]["n_trials"]
threads = config["optimization"]["threads"]
container_path = config["cluster"]["container_path"]
```

### Configuration Sections

#### `optimization`

Global optimization settings:
- `n_trials`: Number of Optuna trials (default: 20000)
- `threads`: Number of threads per alignment (default: 8)
- `mapq_cutoff`: MAPQ cutoff for high-quality alignments (default: 20)
- `qscore_cutoff`: Quality score cutoff for bit vectors (default: 25)
- `parameters`: Parameter search ranges for Optuna (see below)

#### `parameters`

Defines the search space for each parameter:
- **Integer parameters**: `min`, `max`, `step`
- **Categorical parameters**: `options` (list of possible values)

Examples:
- `seed_length`: Integer range 10-22, step 2
- `score_min`: Categorical with predefined options
- `sensitivity_mode`: Categorical with mode options

#### `cluster`

Cluster job submission settings:
- `time`: Time limit per job (format: HH:MM:SS)
- `memory`: Memory per job (e.g., "16G")
- `cpus`: CPUs per job
- `email`: Email for notifications (optional)
- `job_name_prefix`: Prefix for job names
- `container_path`: Path to Singularity/Apptainer container (.sif file)

#### `output`

Output directory and file handling:
- `base_dir`: Base directory for all results
- `keep_intermediates`: Whether to keep SAM files and other intermediates

#### `case_overrides`

Optional per-case settings to override defaults:

```yaml
case_overrides:
  case_1:
    n_trials: 500
    time: "48:00:00"
    memory: "32G"
```

## Building Containers

### Docker

```bash
cd scripts/container
docker build -f Dockerfile -t rna-map-optimization:latest ..
```

### Singularity/Apptainer

```bash
cd scripts/container
apptainer build --fakeroot rna-map-optimization.sif rna-map-optimization.def
```

Or from the project root:

```bash
apptainer build --fakeroot rna-map-optimization.sif scripts/container/rna-map-optimization.def
```

### Verify Container

```bash
CONTAINER_PATH=/path/to/rna-map-optimization.sif

# Test imports
apptainer exec $CONTAINER_PATH python -c "import rna_map_mini; print('✓ rna-map-mini installed')"
apptainer exec $CONTAINER_PATH python -c "import optuna; print('✓ optuna installed')"
apptainer exec $CONTAINER_PATH bowtie2 --version
```

## Using Container with Optimization

Set the container path in the config file or via environment variable:

```bash
export CONTAINER_PATH=/path/to/rna-map-optimization.sif
```

The optimization scripts will automatically use the container if `CONTAINER_PATH` is set or if `cluster.container_path` is configured in the config file.

