#!/bin/bash
# Submit optimization jobs for all test cases in data/ directory

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
DATA_DIR="${PROJECT_ROOT}/data"
CONFIG_FILE="${SCRIPT_DIR}/cluster_optimization_config.yml"
ENV_NAME="rna-map-optimization"

# Default values (can be overridden by config)
MAX_COMBINATIONS=200
READ_LENGTH=150
THREADS=8
MAPQ_CUTOFF=20
TIME="24:00:00"
MEMORY="16G"
CPUS=8
EMAIL=""
JOB_NAME_PREFIX="bt2_opt"
OUTPUT_BASE_DIR="optimization_results"

echo "=========================================="
echo "Bowtie2 Optimization Job Submission"
echo "=========================================="
echo ""
echo "Project root: ${PROJECT_ROOT}"
echo "Data directory: ${DATA_DIR}"
echo "Config file: ${CONFIG_FILE}"
echo ""

# Check if data directory exists
if [ ! -d "${DATA_DIR}" ]; then
    echo "ERROR: Data directory not found: ${DATA_DIR}"
    echo "Please create the data/ directory and add test cases"
    exit 1
fi

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "WARNING: Config file not found: ${CONFIG_FILE}"
    echo "Using default settings"
else
    echo "Loading configuration from ${CONFIG_FILE}..."
    # Parse YAML config (simple parsing - assumes Python is available)
    if command -v python3 &> /dev/null; then
        eval "$(python3 << EOF
import yaml
import sys
with open('${CONFIG_FILE}', 'r') as f:
    config = yaml.safe_load(f)
    
opt = config.get('optimization', {})
cluster = config.get('cluster', {})
output = config.get('output', {})

print(f"MAX_COMBINATIONS={opt.get('max_combinations', 200)}")
print(f"READ_LENGTH={opt.get('read_length', 150)}")
print(f"THREADS={opt.get('threads', 8)}")
print(f"MAPQ_CUTOFF={opt.get('mapq_cutoff', 20)}")
print(f"TIME='{cluster.get('time', '24:00:00')}'")
print(f"MEMORY='{cluster.get('memory', '16G')}'")
print(f"CPUS={cluster.get('cpus', 8)}")
print(f"EMAIL='{cluster.get('email', '') or ''}'")
print(f"JOB_NAME_PREFIX='{cluster.get('job_name_prefix', 'bt2_opt')}'")
print(f"OUTPUT_BASE_DIR='{output.get('base_dir', 'optimization_results')}'")
container_path = cluster.get('container_path')
if container_path:
    print(f"CONTAINER_PATH_CONFIG='{container_path}'")
EOF
)"
    fi
fi

# Check for container or conda environment
# Use CONTAINER_PATH from environment if set, otherwise from config file
CONTAINER_PATH="${CONTAINER_PATH:-${CONTAINER_PATH_CONFIG:-}}"
USE_CONTAINER=false

if [ -n "${CONTAINER_PATH}" ] && [ -f "${CONTAINER_PATH}" ]; then
    USE_CONTAINER=true
    echo "Using container: ${CONTAINER_PATH}"
    
    # Check if container command is available
    if ! command -v apptainer &> /dev/null && ! command -v singularity &> /dev/null; then
        echo "ERROR: Container specified but Apptainer/Singularity not found"
        exit 1
    fi
elif ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "ERROR: Neither container nor conda environment '${ENV_NAME}' found"
    echo "Please either:"
    echo "  1. Build container: bash ${SCRIPT_DIR}/build_optimization_container.sh"
    echo "  2. Or setup conda env: bash ${SCRIPT_DIR}/setup_optimization_env.sh"
    exit 1
else
    echo "Using conda environment: ${ENV_NAME}"
fi

# Find all test cases in data/
echo "Scanning for test cases in ${DATA_DIR}..."
echo ""

TEST_CASES=()
while IFS= read -r -d '' case_dir; do
    case_name=$(basename "${case_dir}")
    
    # Check for required files
    fasta_file=$(find "${case_dir}" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -1)
    fastq1_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R1*.fastq*" -o -name "*_1.fastq*" -o -name "*.fastq*" \) | head -1)
    fastq2_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R2*.fastq*" -o -name "*_2.fastq*" \) | head -1)
    
    if [ -z "${fasta_file}" ] || [ -z "${fastq1_file}" ]; then
        echo "⚠️  Skipping ${case_name}: missing required files (fasta or fastq1)"
        continue
    fi
    
    TEST_CASES+=("${case_name}")
    echo "  ✓ Found case: ${case_name}"
    echo "    FASTA: $(basename "${fasta_file}")"
    echo "    FASTQ1: $(basename "${fastq1_file}")"
    if [ -n "${fastq2_file}" ]; then
        echo "    FASTQ2: $(basename "${fastq2_file}")"
    fi
    echo ""
done < <(find "${DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

if [ ${#TEST_CASES[@]} -eq 0 ]; then
    echo "ERROR: No valid test cases found in ${DATA_DIR}"
    echo "Each case should be a subdirectory with:"
    echo "  - *.fasta or *.fa (reference sequence)"
    echo "  - *_R1*.fastq* or *_1.fastq* (read 1)"
    echo "  - *_R2*.fastq* or *_2.fastq* (read 2, optional)"
    exit 1
fi

echo "Found ${#TEST_CASES[@]} test case(s)"
echo ""

# Create output directory
OUTPUT_DIR="${PROJECT_ROOT}/${OUTPUT_BASE_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Create job submission directory
JOB_DIR="${OUTPUT_DIR}/job_scripts"
mkdir -p "${JOB_DIR}"

# Submit jobs
echo "Submitting optimization jobs..."
echo ""

for case_name in "${TEST_CASES[@]}"; do
    case_dir="${DATA_DIR}/${case_name}"
    fasta_file=$(find "${case_dir}" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -1)
    fastq1_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R1*.fastq*" -o -name "*_1.fastq*" -o -name "*.fastq*" \) | head -1)
    fastq2_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R2*.fastq*" -o -name "*_2.fastq*" \) | head -1)
    
    output_case_dir="${OUTPUT_DIR}/${case_name}"
    job_script="${JOB_DIR}/${case_name}_optimization.sh"
    job_name="${JOB_NAME_PREFIX}_${case_name}"
    
    # Determine container command
    if [ "$USE_CONTAINER" = true ]; then
        if command -v apptainer &> /dev/null; then
            CONTAINER_CMD="apptainer"
        else
            CONTAINER_CMD="singularity"
        fi
    fi
    
    # Build SLURM script
    cat > "${job_script}" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEMORY}
#SBATCH --output=${output_case_dir}/%j.out
#SBATCH --error=${output_case_dir}/%j.err
$( [ -n "${EMAIL}" ] && echo "#SBATCH --mail-user=${EMAIL}" && echo "#SBATCH --mail-type=ALL" )

# Print job info
echo "=========================================="
echo "Job: ${job_name}"
echo "Case: ${case_name}"
echo "Started: \$(date)"
echo "Node: \$(hostname)"
echo "=========================================="
echo ""

# Create output directory
mkdir -p ${output_case_dir}

$(if [ "$USE_CONTAINER" = true ]; then cat << CONTAINER_EOF
# Using container: ${CONTAINER_PATH}
CONTAINER_CMD="${CONTAINER_CMD}"
CONTAINER_PATH="${CONTAINER_PATH}"

# Mount data and results directories
BIND_ARGS="-B ${DATA_DIR}:/data -B ${OUTPUT_DIR}:/results -B ${PROJECT_ROOT}:/work"

# Run optimization inside container
echo "Running optimization for ${case_name} in container..."
echo "FASTA: ${fasta_file}"
echo "FASTQ1: ${fastq1_file}"
$( [ -n "${fastq2_file}" ] && echo "echo \"FASTQ2: ${fastq2_file}\"" )
echo ""

\${CONTAINER_CMD} exec \${BIND_ARGS} \${CONTAINER_PATH} python /work/scripts/optimize_bowtie2_params.py \\
    --fasta /data/${case_name}/$(basename "${fasta_file}") \\
    --fastq1 /data/${case_name}/$(basename "${fastq1_file}") \\
    $( [ -n "${fastq2_file}" ] && echo "--fastq2 /data/${case_name}/$(basename "${fastq2_file}") \\" ) \\
    --output-dir /results/${case_name} \\
    --read-length ${READ_LENGTH} \\
    --threads ${THREADS} \\
    --mapq-cutoff ${MAPQ_CUTOFF} \\
    --max-combinations ${MAX_COMBINATIONS}
CONTAINER_EOF
else cat << CONDA_EOF
# Activate conda environment
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ${ENV_NAME}

# Change to project directory
cd ${PROJECT_ROOT}

# Ensure rna_map is installed from src/rna_map
cd src/rna_map && pip install -e . --quiet && cd ../..

# Run optimization
echo "Running optimization for ${case_name}..."
echo "FASTA: ${fasta_file}"
echo "FASTQ1: ${fastq1_file}"
$( [ -n "${fastq2_file}" ] && echo "echo \"FASTQ2: ${fastq2_file}\"" )
echo ""

python scripts/optimize_bowtie2_params.py \\
    --fasta ${fasta_file} \\
    --fastq1 ${fastq1_file} \\
    $( [ -n "${fastq2_file}" ] && echo "--fastq2 ${fastq2_file} \\" ) \\
    --output-dir ${output_case_dir} \\
    --read-length ${READ_LENGTH} \\
    --threads ${THREADS} \\
    --mapq-cutoff ${MAPQ_CUTOFF} \\
    --max-combinations ${MAX_COMBINATIONS}
CONDA_EOF
fi)

EXIT_CODE=\$?

echo ""
echo "=========================================="
echo "Job completed: \$(date)"
echo "Exit code: \${EXIT_CODE}"
echo "=========================================="

exit \${EXIT_CODE}
EOF

    chmod +x "${job_script}"
    
    # Submit job
    echo "Submitting job for ${case_name}..."
    JOB_ID=$(sbatch "${job_script}" | awk '{print $4}')
    echo "  Job ID: ${JOB_ID}"
    echo "  Script: ${job_script}"
    echo "  Output: ${output_case_dir}"
    echo ""
done

echo "=========================================="
echo "✅ All jobs submitted!"
echo "=========================================="
echo ""
echo "Submitted ${#TEST_CASES[@]} job(s)"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo ""
echo "Check results in:"
echo "  ${OUTPUT_DIR}"
echo ""

