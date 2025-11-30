#!/bin/bash
# Validate data directory structure before submitting jobs

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
DATA_DIR="${PROJECT_ROOT}/data"

echo "=========================================="
echo "Validating Data Directory Structure"
echo "=========================================="
echo ""
echo "Data directory: ${DATA_DIR}"
echo ""

if [ ! -d "${DATA_DIR}" ]; then
    echo "❌ Data directory not found: ${DATA_DIR}"
    echo ""
    echo "Please create the data directory:"
    echo "  mkdir -p ${DATA_DIR}"
    echo ""
    echo "Then add test cases as subdirectories:"
    echo "  ${DATA_DIR}/"
    echo "  ├── case_1/"
    echo "  │   ├── reference.fasta"
    echo "  │   ├── reads_R1.fastq.gz"
    echo "  │   └── reads_R2.fastq.gz"
    echo "  └── case_2/"
    echo "      └── ..."
    exit 1
fi

echo "Scanning for test cases..."
echo ""

VALID_CASES=0
INVALID_CASES=0

while IFS= read -r -d '' case_dir; do
    case_name=$(basename "${case_dir}")
    issues=()
    
    # Check for FASTA file
    fasta_file=$(find "${case_dir}" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -1)
    if [ -z "${fasta_file}" ]; then
        issues+=("Missing FASTA file (*.fasta or *.fa)")
    fi
    
    # Check for FASTQ1 file
    fastq1_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R1*.fastq*" -o -name "*_1.fastq*" -o -name "*.fastq*" \) | head -1)
    if [ -z "${fastq1_file}" ]; then
        issues+=("Missing FASTQ1 file (*_R1*.fastq* or *_1.fastq*)")
    fi
    
    # Check for FASTQ2 file (optional but recommended for paired-end)
    fastq2_file=$(find "${case_dir}" -maxdepth 1 \( -name "*_R2*.fastq*" -o -name "*_2.fastq*" \) | head -1)
    
    if [ ${#issues[@]} -eq 0 ]; then
        VALID_CASES=$((VALID_CASES + 1))
        echo "✅ ${case_name}"
        echo "   FASTA: $(basename "${fasta_file}")"
        echo "   FASTQ1: $(basename "${fastq1_file}")"
        if [ -n "${fastq2_file}" ]; then
            echo "   FASTQ2: $(basename "${fastq2_file}") (paired-end)"
        else
            echo "   FASTQ2: Not found (single-end mode)"
        fi
        echo ""
    else
        INVALID_CASES=$((INVALID_CASES + 1))
        echo "❌ ${case_name}"
        for issue in "${issues[@]}"; do
            echo "   - ${issue}"
        done
        echo ""
    fi
done < <(find "${DATA_DIR}" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Valid cases: ${VALID_CASES}"
echo "Invalid cases: ${INVALID_CASES}"
echo ""

if [ ${VALID_CASES} -eq 0 ]; then
    echo "❌ No valid test cases found!"
    exit 1
elif [ ${INVALID_CASES} -gt 0 ]; then
    echo "⚠️  Some cases have issues. Fix them before submitting jobs."
    exit 1
else
    echo "✅ All cases are valid. Ready to submit jobs!"
    exit 0
fi

