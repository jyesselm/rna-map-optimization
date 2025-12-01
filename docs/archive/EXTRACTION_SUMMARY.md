# Optimization Repository Extraction Summary

## What Was Done

Successfully extracted the optimization code from the main `rna-map-nextflow` repository into a standalone `rna-map-optimization` repository.

## Repository Structure

```
rna-map-optimization/
├── README.md                      # Updated for standalone repo
├── INSTALL.md                     # Installation guide
├── LICENSE                        # Apache 2.0 (from main repo)
├── pyproject.toml                 # Python package config
├── environment.yml                # Conda environment
├── .gitignore                     # Git ignore rules
│
├── scripts/                       # Optimization scripts
│   ├── optimize_bowtie2_params_optuna.py  # Optuna optimization
│   ├── optimize_bowtie2_params.py         # Grid search
│   ├── run_baseline_test.py               # Baseline comparison
│   ├── collect_top_results.py            # Result aggregation
│   ├── cluster/                           # Cluster scripts
│   └── container/                         # Container scripts
│
├── docs/                          # Documentation
│   ├── BEST_PARAMETERS.md
│   ├── TOP_100_PARAMETER_ANALYSIS.md
│   ├── recommended_params/
│   └── archive/
│
├── config/                        # Configuration
│   └── slurm_optimized.config
│
├── src/                           # Package structure
│   └── rna_map_optimization/
│       └── __init__.py
│
└── test/                          # Test cases and results
    └── test_cases/
```

## Key Changes

### 1. Updated Imports
- Removed `PROJECT_ROOT` references
- Updated import comments to reference rna-map-mini package installation
- All scripts now import from installed `rna-map-mini` package

### 2. Package Structure
- Created `src/rna_map_optimization/` package structure
- Added `pyproject.toml` with dependency on `rna-map-mini>=0.1.0`
- Set up for potential PyPI publishing

### 3. Documentation
- Updated README for standalone repository
- Added INSTALL.md with detailed installation instructions
- Updated environment.yml with installation notes

### 4. Organization
- Moved container scripts to `scripts/container/`
- Moved cluster scripts to `scripts/cluster/`
- Organized test cases under `test/test_cases/`

### 5. CI/CD
- Added GitHub Actions workflow for basic testing
- Tests import verification and script syntax

## Dependencies

The optimization repository depends on:
- **rna-map-mini** package (required)
- **optuna** for Bayesian optimization
- **plotly** for visualization
- **pandas, numpy, pyyaml** for data handling
- **bowtie2** for alignment

## Installation

Users must install the `rna-map-mini` package first:

```bash
pip install git+https://github.com/jyesselm/rna-map-mini.git
```

Then install optimization toolkit:

```bash
conda env create -f environment.yml
conda activate rna-map-optimization
```

## Next Steps

### For Main Repository

1. **Update README.md** to reference optimization repository
2. **Remove `optimization/` directory** (after confirming extraction is complete)
3. **Update documentation** to point to optimization repo
4. **Add note** about optimization being in separate repo

### For Optimization Repository

1. **Create GitHub repository** (if not already created)
2. **Push initial commit** to remote
3. **Set up GitHub Pages** (optional, for documentation)
4. **Add badges** to README (build status, etc.)
5. **Create releases** for versioning

## Interaction Model

```
┌─────────────────────────────┐
│  rna-map-optimization       │
│  (This repository)          │
│                             │
│  - Optimization scripts     │
│  - Parameter analysis       │
│  - Recommended params       │
└──────────┬──────────────────┘
           │
           │ Depends on
           ▼
┌─────────────────────────────┐
│  rna-map-nextflow           │
│  (Main repository)           │
│                             │
│  - Python package (rna-map) │
│  - Nextflow workflow        │
│  - Pipeline execution        │
└─────────────────────────────┘
```

## Benefits

1. **Clean Separation**: Main repo focuses on production workflow
2. **Independent Versioning**: Optimization can evolve separately
3. **Reduced Complexity**: Main repo stays focused
4. **Better Collaboration**: Different teams can work on different repos
5. **Easier Distribution**: Users who don't need optimization don't need to install it

## Files Status

- ✅ All optimization scripts copied and updated
- ✅ Documentation updated for standalone use
- ✅ Package structure created
- ✅ Dependencies configured
- ✅ Git repository initialized and committed
- ⏳ Ready to push to remote (when GitHub repo is created)

## Notes

- Test results in `test/` are included but large result files are gitignored
- Container and cluster scripts are organized in subdirectories
- All scripts have been updated to work with installed `rna-map` package
- Installation instructions are clear and comprehensive

