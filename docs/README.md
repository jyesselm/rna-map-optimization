# Bowtie2 Parameter Optimization

Documentation for optimizing Bowtie2 alignment parameters for RNA-MAP workflows.

## Quick Start

### Use Recommended Parameters

Recommended parameters from optimization are available in:
- **[best_parameters.txt](./recommended_params/best_parameters.txt)** - Full parameter string
- **[BEST_PARAMETERS.md](./BEST_PARAMETERS.md)** - Detailed breakdown

### Run Your Own Optimization

1. **Install Optuna** (if not already installed):
   ```bash
   pip install optuna plotly
   # conda env create -f environment.yml
   ```

2. **Run optimization**:
   ```bash
   python3 scripts/optimize_bowtie2_params_optuna.py \
       --fasta reference.fasta \
       --fastq1 reads_R1.fastq \
       --fastq2 reads_R2.fastq \
       --n-trials 100 \
       --output-dir optimization_results
   ```

3. **Analyze results**: Use `scripts/analyze_top_parameters.py`

See **[BOWTIE2_OPTIMIZATION.md](../BOWTIE2_OPTIMIZATION.md)** for complete guide.

## Key Results

### Optimal Parameters Summary

From analysis of top 100 parameter combinations:

**Constants (use these always):**
- Seed length: 18
- Mismatch penalty: 6,2
- Gap penalties: read=8,4, ref=5,3
- Sensitivity: fast-local
- Min insert size: 50

**Recommended defaults:**
- Max insert size: 200
- Seed mismatches: 1
- Score minimum: L,10,0.2
- Repetitive effort: 4

See **[TOP_100_PARAMETER_ANALYSIS.md](./TOP_100_PARAMETER_ANALYSIS.md)** for complete analysis.

## Files

- **[BOWTIE2_OPTIMIZATION.md](../BOWTIE2_OPTIMIZATION.md)** - Main optimization guide
- **[BEST_PARAMETERS.md](./BEST_PARAMETERS.md)** - Recommended parameters
- **[TOP_100_PARAMETER_ANALYSIS.md](./TOP_100_PARAMETER_ANALYSIS.md)** - Analysis results
- **[recommended_params/](./recommended_params/)** - Parameter files
- **[examples/](./examples/)** - Example scripts
- **[archive/](./archive/)** - Historical detailed documentation

## Scripts

- `scripts/optimize_bowtie2_params.py` - Grid search
- `scripts/optimize_bowtie2_params_optuna.py` - Bayesian optimization (recommended)
- `scripts/analyze_top_parameters.py` - Analysis tool
