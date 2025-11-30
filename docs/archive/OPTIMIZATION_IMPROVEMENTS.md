# Optimization Script Improvements

This document summarizes the improvements made to the Optuna-based Bowtie2 parameter optimization script based on test results.

## Key Improvements

### 1. **Added Gap Penalties** ✅
   - **Read gap penalty** (`--rdg`): Options include `None`, `5,3`, `6,4`, `8,4`
   - **Reference gap penalty** (`--rfg`): Options include `None`, `5,3`, `6,4`, `8,4`
   - Gap penalties are important for RNA alignments where gaps are common

### 2. **Expanded Parameter Options** ✅
   - **Score minimum**: Now includes **all three Bowtie2 score function types**:
     - **L-type (Linear)**: `L,0,0.2`, `L,0,0.3`, `L,5,0.1`, `L,10,0.2`, `L,15,0.1`
     - **L-type (negative, for compatible cases)**: `L,0,-0.6`, `L,0,-0.4`, `L,0,-0.2`, `L,-0.6,-0.6`
     - **G-type (Max)**: `G,15,0.1`, `G,20,8`, `G,20,15`, `G,25,10`, `G,30,12`, `G,30,15`
     - **S-type (Sum)**: `S,10,0.1`, `S,15,0.2`, `S,20,0.15`, `S,25,0.2`, `S,30,0.1`
     - Total: **20+ score function options** across all three types
   - **Compatibility validation**: Automatically prunes incompatible combinations
     - L-type with negative values incompatible with mismatch penalty > 0
   - **Mismatch penalty**: Expanded to include:
     - `None`, `6,2`, `4,2`, `2,2`
   - **Sensitivity modes**: Added `very-sensitive-local` option

### 3. **Optimizable Threads** ✅
   - Added `--optimize-threads` flag to optionally optimize thread count
   - When enabled, threads are optimized in range [2, min(threads*2, 16)]
   - Note: More threads = faster alignment, but doesn't affect quality metrics
   - Default: Use fixed thread count (faster optimization)

### 4. **Constrained Non-Varying Parameters** ✅
   - **Seed mismatches**: Constrained to 0-1 (was 0-2, 0 was always used)
   - Only added to params if > 0 (reduces parameter space)
   - **Seed length**: Expanded range to 10-22 (from 10-20)
   - **Max insert size**: Expanded range to 200-1200 (from 200-1000)

### 5. **Improved Early Pruning** ✅
   - Added multiple early pruning checks to avoid wasting time on bad trials:
     - No reads detected
     - Low alignment rate (< 30%, was 50%)
     - No valid alignments in SAM file
     - No accepted bit vectors
   - Added `failure_reason` attribute to failed trials for debugging

### 6. **Better Parameter Tracking** ✅
   - All optional parameters are properly tracked in results
   - CSV summary includes all parameters (including None values)
   - Gap penalties and other new parameters are included

## Usage Examples

### Basic usage (fixed threads):
```bash
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 100 \
    --threads 4 \
    --output-dir optimization_optuna
```

### With thread optimization:
```bash
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 100 \
    --threads 8 \
    --optimize-threads \
    --output-dir optimization_optuna
```

## Parameter Space

The optimization now explores:

1. **seed_length**: 10-22 (step 2)
2. **seed_mismatches**: 0-1 (only if > 0)
3. **maxins**: 200-1200 (step 100)
4. **score_min**: **20+ categorical options** across all three types:
   - L-type (Linear): 5 positive + 4 negative options
   - G-type (Max): 6 options
   - S-type (Sum): 5 options
5. **mismatch_penalty**: 4 categorical options
6. **gap_penalty_read**: 4 categorical options
7. **gap_penalty_ref**: 4 categorical options
8. **sensitivity_mode**: 5 categorical options
9. **threads**: 2-16 (if optimizing threads)

**Total parameter combinations**: Significantly expanded with comprehensive score function coverage

**Total parameter space**: Significantly expanded with more meaningful options for RNA alignment optimization.

## Improvements Based on Test Results

From the initial test run:
- **Trial 0 succeeded**: Used `fast-local`, `mismatch_penalty=4,2`, `seed_length=12`
- **Trial 1 failed**: Used `score_min=L,0,-0.4` which may have been too restrictive

The improvements address:
1. **More options**: Expanded categorical choices allow Optuna to find better combinations
2. **Gap penalties**: Important for RNA but were missing
3. **Better pruning**: Failed trials are identified earlier
4. **Constrained space**: Removed parameters that don't vary (seed_mismatches=0)

## Expected Impact

- **Better results**: More parameter options should lead to better signal-to-noise ratios
- **Faster optimization**: Early pruning reduces time spent on failed trials
- **More informative**: Failure reasons help understand what doesn't work
- **Flexible threading**: Option to optimize threads if desired (though quality doesn't depend on threads)

