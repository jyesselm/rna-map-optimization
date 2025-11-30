# Added Parameters to Optimization Script

All high-priority and medium-priority parameters have been successfully added to the optimization script!

## ✅ Newly Added Parameters

### High Priority (RNA-Specific)

1. **Seed Interval** (`-i`)
   - Controls frequency of seed extraction from reads
   - Options: `None`, `S,1,0.5`, `S,1,0.75`, `S,1,1.15`, `S,1,1.5`, `S,1,2.0`
   - Lower values = more seeds = more sensitive but slower

2. **Non-A/C/G/T Penalty** (`--np`)
   - Penalty for non-standard nucleotides (N, I, etc.)
   - Options: `None`, `0`, `1`, `2`, `3`
   - Critical for RNA data with modified bases

3. **Non-A/C/G/T Ceiling** (`--n-ceil`)
   - Maximum number of non-standard bases allowed
   - Options: `None`, `L,0,0.1`, `L,0,0.15`, `L,0,0.2`, `L,0,0.3`
   - Important for handling RNA modifications

4. **Gap Barrier** (`--gbar`)
   - Prevents gaps within N nucleotides of read extremes
   - Options: `None`, `0`, `2`, `4`, `6`, `8`
   - Filters edge gaps in RNA reads

### Medium Priority

5. **Match Bonus** (`--ma`)
   - Bonus score for matches
   - Options: `None`, `0`, `1`, `2`, `3`
   - Affects overall alignment scoring

6. **Extension Effort** (`-D`)
   - How hard Bowtie2 tries to extend alignments
   - Options: `None`, `5`, `10`, `15`, `20`, `25`
   - Higher = more effort = slower but more sensitive

7. **Repetitive Seed Effort** (`-R`)
   - Number of seed sets to try for repetitive sequences
   - Options: `None`, `1`, `2`, `3`, `4`
   - Important for RNA with repetitive regions

8. **Minimum Insert Size** (`-I` / `--minins`)
   - Minimum fragment length for paired-end (only for paired-end data)
   - Options: `None`, `0`, `50`, `100`, `150`
   - Can filter out incorrect pairs

## Parameter Summary

### Total Parameters Now Optimized: **17**

**Core Parameters:**
1. seed_length
2. seed_mismatches
3. maxins
4. minins (paired-end only)

**Scoring Parameters:**
5. score_min (20+ options)
6. mismatch_penalty
7. match_bonus
8. gap_penalty_read
9. gap_penalty_ref
10. np_penalty
11. n_ceil

**Alignment Strategy:**
12. sensitivity_mode (presets)
13. seed_interval
14. extension_effort
15. repetitive_effort
16. gbar

**Performance:**
17. threads (optional)

## Parameter Space Expansion

- **Before**: ~9 parameters
- **After**: ~17 parameters
- **Score function options**: Expanded from 6 to 20+ options
- **Total combinations**: Significantly increased for more comprehensive optimization

## Usage

All parameters are automatically optimized when you run:

```bash
conda activate rna-map-optuna

python scripts/optimize_bowtie2_params_optuna.py \
    --fasta test/resources/case_1/test.fasta \
    --fastq1 test/resources/case_1/test_mate1.fastq \
    --fastq2 test/resources/case_1/test_mate2.fastq \
    --n-trials 100 \
    --output-dir optimization_results
```

## Output

All new parameters are included in:
- `optuna_study.json` - Best parameters
- `optuna_summary.csv` - All trial parameters
- `optuna_study.pkl` - Complete study object
- Best Bowtie2 arguments string

## Expected Impact

With these additions, the optimization should:
1. **Better handle RNA modifications** (np_penalty, n_ceil)
2. **Find optimal seed strategies** (seed_interval)
3. **Better gap handling** (gbar, gap penalties)
4. **More comprehensive scoring** (match_bonus)
5. **Better repeat handling** (repetitive_effort)
6. **Improved sensitivity control** (extension_effort)

The expanded parameter space should lead to better signal-to-noise ratios and more optimal parameter sets for RNA alignment!

