# Additional Bowtie2 Parameters for Optimization

This document lists all potentially relevant Bowtie2 parameters that could be optimized beyond what's currently included.

## Currently Optimized Parameters ✅

1. **seed_length** (-L): 10-22 (step 2)
2. **seed_mismatches** (-N): 0-1
3. **maxins** (-X): 200-1200 (step 100)
4. **score_min** (--score-min): 20+ options (L, G, S types)
5. **mismatch_penalty** (--mp): 4 options
6. **gap_penalty_read** (--rdg): 4 options
7. **gap_penalty_ref** (--rfg): 4 options
8. **sensitivity_mode**: 5 preset options
9. **threads**: 2-16 (optional)

## High-Priority Additional Parameters

### 1. **Seed Interval** (`-i`) - HIGH PRIORITY
- **What it does**: Interval between seed substrings relative to read length
- **Format**: Function like `S,1,1.15` (type,min,max)
- **Impact**: Controls how frequently seeds are extracted from reads
- **Relevance**: Very important for alignment sensitivity
- **Options to test**:
  - `S,1,0.5` (very sensitive - more seeds)
  - `S,1,0.75` (sensitive)
  - `S,1,1.15` (default)
  - `S,1,1.5` (less sensitive - fewer seeds)
  - `S,1,2.0` (very fast - few seeds)
- **Note**: Lower values = more seeds = slower but more sensitive

### 2. **Non-A/C/G/T Penalty** (`--np`) - HIGH PRIORITY FOR RNA
- **What it does**: Penalty for non-standard nucleotides (N, I, etc.)
- **Format**: Integer (default: 1)
- **Impact**: Important for RNA data with modified nucleotides
- **Relevance**: Very relevant for DMS-MaP data with modified bases
- **Options to test**: `0`, `1`, `2`, `3`
- **Note**: Lower values = more tolerant of non-standard bases

### 3. **Non-A/C/G/T Ceiling** (`--n-ceil`) - HIGH PRIORITY FOR RNA
- **What it does**: Maximum number of non-A/C/G/Ts permitted in alignment
- **Format**: Function like `L,0,0.15` (type,min,max)
- **Impact**: Controls how many non-standard bases are allowed
- **Relevance**: Very important for RNA with modifications
- **Options to test**:
  - `L,0,0.1` (strict)
  - `L,0,0.15` (default)
  - `L,0,0.2` (permissive)
  - `L,0,0.3` (very permissive)

### 4. **Gap Barrier** (`--gbar`) - MEDIUM PRIORITY
- **What it does**: Disallow gaps within N nucleotides of read extremes
- **Format**: Integer (default: 4)
- **Impact**: Prevents gaps at read ends (often low-quality regions)
- **Relevance**: Good for filtering edge gaps in RNA reads
- **Options to test**: `0`, `2`, `4`, `6`, `8`
- **Note**: Higher = fewer gaps near ends

### 5. **Match Bonus** (`--ma`) - MEDIUM PRIORITY
- **What it does**: Bonus score for matches
- **Format**: Integer (0 for end-to-end, 2 for local, default varies)
- **Impact**: Affects overall alignment scoring
- **Relevance**: Can affect alignment quality
- **Options to test**: `0`, `1`, `2`, `3`
- **Note**: Higher = rewards matches more

### 6. **Extension Effort** (`-D`) - MEDIUM PRIORITY
- **What it does**: Give up extending after N failed extends in a row
- **Format**: Integer (default: 15, preset-dependent)
- **Impact**: Controls how hard Bowtie2 tries to extend alignments
- **Relevance**: Affects sensitivity vs speed
- **Options to test**: `5`, `10`, `15`, `20`, `25`
- **Note**: Higher = more effort = slower but more sensitive

### 7. **Repetitive Seed Effort** (`-R`) - MEDIUM PRIORITY
- **What it does**: For reads with repetitive seeds, try N sets of seeds
- **Format**: Integer (default: 2, preset-dependent)
- **Impact**: Handles repetitive sequences better
- **Relevance**: Important for RNA with repetitive regions
- **Options to test**: `1`, `2`, `3`, `4`
- **Note**: Higher = better handling of repeats but slower

### 8. **DP Table Padding** (`--dpad`) - LOW PRIORITY
- **What it does**: Include N extra reference characters on sides of DP table
- **Format**: Integer (default: 15)
- **Impact**: Affects dynamic programming table size
- **Relevance**: Mostly performance, minimal quality impact
- **Options to test**: `10`, `15`, `20`
- **Note**: Higher = more memory, minimal quality gain

## Paired-End Specific Parameters

### 9. **Minimum Insert Size** (`-I` / `--minins`) - MEDIUM PRIORITY
- **What it does**: Minimum fragment length for paired-end
- **Format**: Integer (default: 0)
- **Impact**: Filters pairs that are too close together
- **Relevance**: Can filter out PCR duplicates or incorrect pairs
- **Options to test**: `0`, `50`, `100`, `150`
- **Note**: Only relevant for paired-end data

### 10. **Paired-End Orientation** (`--fr`/`--rf`/`--ff`) - LOW PRIORITY
- **What it does**: Controls expected orientation of paired reads
- **Format**: Flag (default: `--fr` for forward/reverse)
- **Impact**: Affects which alignments are considered concordant
- **Relevance**: Usually known from library prep, but could optimize
- **Options to test**: `--fr`, `--rf`, `--ff`
- **Note**: Usually determined by library type

### 11. **Concordance Options** - LOW PRIORITY
- **--dovetail**: Allow mates extending past each other
- **--no-contain**: Not concordant when one mate contains the other
- **--no-overlap**: Not concordant when mates overlap
- **Relevance**: Edge cases, usually not needed

## Alignment Strategy Parameters

### 12. **Alignment Mode** (`--local` vs `--end-to-end`) - HIGH PRIORITY
- **What it does**: Local (clipped) vs end-to-end alignment
- **Impact**: Major difference in alignment behavior
- **Relevance**: Currently using `--local`, but could test `--end-to-end`
- **Note**: Would require significant changes to current setup

### 13. **One-MM-Upfront** (`--no-1mm-upfront`) - MEDIUM PRIORITY
- **What it does**: Disable allowing 1 mismatch alignments before optimal seeded alignment scan
- **Impact**: Changes alignment strategy/order
- **Relevance**: Could affect speed/quality tradeoff
- **Options to test**: Include or exclude this flag

## Quality-Related Parameters

### 14. **Ignore Qualities** (`--ignore-quals`) - LOW PRIORITY
- **What it does**: Treat all quality values as 30 (Phred)
- **Impact**: Ignores quality scores completely
- **Relevance**: Usually not recommended, but could test
- **Note**: Generally decreases quality

### 15. **Quality Encoding** (`--phred33`/`--phred64`) - LOW PRIORITY
- **What it does**: Specify quality score encoding
- **Impact**: Correct interpretation of quality scores
- **Relevance**: Usually auto-detected, but could be explicit
- **Note**: Should match your data format

## Read Filtering Parameters

### 16. **Directional Alignment** (`--nofw`/`--norc`) - LOW PRIORITY
- **What it does**: Don't align forward or reverse-complement version
- **Impact**: Limits alignment direction
- **Relevance**: RNA-seq usually bidirectional, but could test
- **Note**: Usually not relevant for RNA

## Recommendations for Implementation

### Phase 1: High-Impact Additions
1. **Seed Interval** (`-i`) - Very important for sensitivity
2. **Non-A/C/G/T Penalty** (`--np`) - Important for RNA modifications
3. **Non-A/C/G/T Ceiling** (`--n-ceil`) - Important for RNA modifications
4. **Gap Barrier** (`--gbar`) - Good for RNA edge quality

### Phase 2: Medium-Impact Additions
5. **Match Bonus** (`--ma`) - Scoring parameter
6. **Extension Effort** (`-D`) - Sensitivity vs speed
7. **Repetitive Seed Effort** (`-R`) - Repeat handling
8. **Minimum Insert Size** (`-I`) - Paired-end filtering

### Phase 3: Low-Impact or Edge Cases
9. **DP Table Padding** (`--dpad`)
10. **Paired-end orientation**
11. **Other alignment strategy flags**

## Parameter Interactions

Be aware of these interactions:
- **Sensitivity presets** automatically set `-D`, `-R`, `-i`, `-L`, `-N`
- If using presets, many parameters are already set
- **Score functions** interact with mismatch penalties (already handled)
- **Gap penalties** work together (read + ref)

## Implementation Considerations

1. **Parameter space explosion**: Adding more parameters increases combinations exponentially
2. **Optuna efficiency**: Bayesian optimization handles large spaces well, but too many parameters can slow convergence
3. **Validation**: Need to ensure parameter combinations are valid
4. **Default conflicts**: Presets override many parameters - may need to choose preset vs individual params

## Suggested Next Steps

1. **Start with seed interval** (`-i`) - high impact, easy to add
2. **Add RNA-specific parameters** (`--np`, `--n-ceil`) - important for your use case
3. **Add gap barrier** (`--gbar`) - quick win for RNA
4. **Consider match bonus** (`--ma`) - simple scoring adjustment
5. **Evaluate effort parameters** (`-D`, `-R`) - balance sensitivity vs speed

## Example: Expanded Parameter Set

If implementing all high-priority parameters:

```python
# Seed interval
seed_interval = trial.suggest_categorical("seed_interval", [
    "S,1,0.5", "S,1,0.75", "S,1,1.15", "S,1,1.5", "S,1,2.0"
])

# Non-A/C/G/T penalty
np_penalty = trial.suggest_int("np_penalty", 0, 3)

# Non-A/C/G/T ceiling
n_ceil = trial.suggest_categorical("n_ceil", [
    "L,0,0.1", "L,0,0.15", "L,0,0.2", "L,0,0.3"
])

# Gap barrier
gbar = trial.suggest_int("gbar", 0, 8, step=2)

# Match bonus
match_bonus = trial.suggest_int("match_bonus", 0, 3)

# Extension effort
extension_effort = trial.suggest_int("extension_effort", 5, 25, step=5)
```

This would add ~6 more parameters, significantly expanding the search space but potentially finding better optima.

