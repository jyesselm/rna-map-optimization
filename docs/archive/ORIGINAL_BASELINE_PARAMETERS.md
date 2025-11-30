# Original/Baseline Bowtie2 Parameters

These are the **original default parameters** used before optimization.

## Original Baseline Parameters

### From `conf/base.config` (Nextflow defaults)
```
--local --no-unal --no-discordant --no-mixed -X 1000 -L 12
```

### From `test/resources/default.yml` (Python defaults)
```
--local --no-unal --no-discordant --no-mixed -X 1000 -L 12 -p 16
```

### Expanded Baseline (from optimization script)
The baseline test uses:
- `--local`
- `--no-unal`
- `--no-discordant`
- `--no-mixed`
- `-L 12` (seed length)
- `-X 1000` (maximum insert size)
- `-p 16` (threads, if specified)

**No additional parameters were set:**
- No score_min
- No mismatch_penalty
- No gap penalties
- No seed interval
- No sensitivity mode presets
- No other advanced parameters

## Comparison: Original vs Optimized

| Parameter | Original (Baseline) | Optimized (Top 100) | Change |
|-----------|-------------------|---------------------|--------|
| **seed_length** | `-L 12` | `-L 18` | ⬆️ Longer seeds |
| **maxins** | `-X 1000` | `-X 200` | ⬇️ Smaller inserts |
| **seed_mismatches** | Not set (default 0) | `-N 1` | ⬆️ Allow 1 mismatch |
| **score_min** | Not set | `--score-min L,10,0.2` | ✅ Added |
| **mismatch_penalty** | Not set | `--mp 6,2` | ✅ Added |
| **gap_penalty_read** | Not set | `--rdg 8,4` | ✅ Added |
| **gap_penalty_ref** | Not set | `--rfg 5,3` | ✅ Added |
| **sensitivity_mode** | Not set | `--fast-local` | ✅ Added |
| **seed_interval** | Not set | `-i S,1,2.0` | ✅ Added |
| **np_penalty** | Not set | `--np 2` | ✅ Added |
| **n_ceil** | Not set | `--n-ceil L,0,0.3` | ✅ Added |
| **gbar** | Not set | `--gbar 6` | ✅ Added |
| **match_bonus** | Not set | `--ma 0` | ✅ Added |
| **extension_effort** | Not set | `-D 15` | ✅ Added |
| **repetitive_effort** | Not set | `-R 4` | ✅ Added |
| **minins** | Not set | `-I 50` | ✅ Added |

## Original Parameters Summary

### Minimal Set (6 parameters)
1. `--local` - Local alignment mode
2. `--no-unal` - Suppress unaligned reads
3. `--no-discordant` - Suppress discordant alignments
4. `--no-mixed` - Suppress mixed alignments
5. `-L 12` - Seed length = 12
6. `-X 1000` - Maximum insert size = 1000

### What Was Missing

The original baseline was **very minimal** - it only set:
- Alignment mode (local)
- Read filtering flags
- Seed length
- Maximum insert size

**Missing important parameters:**
- ❌ Score functions
- ❌ Mismatch penalties
- ❌ Gap penalties
- ❌ Sensitivity modes
- ❌ Seed interval
- ❌ RNA-specific parameters (np, n-ceil)
- ❌ Effort parameters

## Key Improvements from Optimization

### 1. Seed Length: 12 → 18
- Longer seeds improve alignment quality
- **Impact**: Better specificity

### 2. Max Insert Size: 1000 → 200
- Smaller inserts filter out incorrect pairs
- **Impact**: Better paired-end filtering

### 3. Added Score Minimum
- `--score-min L,10,0.2` filters low-quality alignments
- **Impact**: Improved alignment quality

### 4. Added Gap Penalties
- Read gap: `8,4` (higher penalty)
- Ref gap: `5,3` (lower penalty)
- **Impact**: Better RNA alignment handling

### 5. Added RNA-Specific Parameters
- `--np 2` - Non-standard nucleotide penalty
- `--n-ceil L,0,0.3` - Non-standard nucleotide ceiling
- **Impact**: Better handling of RNA modifications

### 6. Added Sensitivity Control
- `--fast-local` mode
- Seed interval optimization
- Effort parameters
- **Impact**: Better balance of speed and sensitivity

## Original Baseline Command

```bash
bowtie2 \
    --local \
    --no-unal \
    --no-discordant \
    --no-mixed \
    -L 12 \
    -X 1000 \
    -x index \
    -1 reads_R1.fastq \
    -2 reads_R2.fastq \
    -S output.sam
```

## Optimized Command (for comparison)

```bash
bowtie2 \
    --local \
    --no-unal \
    --no-discordant \
    --no-mixed \
    -L 18 \
    -X 200 \
    -N 1 \
    --score-min L,10,0.2 \
    --mp 6,2 \
    --rdg 8,4 \
    --rfg 5,3 \
    --fast-local \
    -i S,1,2.0 \
    --np 2 \
    --n-ceil L,0,0.3 \
    --gbar 6 \
    --ma 0 \
    -D 15 \
    -R 4 \
    -I 50 \
    -x index \
    -1 reads_R1.fastq \
    -2 reads_R2.fastq \
    -S output.sam
```

## Impact Summary

The optimization added **11 additional parameters** that weren't in the original baseline:
- Most importantly: score functions, gap penalties, and RNA-specific parameters
- These additions significantly improved signal-to-noise ratio

**Original**: 6 parameters (minimal)  
**Optimized**: 17 parameters (comprehensive)

The optimization process discovered that many additional parameters are crucial for optimal RNA alignment quality!

