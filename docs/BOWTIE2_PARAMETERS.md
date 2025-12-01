# Bowtie2 Parameter Reference

This document explains what each Bowtie2 parameter does and how it affects alignment behavior. These parameters are used in the optimization process to find the best settings for RNA-MAP experiments.

## Overview

Bowtie2 is a fast and sensitive aligner for sequencing reads. The parameters control:
- **Sensitivity vs. Speed**: How thoroughly to search for alignments
- **Stringency**: How strict to be about mismatches and gaps
- **Seed-based alignment**: Initial search strategy
- **Scoring**: How to evaluate alignment quality

## Fixed Parameters (Always Used)

These parameters are always set in the optimization and are not varied:

- **`--local`**: Local alignment mode (allows soft-clipping at read ends)
  - Essential for RNA-MAP because reads may not align end-to-end
  - Allows partial alignments which is important for mutation detection

- **`--no-unal`**: Don't report unaligned reads in SAM output
  - Reduces output file size
  - Only aligned reads are needed for analysis

- **`--no-discordant`**: Don't report discordant paired-end alignments
  - Discordant = pairs that don't align in expected orientation/distance
  - Filters out likely mis-alignments

- **`--no-mixed`**: Don't report mixed alignments (one end aligned, one not)
  - Ensures both ends of paired reads align
  - Important for paired-end analysis

## Core Alignment Parameters

### `seed_length` (-L)

**Type**: Integer (10-22, step 2)  
**Default**: 22  
**What it does**: Length of the seed (initial exact match) used to start alignment search.

**Effect**:
- **Longer seeds (18-22)**: Faster, more specific, fewer false positives
  - Better for: Multiple reference sequences, high-quality data
  - Trade-off: May miss alignments with mutations in seed region
- **Shorter seeds (10-14)**: Slower, more sensitive, finds more alignments
  - Better for: Single sequence, mutation detection, lower quality data
  - Trade-off: More false positives, slower alignment

**For RNA-MAP**: Shorter seeds (10-14) are often better because they allow mutations in the seed region, which is important for detecting DMS modifications.

### `maxins` (-X)

**Type**: Integer (200-1200, step 100)  
**Default**: 500  
**What it does**: Maximum insert size for paired-end reads (distance between R1 and R2).

**Effect**:
- **Smaller values (200-400)**: Stricter pairing, fewer discordant pairs
  - Better for: Short RNA fragments (~150bp), precise fragment size
- **Larger values (800-1200)**: More permissive, allows longer fragments
  - Better for: Variable fragment sizes, degraded RNA

**For RNA-MAP**: Typically 200-500bp works well for ~150bp RNA fragments.

### `seed_mismatches` (-N)

**Type**: Integer (0-1)  
**Default**: 0  
**What it does**: Number of mismatches allowed in the seed.

**Effect**:
- **0 mismatches**: Strict seed matching, faster, fewer false positives
- **1 mismatch**: More sensitive, finds alignments with mutations in seed, slower

**For RNA-MAP**: 1 mismatch is often beneficial to detect mutations in the seed region.

## Scoring Parameters

### `score_min` (--score-min)

**Type**: Categorical (21 options)  
**Default**: None (uses Bowtie2 defaults)  
**What it does**: Minimum alignment score threshold. Format: `TYPE,MIN,FRAC`

**Types**:
- **L (Linear)**: `L,MIN,FRAC` - Score must be ≥ MIN + FRAC × read_length
  - Examples: `L,10,0.2` means score ≥ 10 + 0.2×read_length
  - Negative values (e.g., `L,0,-0.6`) allow lower scores (more permissive)
- **G (Function)**: `G,MIN,FRAC` - Score must be ≥ MIN + FRAC × (read_length - MIN)
- **S (Sum)**: `S,MIN,FRAC` - Score must be ≥ MIN + FRAC × read_length

**Effect**:
- **Higher thresholds** (e.g., `L,15,0.1`): Stricter, fewer alignments, higher quality
- **Lower/negative thresholds** (e.g., `L,0,-0.6`): More permissive, more alignments, allows mutations

**For RNA-MAP**: Negative or low thresholds (e.g., `L,0,-0.6`) are often needed to allow mutations while maintaining alignment quality.

### `mismatch_penalty` (--mp)

**Type**: Categorical (4 options)  
**Default**: None (uses Bowtie2 default: 6,2)  
**What it does**: Penalty for mismatches. Format: `MAX,MIN`

- `MAX`: Penalty for mismatches with highest-quality bases
- `MIN`: Penalty for mismatches with lowest-quality bases

**Options**:
- `6,2`: Default, strict (penalizes mismatches heavily)
- `4,2`: Moderate (allows more mismatches)
- `2,2`: Lenient (allows many mismatches)

**Effect**:
- **Stricter penalties** (`6,2`): Fewer mismatches allowed, higher alignment quality
- **Lenient penalties** (`2,2`): More mismatches allowed, better for mutation detection

**For RNA-MAP**: Moderate penalties (`4,2` or `6,2`) balance quality with mutation detection.

**Note**: Incompatible with negative `score_min` values - if using `L,0,-0.6`, don't use strict mismatch penalties.

### `gap_penalty_read` (--rdg) and `gap_penalty_ref` (--rfg)

**Type**: Categorical (4 options each)  
**Default**: None  
**What it does**: Penalties for gaps (insertions/deletions). Format: `OPEN,EXTEND`

- `OPEN`: Penalty for opening a gap
- `EXTEND`: Penalty for extending a gap

**Options**:
- `8,4`: Strict (fewer gaps allowed)
- `6,4`: Moderate
- `5,3`: Lenient (allows more gaps)

**Effect**:
- **Stricter penalties**: Fewer gaps, cleaner alignments
- **Lenient penalties**: More gaps allowed, better for sequences with insertions/deletions

**For RNA-MAP**: Moderate penalties (`5,3` or `6,4`) work well for most cases.

## Sensitivity Modes

### `sensitivity_mode`

**Type**: Categorical (5 options)  
**Default**: None (uses default local mode)  
**What it does**: Preset sensitivity/performance trade-offs.

**Options**:
- `very-fast-local`: Fastest, least sensitive
- `fast-local`: Fast, moderate sensitivity (often good for RNA-MAP)
- `sensitive-local`: Balanced (default for many cases)
- `very-sensitive-local`: Slowest, most sensitive

**Effect**:
- **Faster modes**: Quick alignment, may miss some alignments
- **Sensitive modes**: Thorough search, finds more alignments, slower

**For RNA-MAP**: `fast-local` or `sensitive-local` are good starting points.

## Seed Interval

### `seed_interval` (-i)

**Type**: Categorical (6 options)  
**Default**: None  
**What it does**: Controls how frequently seeds are extracted from the read. Format: `FUNC,START,INC`

- `FUNC`: Function type (usually `S` for step)
- `START`: Starting position
- `INC`: Increment (lower = more seeds)

**Options**:
- `S,1,0.5`: Many seeds (every 0.5×read_length)
- `S,1,0.75`: Moderate seeds
- `S,1,1.15`: Default (fewer seeds)
- `S,1,1.5`: Few seeds
- `S,1,2.0`: Very few seeds

**Effect**:
- **More seeds** (`S,1,0.5`): Better for multiple references, more thorough search, slower
- **Fewer seeds** (`S,1,2.0`): Faster, fewer seeds to check

**For RNA-MAP**: More seeds (`S,1,0.5-0.75`) help when aligning to multiple reference sequences.

## Advanced Parameters

### `np_penalty` (--np)

**Type**: Categorical (5 options)  
**Default**: None (uses Bowtie2 default: 1)  
**What it does**: Penalty for non-A/C/G/T bases (N, ambiguous bases).

**Options**: 0, 1, 2, 3

**Effect**:
- **Lower penalty** (0-1): Allows more N bases in alignments
- **Higher penalty** (2-3): Stricter about N bases

**For RNA-MAP**: Usually not critical, default (1) is fine.

### `n_ceil` (--n-ceil)

**Type**: Categorical (5 options)  
**Default**: None  
**What it does**: Ceiling on number of N bases allowed. Format: `L,0,FRAC`

**Options**: `L,0,0.1`, `L,0,0.15`, `L,0,0.2`, `L,0,0.3`

**Effect**:
- **Lower ceiling** (`L,0,0.1`): Stricter, fewer N bases allowed
- **Higher ceiling** (`L,0,0.3`): More permissive

**For RNA-MAP**: Usually not critical unless reads have many N bases.

### `gbar` (--gbar)

**Type**: Categorical (6 options)  
**Default**: None (uses Bowtie2 default: 4)  
**What it does**: Gap barrier - prevents gaps within this many bases of read ends.

**Options**: 0, 2, 4, 6, 8

**Effect**:
- **Higher values** (6-8): Prevents gaps near read ends (cleaner alignments)
- **Lower values** (0-2): Allows gaps anywhere

**For RNA-MAP**: Default (4) or higher (6-8) helps with read quality near ends.

### `match_bonus` (--ma)

**Type**: Categorical (5 options)  
**Default**: None (uses Bowtie2 default: 2)  
**What it does**: Bonus score for matches.

**Options**: 0, 1, 2, 3

**Effect**:
- **Higher bonus** (2-3): Rewards matches more, prefers high-identity alignments
- **Lower bonus** (0-1): Less reward for matches

**For RNA-MAP**: Default (2) is usually fine.

### `extension_effort` (-D)

**Type**: Categorical (6 options)  
**Default**: None (uses Bowtie2 default: 15)  
**What it does**: Maximum number of times to extend alignment attempts.

**Options**: 5, 10, 15, 20, 25

**Effect**:
- **Higher effort** (20-25): More thorough search, finds more alignments, slower
- **Lower effort** (5-10): Faster, may miss some alignments

**For RNA-MAP**: Moderate to high effort (15-25) helps find best alignments, especially with multiple references.

### `repetitive_effort` (-R)

**Type**: Categorical (5 options)  
**Default**: None (uses Bowtie2 default: 2)  
**What it does**: Maximum number of times to re-seed when encountering repetitive sequences.

**Options**: 1, 2, 3, 4

**Effect**:
- **Higher effort** (3-4): Better handling of repetitive sequences, slower
- **Lower effort** (1-2): Faster, may have issues with repetitive regions

**For RNA-MAP**: Higher effort (3-4) helps when aligning to multiple similar reference sequences.

### `minins` (-I)

**Type**: Categorical (5 options, paired-end only)  
**Default**: None  
**What it does**: Minimum insert size for paired-end reads.

**Options**: 0, 50, 100, 150

**Effect**:
- **Higher minimum** (100-150): Stricter pairing, requires larger fragments
- **Lower minimum** (0-50): More permissive, allows shorter fragments

**For RNA-MAP**: 50-100bp works well for typical RNA fragments.

## Quality Filtering Parameters (Post-Alignment)

These parameters are used after alignment to filter results:

### `mapq_cutoff`

**Type**: Integer (10-40, step 5)  
**Default**: 20  
**What it does**: Minimum MAPQ (mapping quality) score for accepting alignments.

**MAPQ Scale**: 0-60 (higher = more confident)

**Effect**:
- **Lower cutoff** (10-15): Includes more alignments, may have lower quality
- **Higher cutoff** (25-40): Stricter, only high-confidence alignments

**For RNA-MAP**: 20-25 is a good balance between quality and sensitivity.

### `qscore_cutoff`

**Type**: Integer (20-30, step 5)  
**Default**: 25  
**What it does**: Quality score cutoff for base calling in bit vector generation.

**Quality Score Scale**: 0-40 (Phred scale)

**Effect**:
- **Lower cutoff** (20): Includes more bases, may have more errors
- **Higher cutoff** (30): Stricter, only high-quality bases

**For RNA-MAP**: 25 is standard, balances quality and coverage.

## Parameter Interactions

### Important Combinations

1. **Mutation Detection** (single sequence):
   - Short seed (`seed_length: 10-12`)
   - Allow seed mismatches (`seed_mismatches: 1`)
   - Permissive score (`score_min: L,0,-0.6`)
   - Moderate mismatch penalty (`mismatch_penalty: 4,2`)

2. **Multiple Sequences** (discrimination):
   - Longer seed (`seed_length: 18-20`)
   - No seed mismatches (`seed_mismatches: 0`)
   - Stricter score (`score_min: L,10,0.2`)
   - Stricter mismatch penalty (`mismatch_penalty: 6,2`)
   - More seeds (`seed_interval: S,1,0.5`)

3. **Speed vs. Quality**:
   - Fast: `very-fast-local`, fewer seeds, lower effort
   - Quality: `sensitive-local`, more seeds, higher effort

### Incompatible Combinations

- **Negative score_min + Strict mismatch penalty**: Don't use `L,0,-0.6` with `mismatch_penalty: 6,2` (they conflict)
- **Very short seed + No mismatches**: `seed_length: 10` with `seed_mismatches: 0` may be too restrictive

## Optimization Strategy

The optimization process tests different combinations to find the best balance:

1. **Signal-to-Noise Ratio**: Maximize AC/GU ratio (mutation signal)
2. **Alignment Rate**: Maintain good coverage
3. **MAPQ Quality**: Ensure high-confidence alignments
4. **Bit Vector Acceptance**: Maximize usable reads

The best parameters depend on:
- Number of reference sequences (single vs. multiple)
- Read quality
- Expected mutation rate
- Fragment length

## Practical Examples

### Example 1: Single Sequence with Mutations

```bash
bowtie2 --local --no-unal --no-discordant --no-mixed \
    -L 10 -N 1 -X 300 \
    --score-min L,0,-0.6 \
    --mp 4,2 \
    -x index -1 R1.fastq -2 R2.fastq -S output.sam
```

**Why these parameters?**
- Short seed (`-L 10`) allows mutations in seed region
- Seed mismatches (`-N 1`) permits mutations
- Permissive score (`L,0,-0.6`) allows lower-scoring alignments with mutations
- Moderate mismatch penalty (`4,2`) balances quality and mutation detection

### Example 2: Multiple Reference Sequences

```bash
bowtie2 --local --no-unal --no-discordant --no-mixed \
    -L 18 -N 0 -X 300 \
    --score-min L,10,0.2 \
    --mp 6,2 \
    --rdg 8,4 --rfg 5,3 \
    -i S,1,0.5 -D 25 -R 4 \
    --sensitive-local \
    -x index -1 R1.fastq -2 R2.fastq -S output.sam
```

**Why these parameters?**
- Longer seed (`-L 18`) for better discrimination
- No seed mismatches (`-N 0`) for strict matching
- Stricter score (`L,10,0.2`) prefers better matches
- More seeds (`-i S,1,0.5`) for thorough search
- High effort (`-D 25 -R 4`) to find best alignment

### Example 3: Fast Alignment (Lower Quality Acceptable)

```bash
bowtie2 --local --no-unal --no-discordant --no-mixed \
    -L 12 -X 500 \
    --very-fast-local \
    -x index -1 R1.fastq -2 R2.fastq -S output.sam
```

**Why these parameters?**
- Medium seed length for balance
- Fast preset for speed
- Minimal other parameters for simplicity

## Common Parameter Patterns

### Pattern 1: Maximum Sensitivity
- `seed_length: 10`
- `seed_mismatches: 1`
- `score_min: L,0,-0.6`
- `mismatch_penalty: 2,2`
- `sensitivity_mode: very-sensitive-local`
- `extension_effort: 25`
- `repetitive_effort: 4`

### Pattern 2: Maximum Speed
- `seed_length: 22`
- `seed_mismatches: 0`
- `sensitivity_mode: very-fast-local`
- `extension_effort: 5`
- `repetitive_effort: 1`

### Pattern 3: Balanced (Recommended Starting Point)
- `seed_length: 12-15`
- `seed_mismatches: 0-1`
- `score_min: L,0,-0.4` or `L,10,0.2`
- `mismatch_penalty: 6,2` or `4,2`
- `sensitivity_mode: fast-local` or `sensitive-local`
- `extension_effort: 15-20`
- `repetitive_effort: 2-3`

## References

- [Bowtie2 Manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [Bowtie2 Paper](https://www.nature.com/articles/nmeth.1923)
- [Bowtie2 GitHub](https://github.com/BenLangmead/bowtie2)

## Quick Reference Table

| Parameter | Flag | Range | Effect |
|-----------|------|-------|--------|
| `seed_length` | `-L` | 10-22 | Shorter = more sensitive, longer = faster |
| `maxins` | `-X` | 200-1200 | Max paired-end insert size |
| `seed_mismatches` | `-N` | 0-1 | Allow mutations in seed |
| `score_min` | `--score-min` | Various | Alignment score threshold |
| `mismatch_penalty` | `--mp` | `6,2`, `4,2`, `2,2` | Penalty for mismatches |
| `gap_penalty_read` | `--rdg` | `5,3`, `6,4`, `8,4` | Penalty for gaps in read |
| `gap_penalty_ref` | `--rfg` | `5,3`, `6,4`, `8,4` | Penalty for gaps in reference |
| `sensitivity_mode` | `--fast-local`, etc. | 5 options | Speed vs. sensitivity preset |
| `seed_interval` | `-i` | `S,1,0.5` to `S,1,2.0` | Seed extraction frequency |
| `extension_effort` | `-D` | 5-25 | Alignment extension attempts |
| `repetitive_effort` | `-R` | 1-4 | Re-seeding for repetitive sequences |
| `minins` | `-I` | 0-150 | Min paired-end insert size |

