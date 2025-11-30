# Parameter Comparison: Original vs Optimized

## Original Baseline Parameters (Before Optimization)

```
--local --no-unal --no-discordant --no-mixed -X 1000 -L 12
```

**Only 6 parameters were set** - very minimal configuration.

### Original Parameters Explained

1. `--local` - Local alignment mode (allows clipping)
2. `--no-unal` - Don't output unaligned reads
3. `--no-discordant` - Don't output discordant pairs
4. `--no-mixed` - Don't output unpaired reads from pairs
5. `-X 1000` - Maximum insert size = 1000bp
6. `-L 12` - Seed length = 12bp

### What Was Missing

The original baseline **did not set**:
- ❌ Score minimum (used Bowtie2 default)
- ❌ Mismatch penalties (used Bowtie2 default)
- ❌ Gap penalties (used Bowtie2 default)
- ❌ Sensitivity modes
- ❌ Seed interval
- ❌ RNA-specific parameters (np, n-ceil)
- ❌ Seed mismatches
- ❌ Any effort parameters

## Optimized Parameters (After Optimization)

```
--local --no-unal --no-discordant --no-mixed -L 18 -X 200 -N 1 
--score-min L,10,0.2 --mp 6,2 --rdg 8,4 --rfg 5,3 --fast-local 
-i S,1,2.0 --np 2 --n-ceil L,0,0.3 --gbar 6 --ma 0 -D 15 -R 4 -I 50
```

**17 parameters** - comprehensive optimized configuration.

## Side-by-Side Comparison

| Parameter | Original | Optimized | Change |
|-----------|----------|-----------|--------|
| Alignment mode | `--local` | `--local` | ✅ Same |
| Filtering | `--no-unal --no-discordant --no-mixed` | `--no-unal --no-discordant --no-mixed` | ✅ Same |
| **Seed length** | `-L 12` | `-L 18` | ⬆️ **Increased** (12 → 18) |
| **Max insert size** | `-X 1000` | `-X 200` | ⬇️ **Decreased** (1000 → 200) |
| Seed mismatches | *(not set)* | `-N 1` | ✅ **Added** (allow 1 mismatch) |
| Score minimum | *(not set)* | `--score-min L,10,0.2` | ✅ **Added** |
| Mismatch penalty | *(not set)* | `--mp 6,2` | ✅ **Added** |
| Read gap penalty | *(not set)* | `--rdg 8,4` | ✅ **Added** |
| Ref gap penalty | *(not set)* | `--rfg 5,3` | ✅ **Added** |
| Sensitivity mode | *(not set)* | `--fast-local` | ✅ **Added** |
| Seed interval | *(not set)* | `-i S,1,2.0` | ✅ **Added** |
| NP penalty | *(not set)* | `--np 2` | ✅ **Added** |
| N ceiling | *(not set)* | `--n-ceil L,0,0.3` | ✅ **Added** |
| Gap barrier | *(not set)* | `--gbar 6` | ✅ **Added** |
| Match bonus | *(not set)* | `--ma 0` | ✅ **Added** |
| Extension effort | *(not set)* | `-D 15` | ✅ **Added** |
| Repetitive effort | *(not set)* | `-R 4` | ✅ **Added** |
| Min insert size | *(not set)* | `-I 50` | ✅ **Added** |

## Key Changes

### 1. Changed Parameters
- **Seed length**: 12 → **18** (longer = better specificity)
- **Max insert size**: 1000 → **200** (smaller = filters incorrect pairs)

### 2. Added Parameters (11 new)
All of these were **missing** from the original and are now **critical**:
- Score minimum function
- Mismatch penalty
- Gap penalties (read and ref)
- Sensitivity mode preset
- Seed interval
- RNA-specific parameters (np, n-ceil)
- Gap barrier
- Match bonus
- Effort parameters
- Minimum insert size
- Seed mismatches

## Impact

**Original baseline**: 6 parameters (minimal, used defaults for everything else)  
**Optimized**: 17 parameters (comprehensive, all critical parameters tuned)

The optimization process discovered that **many additional parameters** that weren't in the original baseline are actually **critical for optimal performance**!

## Files

- **Original**: `original_baseline_parameters.txt`
- **Optimized**: `best_parameters.txt`
- **Comparison**: This document

