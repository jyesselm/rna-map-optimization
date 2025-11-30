# Parameters for Optimizing Multiple Reference Sequences

When aligning reads to multiple similar reference sequences, you need stricter parameters to ensure reads align to the **best match** and distinguish between similar sequences.

## Currently Optimized ✅

1. **seed_length** (-L): 10, 12, 15, 20
2. **seed_mismatches** (-N): 0, 1, 2
3. **maxins** (-X): 200, 300, 500, 1000
4. **score_min** (--score-min): L,0,-0.6, L,0,-0.4, G,20,15
5. **mismatch_penalty** (--mp): 6,2, 4,2
6. **gap_penalty_read** (--rdg): 8,4, 5,3
7. **gap_penalty_ref** (--rfg): 5,3
8. **sensitivity_mode**: very-fast-local, fast-local, sensitive-local
9. **max_alignments** (-k): 1 (only best alignment)

## Critical Missing Parameters for Multiple Sequences

### 1. **Seed Interval** (`-i`) ⭐⭐⭐ CRITICAL
**Why**: Controls how frequently seeds are extracted from reads. More seeds = better discrimination between similar sequences.

**Options to add**:
- `S,1,0.5` (very sensitive - more seeds, better discrimination)
- `S,1,0.75` (sensitive)
- `S,1,1.15` (default)
- `S,1,1.5` (less sensitive)
- `S,1,2.0` (fast, fewer seeds)

**For multiple sequences**: Prefer lower values (0.5-0.75) for better discrimination.

### 2. **Extension Effort** (`-D`) ⭐⭐ HIGH PRIORITY
**Why**: Controls how hard Bowtie2 tries to extend alignments. More effort = better discrimination between similar sequences.

**Options to add**:
- `10` (moderate effort)
- `15` (default)
- `20` (high effort, better discrimination)
- `25` (very high effort)

**For multiple sequences**: Prefer higher values (20-25) to ensure best match is found.

### 3. **Repetitive Seed Effort** (`-R`) ⭐⭐ HIGH PRIORITY
**Why**: For reads with repetitive seeds, tries multiple seed sets. Important when sequences have similar/repetitive regions.

**Options to add**:
- `2` (default)
- `3` (more effort)
- `4` (high effort, better for similar sequences)

**For multiple sequences**: Prefer higher values (3-4) to handle repetitive regions better.

### 4. **Match Bonus** (`--ma`) ⭐ MEDIUM PRIORITY
**Why**: Bonus score for matches. Can affect which alignment is chosen as "best" when sequences are similar.

**Options to add**:
- `0` (no bonus, default for local)
- `1` (small bonus)
- `2` (moderate bonus)
- `3` (high bonus)

**For multiple sequences**: Test different values to see which helps distinguish best matches.

### 5. **Non-A/C/G/T Penalty** (`--np`) ⭐ MEDIUM PRIORITY (RNA-specific)
**Why**: Penalty for non-standard nucleotides. Important for RNA with modifications.

**Options to add**:
- `0` (tolerant)
- `1` (default)
- `2` (stricter)
- `3` (very strict)

**For multiple sequences**: Stricter values (2-3) can help distinguish between sequences.

### 6. **Non-A/C/G/T Ceiling** (`--n-ceil`) ⭐ MEDIUM PRIORITY (RNA-specific)
**Why**: Maximum number of non-standard bases allowed. Can filter ambiguous alignments.

**Options to add**:
- `L,0,0.1` (strict)
- `L,0,0.15` (default)
- `L,0,0.2` (permissive)
- `L,0,0.3` (very permissive)

**For multiple sequences**: Stricter values (0.1-0.15) can help filter poor alignments.

### 7. **Gap Barrier** (`--gbar`) ⭐ MEDIUM PRIORITY
**Why**: Prevents gaps within N nucleotides of read ends. Filters edge gaps that might be ambiguous.

**Options to add**:
- `0` (no barrier)
- `2` (small barrier)
- `4` (default)
- `6` (larger barrier)
- `8` (very large barrier)

**For multiple sequences**: Higher values (6-8) can help filter edge artifacts.

### 8. **Minimum Insert Size** (`-I`) ⭐ MEDIUM PRIORITY (Paired-end only)
**Why**: Minimum fragment length. Can filter incorrect pairs that might align ambiguously.

**Options to add**:
- `0` (no minimum)
- `50` (small minimum)
- `100` (moderate minimum)
- `150` (large minimum)

**For multiple sequences**: Moderate values (50-100) can help filter incorrect pairs.

## Recommended Parameter Ranges for Multiple Sequences

### Strict Mode (Best Discrimination)
```python
{
    "seed_length": 18,  # Longer seeds = more specific
    "seed_mismatches": 0,  # No mismatches in seed = stricter
    "seed_interval": "S,1,0.5",  # More seeds = better discrimination
    "score_min": "L,10,0.2",  # Higher threshold = prefer better matches
    "mismatch_penalty": "6,2",  # Stricter penalties
    "gap_penalty_read": "8,4",  # Stricter gap penalties
    "gap_penalty_ref": "5,3",
    "extension_effort": 25,  # High effort = find best match
    "repetitive_effort": 4,  # High effort for repeats
    "match_bonus": 0,  # Test different values
    "np_penalty": 2,  # Stricter for RNA
    "n_ceil": "L,0,0.1",  # Stricter ceiling
    "gbar": 6,  # Larger gap barrier
    "max_alignments": 1,  # Only best alignment
}
```

### Balanced Mode (Good Discrimination + Speed)
```python
{
    "seed_length": 15,  # Moderate length
    "seed_mismatches": 0,  # No seed mismatches
    "seed_interval": "S,1,0.75",  # Moderate seed frequency
    "score_min": "L,10,0.2",  # Moderate threshold
    "mismatch_penalty": "6,2",
    "gap_penalty_read": "8,4",
    "gap_penalty_ref": "5,3",
    "extension_effort": 20,  # Moderate-high effort
    "repetitive_effort": 3,  # Moderate-high effort
    "match_bonus": 0,
    "np_penalty": 2,
    "n_ceil": "L,0,0.15",
    "gbar": 4,
    "max_alignments": 1,
}
```

## Implementation Priority

### Phase 1: Add to Optimization (High Impact)
1. **Seed Interval** (`-i`) - Very important for discrimination
2. **Extension Effort** (`-D`) - Helps find best match
3. **Repetitive Seed Effort** (`-R`) - Important for similar sequences

### Phase 2: Add to Optimization (Medium Impact)
4. **Match Bonus** (`--ma`) - Affects scoring
5. **Non-A/C/G/T Penalty** (`--np`) - RNA-specific
6. **Non-A/C/G/T Ceiling** (`--n-ceil`) - RNA-specific
7. **Gap Barrier** (`--gbar`) - Filters edge artifacts

### Phase 3: Consider Adding
8. **Minimum Insert Size** (`-I`) - Paired-end only

## Key Differences: Single vs Multiple Sequences

| Parameter | Single Sequence | Multiple Sequences |
|-----------|---------------|-------------------|
| **seed_length** | 10-12 (shorter, allows mutations) | 15-20 (longer, more specific) |
| **seed_mismatches** | 0-1 (allows mutations) | 0 (strict, no seed mismatches) |
| **seed_interval** | S,1,1.5-2.0 (fewer seeds, faster) | S,1,0.5-0.75 (more seeds, better discrimination) |
| **score_min** | L,0,-0.6 (permissive, allows mutations) | L,10,0.2 (stricter, prefers better matches) |
| **extension_effort** | 10-15 (moderate) | 20-25 (high, find best match) |
| **repetitive_effort** | 2 (default) | 3-4 (high, handle similar sequences) |
| **max_alignments** | 1 (or more for mutation tracking) | 1 (only best match) |

## Expected Impact

Adding these parameters should:
1. **Reduce multi-mapping**: Better discrimination between sequences
2. **Improve MAPQ scores**: More confident alignments to best match
3. **Increase specificity**: Fewer ambiguous alignments
4. **Better signal-to-noise**: Clearer distinction between sequences

## Testing Strategy

1. **Start with Phase 1 parameters** (seed_interval, extension_effort, repetitive_effort)
2. **Run optimization** on case_2 (9 sequences)
3. **Compare multi-mapping rates** - should be 0% with `-k 1`
4. **Compare MAPQ distributions** - should be higher with better parameters
5. **Add Phase 2 parameters** if needed for further improvement

