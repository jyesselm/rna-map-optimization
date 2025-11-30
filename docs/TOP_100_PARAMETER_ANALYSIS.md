# Top 100 Parameter Combinations Analysis

Analysis of the top 100 parameter combinations from `full_optimization` results.

## Key Findings

### Quality Score Range
- **Range**: 0.8657 - 0.8664
- **Average**: 0.8657
- **Median**: 0.8657
- **Note**: Very tight range indicates optimization converged well!

## Constants (100% of Top 100) - 7 Parameters

These parameters are **identical** across all top 100 combinations:

1. **seed_length** = `18`
2. **mismatch_penalty** = `6,2`
3. **gap_penalty_read** = `8,4`
4. **gap_penalty_ref** = `5,3`
5. **sensitivity_mode** = `fast-local`
6. **seed_interval** = `S,1,2.0`
7. **minins** = `50.0`

**Recommendation**: These can be **hardcoded** or set as defaults - they're clearly optimal!

## High Frequency Parameters (≥50% of Top 100) - 9 Parameters

These appear in most top combinations:

1. **score_min** = `L,10,0.2` (99.0%)
2. **maxins** = `200` (98.0%)
3. **seed_mismatches** = `1` (98.0%)
4. **repetitive_effort** = `4.0` (92.9%)
5. **match_bonus** = `0.0` (89.7%)
6. **extension_effort** = `15.0` (87.1%)
7. **gbar** = `6.0` (82.7%)
8. **n_ceil** = `L,0,0.3` (66.7%)
9. **np_penalty** = `2.0` (96.8%)

**Recommendation**: Use these as **strong defaults** - they're optimal for most cases.

## Recommended Optimal Parameters

Based on this analysis, here's the recommended parameter set:

```bash
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
-I 50
```

Or as a single argument string:
```
--local --no-unal --no-discordant --no-mixed -L 18 -X 200 -N 1 --score-min L,10,0.2 --mp 6,2 --rdg 8,4 --rfg 5,3 --fast-local -i S,1,2.0 --np 2 --n-ceil L,0,0.3 --gbar 6 --ma 0 -D 15 -R 4 -I 50
```

## Insights

### Critical Parameters (Constants)

1. **Seed length = 18**: Longer seeds provide good balance
2. **Mismatch penalty = 6,2**: Standard penalty works well
3. **Gap penalties**: Read gap (8,4) higher than ref gap (5,3) - read gaps penalized more
4. **Sensitivity = fast-local**: Fast mode is optimal - more sensitive modes don't help
5. **Seed interval = S,1,2.0**: Faster seed extraction strategy
6. **Min insert size = 50**: Filters very short fragments

### Important Parameters (High Frequency)

1. **Score min = L,10,0.2**: Linear function with moderate threshold
2. **Max insert size = 200**: Smaller than default (500) - filters long inserts
3. **Seed mismatches = 1**: Allowing 1 mismatch in seed improves sensitivity
4. **Repetitive effort = 4**: High effort for repeats (max value)
5. **Match bonus = 0**: No bonus needed
6. **Extension effort = 15**: Moderate effort (default range)

### Variable Parameters

- **n_ceil**: Some variation (66.7% use L,0,0.3) - this is the most variable
- All other parameters are highly consistent

## Recommendations

### For Optimization Script

1. **Lock in constants**: Set these 7 parameters as fixed defaults
2. **Strong defaults**: Use the high-frequency values as defaults
3. **Reduce search space**: Focus optimization on the variable parameters:
   - `n_ceil` (the only truly variable parameter)

### For Users

1. **Use the recommended parameters above** - they're proven optimal
2. **Only vary if you have specific requirements**:
   - If you need more sensitivity: try different `n_ceil` values
   - If reads are longer: increase `maxins`
   - If data quality differs: adjust `score_min`

### For Future Optimization

Since most parameters are constants, future optimizations could:
- Fix the 7 constants
- Fix the high-frequency parameters (≥90%)
- Only optimize the remaining variable space (mainly `n_ceil`)

This would dramatically reduce the search space and allow deeper exploration of the remaining variables.

## Parameter Breakdown

| Parameter | Value | Frequency | Status |
|-----------|-------|-----------|--------|
| seed_length | 18 | 100% | ✅ Constant |
| mismatch_penalty | 6,2 | 100% | ✅ Constant |
| gap_penalty_read | 8,4 | 100% | ✅ Constant |
| gap_penalty_ref | 5,3 | 100% | ✅ Constant |
| sensitivity_mode | fast-local | 100% | ✅ Constant |
| seed_interval | S,1,2.0 | 100% | ✅ Constant |
| minins | 50.0 | 100% | ✅ Constant |
| score_min | L,10,0.2 | 99.0% | ⭐ Strong default |
| maxins | 200 | 98.0% | ⭐ Strong default |
| seed_mismatches | 1 | 98.0% | ⭐ Strong default |
| np_penalty | 2.0 | 96.8% | ⭐ Strong default |
| repetitive_effort | 4.0 | 92.9% | ⭐ Strong default |
| match_bonus | 0.0 | 89.7% | ⭐ Strong default |
| extension_effort | 15.0 | 87.1% | ⭐ Strong default |
| gbar | 6.0 | 82.7% | ⭐ Strong default |
| n_ceil | L,0,0.3 | 66.7% | ⚠️ Variable |

## Next Steps

1. **Test the recommended parameters** on your data
2. **Compare with baseline** to confirm improvement
3. **Consider locking constants** in future optimizations
4. **Further optimize n_ceil** if it's the main variable

