# Bowtie2 Score Function Types

This document explains the different types of score functions available in Bowtie2 and how they're used in the optimization.

## Score Function Types

Bowtie2 supports three types of score functions for the `--score-min` parameter:

### 1. **L-type (Linear)**
- **Format**: `L,A,B`
- **Formula**: `score = A + B × read_length`
- **Use case**: Linear relationship with read length
- **Examples**:
  - `L,0,0.2`: score = 0 + 0.2 × read_length
  - `L,10,0.1`: score = 10 + 0.1 × read_length
  - `L,0,-0.6`: score = 0 - 0.6 × read_length (negative slope)

**Important**: L-type with negative values (like `L,0,-0.6`) is **incompatible** with mismatch penalties > 0. The optimization script automatically prunes these invalid combinations.

### 2. **G-type (Max)**
- **Format**: `G,A,B`
- **Formula**: `score = max(A, B × read_length)`
- **Use case**: Takes the maximum of fixed value or read-length-dependent value
- **Examples**:
  - `G,20,8`: score = max(20, 8 × read_length)
  - `G,20,15`: score = max(20, 15 × read_length)
  - `G,30,12`: score = max(30, 12 × read_length)

**Advantage**: Always positive, works with any mismatch penalty setting.

### 3. **S-type (Sum)**
- **Format**: `S,A,B`
- **Formula**: `score = A + B × read_length` (similar to L-type)
- **Use case**: Alternative to L-type, similar behavior
- **Examples**:
  - `S,20,0.15`: score = 20 + 0.15 × read_length
  - `S,15,0.2`: score = 15 + 0.2 × read_length
  - `S,25,0.2`: score = 25 + 0.2 × read_length

**Advantage**: Always positive, works with any mismatch penalty setting.

## Optimization Coverage

The optimization script now tests:

- **L-type (positive)**: 5 options
- **L-type (negative)**: 4 options (only when compatible)
- **G-type**: 6 options
- **S-type**: 5 options

**Total: 20+ score function options**

## Compatibility Rules

1. **L-type with negative values** (`L,0,-0.6`, etc.):
   - ✅ Compatible with: `mismatch_penalty=None` or `mismatch_penalty="X,0"`
   - ❌ Incompatible with: `mismatch_penalty="6,2"`, `"4,2"`, `"2,2"` (when second value > 0)

2. **L-type with positive values** (`L,0,0.2`, etc.):
   - ✅ Compatible with: All mismatch penalty settings

3. **G-type** (`G,20,15`, etc.):
   - ✅ Compatible with: All mismatch penalty settings

4. **S-type** (`S,20,0.15`, etc.):
   - ✅ Compatible with: All mismatch penalty settings

## Automatic Validation

The optimization script automatically:
1. **Tests all score function types** during optimization
2. **Validates compatibility** between score_min and mismatch_penalty
3. **Prunes invalid combinations** early (returns 0.0 for incompatible pairs)
4. **Records failure reasons** in trial attributes for debugging

## Example Results

Recent optimization runs have shown:
- **S-type functions** (`S,20,0.15`) often perform well
- **G-type functions** (`G,20,15`) are also effective
- **L-type with positive values** can be optimal in some cases
- Invalid combinations are correctly pruned (quality_score = 0.0)

## Usage in Best Parameters

The best parameters often include score functions like:
- `S,20,0.15` - S-type with moderate threshold
- `G,20,15` - G-type with read-length scaling
- `L,0,0.2` - L-type with positive slope

These provide good balance between sensitivity and specificity for RNA alignment.

