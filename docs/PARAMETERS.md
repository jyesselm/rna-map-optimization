# Optimizable Parameters

This document describes all parameters that can be optimized in the Bowtie2 parameter optimization process.

## Bowtie2 Alignment Parameters

These parameters directly control the Bowtie2 alignment behavior:

### Core Parameters

- **seed_length** (int, 10-22, step 2): Length of seed used for initial alignment search
- **maxins** (int, 200-1200, step 100): Maximum insert size for paired-end reads
- **seed_mismatches** (int, 0-1): Number of mismatches allowed in seed

### Scoring Parameters

- **score_min** (categorical, 21 options): Minimum alignment score threshold
  - Options include linear (L), function (G), and sum (S) scoring functions
  - Examples: `L,10,0.2`, `G,20,15`, `S,15,0.2`
  
- **mismatch_penalty** (categorical, 4 options): Penalty for mismatches
  - Options: null, `6,2`, `4,2`, `2,2`
  
- **gap_penalty_read** (categorical, 4 options): Gap penalties for read
  - Options: null, `5,3`, `6,4`, `8,4`
  
- **gap_penalty_ref** (categorical, 4 options): Gap penalties for reference
  - Options: null, `5,3`, `6,4`, `8,4`

### Sensitivity and Performance

- **sensitivity_mode** (categorical, 5 options): Sensitivity preset mode
  - Options: null, `very-fast-local`, `fast-local`, `sensitive-local`, `very-sensitive-local`
  
- **seed_interval** (categorical, 6 options): Seed interval function
  - Options: null, `S,1,0.5`, `S,1,0.75`, `S,1,1.15`, `S,1,1.5`, `S,1,2.0`

### Advanced Parameters

- **np_penalty** (categorical, 5 options): Penalty for non-A/C/G/T bases
  - Options: null, 0, 1, 2, 3
  
- **n_ceil** (categorical, 5 options): Ceiling on non-A/C/G/T bases
  - Options: null, `L,0,0.1`, `L,0,0.15`, `L,0,0.2`, `L,0,0.3`
  
- **gbar** (categorical, 6 options): Gap barrier parameter
  - Options: null, 0, 2, 4, 6, 8
  
- **match_bonus** (categorical, 5 options): Bonus for matches
  - Options: null, 0, 1, 2, 3
  
- **extension_effort** (categorical, 6 options): Extension effort (-D parameter)
  - Options: null, 5, 10, 15, 20, 25
  
- **repetitive_effort** (categorical, 5 options): Repetitive seed effort (-R parameter)
  - Options: null, 1, 2, 3, 4
  
- **minins** (categorical, 5 options): Minimum insert size (paired-end only)
  - Options: null, 0, 50, 100, 150

## Post-Alignment Quality Parameters

These parameters control how alignments are filtered and evaluated after alignment:

### Quality Filtering

- **mapq_cutoff** (int, 10-40, step 5): Minimum MAPQ score for high-quality alignments
  - Used to filter alignments by mapping quality
  - Lower values include more alignments but may reduce quality
  
- **qscore_cutoff** (int, 20-30, step 5): Quality score cutoff for bit vector generation
  - Controls quality threshold for base calling in bit vectors
  - Affects mutation detection sensitivity

### Quality Score Calculation

- **quality_score_weights** (categorical, 7 options): Weights for quality score components
  - Format: `[snr_weight, alignment_rate_weight, mapq_weight, bv_acceptance_weight]`
  - Must sum to 1.0
  - Default: `[0.40, 0.30, 0.20, 0.10]` (prioritize SNR)
  - Options range from SNR-focused to alignment-rate-focused

- **snr_normalization_factor** (float, 5.0-15.0, step 2.5): Factor to normalize SNR
  - Used to scale signal-to-noise ratio for quality score calculation
  - Default: 10.0
  - Affects how much SNR contributes to final score

- **min_alignment_rate** (categorical, 5 options): Minimum alignment rate threshold
  - Used for early pruning of poor alignments (0.20-0.40)
  - Default: 0.30
  - Lower values allow more trials to proceed; higher values prune earlier

## Parameter Count Summary

- **Bowtie2 alignment parameters**: 15 parameters
- **Post-alignment quality parameters**: 5 parameters
- **Total optimizable parameters**: 20 parameters

## Total Search Space

With all parameters:
- Integer ranges: 3 parameters
- Categorical parameters: 17 parameters with varying option counts
- Estimated total combinations: Millions (practically infinite)

Optuna uses Bayesian optimization to efficiently search this space, typically finding good solutions in 100-20,000 trials.

## Notes

1. **Parameter Interactions**: Many parameters interact with each other. For example, `score_min` and `mismatch_penalty` must be compatible.

2. **Early Pruning**: The `min_alignment_rate` threshold can significantly speed up optimization by pruning poor trials early.

3. **Quality Score Weights**: Different weight combinations prioritize different aspects of alignment quality. Choose based on your priorities:
   - High SNR weight: Better mutation detection
   - High alignment rate weight: More aligned reads
   - High MAPQ weight: Higher mapping confidence

4. **MAPQ vs Qscore Cutoffs**: 
   - `mapq_cutoff`: Filters entire alignments
   - `qscore_cutoff`: Filters individual bases in bit vectors
   - Both affect mutation detection but at different stages

