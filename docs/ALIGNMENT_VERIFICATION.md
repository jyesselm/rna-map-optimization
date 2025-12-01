# Alignment Verification Tools

This document describes tools and methods for verifying Bowtie2 alignments as a secondary check. These are useful for validating alignment quality on small subsamples.

## Available Verification Methods

### 1. Biopython Pairwise Alignment (Recommended for Small Samples)

**Pros:**
- Pure Python, easy to install
- Exact alignment algorithms (Needleman-Wunsch, Smith-Waterman)
- Detailed alignment information
- No external dependencies

**Cons:**
- Slow for large datasets (but fine for small subsamples)
- Memory intensive for very long sequences

**Installation:**
```bash
pip install biopython
```

**Usage:**
```bash
python scripts/verify_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --fastq reads.fastq \
    --method biopython \
    --max-reads 100
```

**What it does:**
- Performs local pairwise alignment for each read
- Calculates identity, mismatches, gaps
- Compares with Bowtie2 alignment results
- Reports disagreements and alignment quality

### 2. Parasail (SIMD-Accelerated)

**Pros:**
- Fast (SIMD-accelerated)
- Still accurate
- Good for medium-sized samples

**Cons:**
- Requires compilation (but pip installs pre-built wheels)
- Less detailed than Biopython

**Installation:**
```bash
pip install parasail
```

**Usage:**
```bash
python scripts/verify_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --fastq reads.fastq \
    --method parasail \
    --max-reads 500
```

### 3. Minimap2 (Alternative Aligner)

**Pros:**
- Very fast
- Designed for short reads
- Can align entire dataset quickly
- Industry-standard tool

**Cons:**
- External binary dependency
- Different algorithm (may give different results)

**Installation:**
```bash
conda install -c bioconda minimap2
# or
brew install minimap2  # macOS
```

**Usage:**
```bash
python scripts/verify_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --fastq reads.fastq \
    --method minimap2
```

**Direct usage:**
```bash
# Align with minimap2
minimap2 -x sr reference.fasta reads.fastq > minimap2_aligned.sam

# Compare with Bowtie2 results
```

### 4. BLAST (Nucleotide BLAST)

**Pros:**
- Very sensitive
- Industry standard
- Detailed statistics

**Cons:**
- Slowest option
- Requires database building
- Overkill for most cases

**Installation:**
```bash
conda install -c bioconda blast
```

**Usage:**
```bash
# Build BLAST database
makeblastdb -in reference.fasta -dbtype nucl -out ref_db

# Align reads
blastn -query reads.fasta -db ref_db -outfmt 6 > blast_results.txt
```

## Verification Script

The `scripts/verify_alignments.py` script provides a unified interface for verification:

### Features

1. **Read-by-read verification**: Checks each aligned read individually
2. **Identity calculation**: Computes alignment identity
3. **Disagreement detection**: Finds cases where methods disagree
4. **Statistics**: Reports agreement rates, identity distributions

### Example Output

```
Verification Results:
  Total reads checked: 100
  Aligned in Bowtie2: 95
  Aligned in biopython: 93
  Agreement: 91
  Disagreement: 2
  Bowtie2 only: 2
  Avg identity: 94.5%
  Min identity: 87.2%
  Max identity: 99.8%

  Disagreements found: 2
  First 5 disagreements:
    - read_001: MAPQ=15, Identity=65.3%
    - read_042: MAPQ=20, Identity=72.1%
```

## When to Use Verification

### Recommended Use Cases

1. **Parameter Validation**: After optimization, verify best parameters on a small subset
2. **Quality Control**: Check alignment quality on representative samples
3. **Troubleshooting**: Investigate low-quality alignments
4. **Method Comparison**: Compare Bowtie2 with other aligners

### Sample Size Recommendations

- **Biopython**: 10-100 reads (slow but thorough)
- **Parasail**: 100-1000 reads (balanced)
- **Minimap2**: Full dataset (fast enough for verification)

## Interpretation

### Good Agreement

- **Identity > 90%**: Excellent alignment
- **Agreement rate > 95%**: Methods agree, alignment is reliable
- **Low disagreement**: Bowtie2 parameters are appropriate

### Potential Issues

- **Low identity (< 80%)**: May indicate:
  - Poor quality reads
  - Wrong reference sequence
  - Too permissive alignment parameters
  
- **High disagreement**: May indicate:
  - Suboptimal Bowtie2 parameters
  - Different alignment strategies needed
  - Reference sequence issues

### Action Items

1. **If identity is low**: Check read quality, verify reference sequence
2. **If disagreement is high**: Consider adjusting Bowtie2 parameters
3. **If many reads don't align**: Check if they should align (verify with minimap2/BLAST)

## Integration with Optimization

You can add verification as a post-processing step:

```python
# After optimization, verify best parameters
best_sam = "optimization_results/best_trial/aligned.sam"

# Verify on small subset
subprocess.run([
    "python", "scripts/verify_alignments.py",
    "--sam", best_sam,
    "--fasta", "reference.fasta",
    "--fastq", "reads.fastq",
    "--method", "biopython",
    "--max-reads", "50",
    "--output", "verification_results.json"
])
```

## Alternative Tools

### For Very Small Samples (< 10 reads)

**Manual inspection** with tools like:
- [NCBI BLAST web interface](https://blast.ncbi.nlm.nih.gov/)
- [EMBOSS Water](http://emboss.sourceforge.net/apps/cvs/emboss/apps/water.html)
- [Pairwise alignment tools](https://www.ebi.ac.uk/Tools/psa/)

### For Medium Samples (100-1000 reads)

- **Minimap2**: Fast and accurate
- **Parasail**: Good balance
- **LAST**: Another fast aligner

### For Large Samples (> 1000 reads)

- **Minimap2**: Recommended
- **BWA-MEM**: Alternative aligner
- **STAR**: For RNA-seq (overkill for simple verification)

## Performance Comparison

| Method | Speed | Accuracy | Best For |
|--------|-------|----------|----------|
| Biopython | Slow | Very High | < 100 reads |
| Parasail | Medium | High | 100-1000 reads |
| Minimap2 | Fast | High | > 1000 reads |
| BLAST | Very Slow | Very High | < 10 reads |

## Example Workflow

```bash
# 1. Run optimization
python scripts/optimize_bowtie2_params_optuna.py \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --n-trials 50 \
    --output-dir optimization_results

# 2. Get best alignment
BEST_SAM="optimization_results/results/trial_0/aligned.sam"

# 3. Verify on small subset
python scripts/verify_alignments.py \
    --sam "$BEST_SAM" \
    --fasta reference.fasta \
    --fastq reads_R1.fastq \
    --method biopython \
    --max-reads 50 \
    --output verification_results.json

# 4. Review results
cat verification_results.json | python -m json.tool
```

## References

- [Biopython Pairwise Alignment](https://biopython.org/docs/1.75/api/Bio.pairwise2.html)
- [Parasail Documentation](https://github.com/jeffdaily/parasail)
- [Minimap2 Manual](https://lh3.github.io/minimap2/minimap2.html)
- [BLAST Documentation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs)

