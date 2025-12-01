# Viewing Alignments in Human-Readable Format

The `view_alignments.py` script provides a beautiful, human-readable visualization of SAM alignments, making it easy to inspect how reads align to reference sequences.

## Features

- **Visual alignment display**: See reference and read sequences aligned side-by-side
- **Match indicators**: `|` for matches, `X` for mismatches
- **Color coding**: Green for matches, red for mismatches, yellow for gaps
- **Quality scores**: Display quality scores below alignments
- **Statistics**: Identity percentage, match/mismatch counts, insertions/deletions
- **Position markers**: Reference positions shown every 10 bases
- **HTML export**: Save colored output as HTML file

## Installation

The script works with or without the `rich` library:

```bash
# Basic (plain text output)
# No additional installation needed

# Enhanced (colored output)
pip install rich
```

## Basic Usage

### View First 10 Alignments

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta
```

### View Specific Number of Reads

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 5
```

### Filter by MAPQ Score

```bash
# Only show high-quality alignments (MAPQ >= 20)
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --min-mapq 20
```

### View Specific Read

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --read-name "FS10000899:22:BPG61606-0731:1:1101:12450:1000"
```

### Save to HTML File

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 10 \
    --output alignments.html
```

### Hide Quality Scores

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --no-quality
```

## Output Format

The script displays alignments in the following format:

```
╭──────────────────────────────────────────────────────────────────────────────╮
│ Read: read_name | Reference: ref_name | Position: 1 | MAPQ: 23              │
╰──────────────────────────────────────────────────────────────────────────────╯
Ref       1 CGGA------ 11    AGATCGAGTA 21    GATCAAAGGA 31    CGTATGGCGG 41    
       ||XX             XXXXXXXXXX       XXXX||X|XX       XXXXX|XX||       
Read     CGATCTGGAA       GATCGAGTAG       ATCAAAGGAC       GTATGGCGGG       
Qual    1775117177       5775151775       5711771771       7151157177       
Identity: 30.8% | Matches: 41 | Mismatches: 92 | Insertions: 1 | Deletions: 6
```

### Understanding the Display

1. **Header**: Shows read name, reference name, alignment position, and MAPQ score
2. **Ref line**: Reference sequence with position numbers every 10 bases
3. **Match line**: `|` for matches, `X` for mismatches, spaces for gaps/clips
4. **Read line**: Read sequence aligned to reference
5. **Qual line**: Quality scores (0-9, where 7-9 = high quality)
6. **Statistics**: Summary of alignment quality

### Color Coding (with `rich` installed)

- **Green**: Matches (correct alignments)
- **Red**: Mismatches (substitutions)
- **Yellow**: Gaps (insertions/deletions)
- **Dim**: Soft-clipped bases (not aligned)

## Use Cases

### 1. Quality Control

Check alignment quality on a sample of reads:

```bash
python scripts/view_alignments.py \
    --sam optimization_results/best_trial/aligned.sam \
    --fasta reference.fasta \
    --max-reads 20 \
    --min-mapq 20 \
    --output qc_check.html
```

### 2. Troubleshooting Low-Quality Alignments

Investigate why certain reads have low MAPQ:

```bash
# View low-quality alignments
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 10 \
    --min-mapq 0 \
    --output low_quality.html
```

### 3. Parameter Validation

After optimization, verify that alignments look correct:

```bash
python scripts/view_alignments.py \
    --sam optimized_aligned.sam \
    --fasta reference.fasta \
    --max-reads 50 \
    --min-mapq 20
```

### 4. Compare Different Alignments

View alignments from different parameter sets side-by-side:

```bash
# Baseline
python scripts/view_alignments.py \
    --sam baseline/aligned.sam \
    --fasta reference.fasta \
    --max-reads 5 \
    --output baseline_view.html

# Optimized
python scripts/view_alignments.py \
    --sam optimized/aligned.sam \
    --fasta reference.fasta \
    --max-reads 5 \
    --output optimized_view.html
```

## Integration with Verification

Combine with the verification script for comprehensive quality checks:

```bash
# 1. Verify alignments
python scripts/verify_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --fastq reads.fastq \
    --method biopython \
    --max-reads 50

# 2. View specific reads that had issues
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 10 \
    --min-mapq 15  # Lower threshold to see problematic reads
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--sam` | Path to SAM file (required) | - |
| `--fasta` | Path to reference FASTA (required) | - |
| `--max-reads` | Maximum alignments to display | 10 |
| `--min-mapq` | Minimum MAPQ to display | 0 |
| `--read-name` | Display specific read by name | - |
| `--no-quality` | Don't show quality scores | False |
| `--output` | Save output to file (HTML or text) | - |

## Tips

1. **Start small**: Use `--max-reads 5` to get a quick overview
2. **Filter by quality**: Use `--min-mapq 20` to focus on high-quality alignments
3. **Save HTML**: Use `--output file.html` to save colored output for sharing
4. **Combine with grep**: Find specific reads first, then view them:
   ```bash
   grep "read_name" aligned.sam | head -1 > temp.sam
   python scripts/view_alignments.py --sam temp.sam --fasta reference.fasta
   ```

## Limitations

- **Large files**: For very large SAM files, use `--max-reads` to limit output
- **CIGAR parsing**: Complex CIGAR strings with many operations may display slowly
- **Memory**: Very long reads (>1000bp) may be truncated in display

## Examples

### Example 1: Quick Quality Check

```bash
python scripts/view_alignments.py \
    --sam test_verification/test_aligned.sam \
    --fasta test_cases/case_1/test.fasta \
    --max-reads 3 \
    --min-mapq 20
```

### Example 2: Detailed Analysis

```bash
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 50 \
    --min-mapq 20 \
    --output detailed_analysis.html
```

### Example 3: Find Problematic Reads

```bash
# View reads with low MAPQ
python scripts/view_alignments.py \
    --sam aligned.sam \
    --fasta reference.fasta \
    --max-reads 20 \
    --min-mapq 0 \
    --output low_mapq_reads.html
```

## Related Tools

- **`verify_alignments.py`**: Verify alignments using alternative methods
- **`run_baseline_test.py`**: Run baseline alignment tests
- **`optimize_bowtie2_params_optuna.py`**: Optimize alignment parameters

## Troubleshooting

### No alignments displayed

- Check that SAM file has alignments (not just headers)
- Lower `--min-mapq` threshold
- Verify reference FASTA matches SAM file references

### Alignment looks wrong

- Check that reference sequence in FASTA matches SAM reference name
- Verify CIGAR string is valid
- Check if read is reverse complemented (flag in SAM)

### Slow performance

- Reduce `--max-reads` value
- Use `--min-mapq` to filter low-quality alignments
- For very large files, consider using `samtools view` to extract subset first

