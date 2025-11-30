# Best Bowtie2 Parameters (From Top 100 Analysis)

These are the optimal parameters derived from analyzing the top 100 parameter combinations from `full_optimization`.

## Optimal Parameters

### Single-Line Format (for Nextflow)
```
--local --no-unal --no-discordant --no-mixed -L 18 -X 200 -N 1 --score-min L,10,0.2 --mp 6,2 --rdg 8,4 --rfg 5,3 --fast-local -i S,1,2.0 --np 2 --n-ceil L,0,0.3 --gbar 6 --ma 0 -D 15 -R 4 -I 50
```

### Multi-Line Format (for readability)
```
--local
--no-unal
--no-discordant
--no-mixed
-L 18
-X 200
-N 1
--score-min L,10,0.2
--mp 6,2
--rdg 8,4
--rfg 5,3
--fast-local
-i S,1,2.0
--np 2
--n-ceil L,0,0.3
--gbar 6
--ma 0
-D 15
-R 4
-I 50
```

## Parameter Breakdown

### Constants (100% of top 100)
- `-L 18` - Seed length
- `--mp 6,2` - Mismatch penalty
- `--rdg 8,4` - Read gap penalty
- `--rfg 5,3` - Reference gap penalty
- `--fast-local` - Sensitivity mode
- `-i S,1,2.0` - Seed interval
- `-I 50` - Minimum insert size

### High Frequency (≥95% of top 100)
- `-X 200` - Maximum insert size (98.0%)
- `-N 1` - Seed mismatches (98.0%)
- `--score-min L,10,0.2` - Score minimum (99.0%)
- `--np 2` - Non-A/C/G/T penalty (96.8%)
- `-R 4` - Repetitive seed effort (92.9%)
- `--ma 0` - Match bonus (89.7%)
- `-D 15` - Extension effort (87.1%)
- `--gbar 6` - Gap barrier (82.7%)

### Variable (most variable parameter)
- `--n-ceil L,0,0.3` - Non-A/C/G/T ceiling (66.7% - most variable)

## Usage in Nextflow

```bash
nextflow run main.nf \
    --fasta reference.fasta \
    --fastq1 reads_R1.fastq \
    --fastq2 reads_R2.fastq \
    --bt2_alignment_args "--local;--no-unal;--no-discordant;--no-mixed;-L 18;-X 200;-N 1;--score-min L,10,0.2;--mp 6,2;--rdg 8,4;--rfg 5,3;--fast-local;-i S,1,2.0;--np 2;--n-ceil L,0,0.3;--gbar 6;--ma 0;-D 15;-R 4;-I 50"
```

## Expected Results (Case 1)

Based on optimization results:
- **Quality Score**: ~0.8657
- **Signal-to-Noise**: ~8.37-8.40
- **Alignment Rate**: ~86-99%
- **Average MAPQ**: ~39-43

## Notes

- These parameters were optimized for case_1 data
- Results may vary for different datasets
- Test on your data before production use
- Consider re-optimizing for different data types

