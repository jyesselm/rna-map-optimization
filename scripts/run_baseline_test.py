#!/usr/bin/env python3
"""Run baseline test with original parameters for comparison.

This script runs a single alignment with baseline parameters and saves
the results for comparison with optimized parameters.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

# rna-map-mini package should be installed:
# pip install git+https://github.com/jyesselm/rna-map-mini.git

# Try to import from rna_map_mini
try:
    from rna_map_mini.io.fasta import fasta_to_dict
except ImportError:
    try:
        from rna_map_mini.core.fasta import fasta_to_dict
    except ImportError:
        # Fallback: implement locally or use utility
        from rna_map_optimization.utils import fasta_to_dict

try:
    from rna_map_mini.logger import get_logger
except ImportError:
    try:
        from rna_map_mini.core.logger import get_logger
    except ImportError:
        import logging
        logging.basicConfig(level=logging.INFO)
        def get_logger(name: str):
            logger = logging.getLogger(name)
            if not logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                logger.addHandler(handler)
                logger.setLevel(logging.INFO)
            return logger

try:
    from rna_map_mini.analysis.bit_vector_iterator import BitVectorIterator
    from rna_map_mini.analysis.mutation_histogram import MutationHistogram
except ImportError:
    try:
        from rna_map_mini.core.bit_vector_iterator import BitVectorIterator
        from rna_map_mini.core.mutation_histogram import MutationHistogram
    except ImportError:
        try:
            from rna_map_mini.bit_vector_iterator import BitVectorIterator
            from rna_map_mini.mutation_histogram import MutationHistogram
        except ImportError:
            raise ImportError(
                "Could not import BitVectorIterator or MutationHistogram from rna_map_mini. "
                "Please install rna-map-mini: pip install git+https://github.com/jyesselm/rna-map-mini.git"
            )

log = get_logger("BASELINE_TEST")


def get_baseline_parameters(num_sequences: int = 1) -> Dict:
    """Get baseline parameter set.
    
    Args:
        num_sequences: Number of reference sequences (1 = single, >1 = multiple)
    
    Returns:
        Dictionary with baseline parameters
    """
    base_params = {
        "local": True,
        "no_unal": True,
        "no_discordant": True,
        "no_mixed": True,
        "maxins": 1000,
    }
    
    # Always use the true baseline (single sequence parameters)
    # This is the original default: --local --no-unal --no-discordant --no-mixed -X 1000 -L 12
    base_params.update({
        "seed_length": 12,
    })
    
    return base_params


def params_to_bowtie2_args(params: Dict) -> List[str]:
    """Convert parameter dictionary to Bowtie2 argument list."""
    args = []
    
    flag_map = {
        "local": "--local",
        "no_unal": "--no-unal",
        "no_discordant": "--no-discordant",
        "no_mixed": "--no-mixed",
        "seed_length": "-L",
        "maxins": "-X",
        "score_min": "--score-min",
        "mismatch_penalty": "--mp",
        "gap_penalty_read": "--rdg",
        "gap_penalty_ref": "--rfg",
    }
    
    for key, value in params.items():
        if key in flag_map and value:
            if isinstance(value, bool) and value:
                args.append(flag_map[key])
            elif isinstance(value, (int, str)):
                args.extend([flag_map[key], str(value)])
    
    return args


def build_bowtie2_index(fasta: Path, index_dir: Path) -> Path:
    """Build Bowtie2 index from FASTA file."""
    index_dir.mkdir(parents=True, exist_ok=True)
    index_name = index_dir / fasta.stem
    
    if (index_name.with_suffix(".1.bt2")).exists():
        log.info(f"Index already exists: {index_name}")
        return index_name
    
    log.info(f"Building Bowtie2 index: {index_name}")
    cmd = ["bowtie2-build", str(fasta), str(index_name)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Failed to build index: {result.stderr}")
    
    return index_name


def run_bowtie2_alignment(
    index: Path,
    fastq1: Path,
    fastq2: Path | None,
    output_sam: Path,
    bt2_args: List[str],
    threads: int = 4,
) -> Tuple[Path, float]:
    """Run Bowtie2 alignment with given parameters."""
    stderr_file = output_sam.parent / "bowtie2_stderr.txt"
    
    cmd = ["bowtie2"] + bt2_args + ["-x", str(index), "-S", str(output_sam)]
    
    if fastq2 and fastq2.exists():
        cmd.extend(["-1", str(fastq1), "-2", str(fastq2)])
    else:
        cmd.extend(["-U", str(fastq1)])
    
    if not any(arg.startswith("-p") for arg in bt2_args):
        cmd.extend(["-p", str(threads)])
    
    start_time = time.time()
    with open(stderr_file, "w") as stderr:
        result = subprocess.run(
            cmd, stderr=stderr, stdout=subprocess.PIPE, text=True
        )
    
    elapsed_time = time.time() - start_time
    
    if result.returncode != 0:
        log.warning(f"Bowtie2 returned non-zero exit code: {result.returncode}")
    
    return stderr_file, elapsed_time


def parse_bowtie2_stats(stderr_file: Path) -> Dict:
    """Parse Bowtie2 alignment statistics from stderr output."""
    stats = {
        "total_reads": 0,
        "aligned_0_times": 0,
        "aligned_exactly_1_time": 0,
        "aligned_more_than_1_time": 0,
        "overall_alignment_rate": 0.0,
    }
    
    if not stderr_file.exists():
        return stats
    
    with open(stderr_file) as f:
        for line in f:
            line = line.strip()
            if "reads; of these:" in line:
                stats["total_reads"] = int(line.split()[0])
            elif "aligned 0 times" in line:
                stats["aligned_0_times"] = int(line.split()[0])
            elif "aligned exactly 1 time" in line:
                stats["aligned_exactly_1_time"] = int(line.split()[0])
            elif "aligned >1 times" in line:
                stats["aligned_more_than_1_time"] = int(line.split()[0])
            elif "overall alignment rate" in line:
                rate_str = line.split("%")[0].split()[-1]
                stats["overall_alignment_rate"] = float(rate_str)
    
    return stats


def parse_sam_quality(sam_file: Path, mapq_cutoff: int = 20) -> Dict:
    """Parse SAM file to calculate quality metrics."""
    metrics = {
        "total_alignments": 0,
        "high_quality_alignments": 0,
        "avg_mapq": 0.0,
        "median_mapq": 0.0,
    }
    
    if not sam_file.exists():
        return metrics
    
    mapq_scores = []
    
    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            
            metrics["total_alignments"] += 1
            mapq = int(fields[4])
            mapq_scores.append(mapq)
            
            if mapq >= mapq_cutoff:
                metrics["high_quality_alignments"] += 1
    
    if mapq_scores:
        metrics["avg_mapq"] = sum(mapq_scores) / len(mapq_scores)
        sorted_scores = sorted(mapq_scores)
        mid = len(sorted_scores) // 2
        metrics["median_mapq"] = (
            sorted_scores[mid]
            if len(sorted_scores) % 2 == 1
            else (sorted_scores[mid - 1] + sorted_scores[mid]) / 2
        )
    
    return metrics


def generate_bit_vectors_and_analyze(
    sam_file: Path,
    ref_seqs: dict[str, str],
    paired: bool,
    qscore_cutoff: int = 25,
    mapq_cutoff: int = 20,
) -> Dict:
    """Generate bit vectors and extract quality metrics."""
    metrics = {
        "total_bit_vectors": 0,
        "accepted_bit_vectors": 0,
        "rejected_low_mapq": 0,
        "signal_to_noise": 0.0,
        "constructs": {},
    }
    
    if not sam_file.exists():
        return metrics
    
    try:
        mut_histos = {}
        for ref_name, ref_seq in ref_seqs.items():
            mut_histos[ref_name] = MutationHistogram(
                ref_name, ref_seq, "DMS", 1, len(ref_seq)
            )
        
        iterator = BitVectorIterator(
            sam_file, ref_seqs, paired=paired, use_pysam=False,
            qscore_cutoff=qscore_cutoff
        )
        
        for bit_vector in iterator:
            metrics["total_bit_vectors"] += 1
            
            if not bit_vector.reads:
                continue
            
            if any(read.mapq < mapq_cutoff for read in bit_vector.reads):
                metrics["rejected_low_mapq"] += 1
                continue
            
            if not bit_vector.data:
                continue
            
            metrics["accepted_bit_vectors"] += 1
            
            ref_name = bit_vector.reads[0].rname
            if ref_name in mut_histos:
                mh = mut_histos[ref_name]
                mh.num_reads += 1
                mh.num_aligned += 1
                
                for pos in mh.get_nuc_coords():
                    if pos not in bit_vector.data:
                        continue
                    read_bit = bit_vector.data[pos]
                    mh.cov_bases[pos] += 1
                    
                    if read_bit in ["A", "C", "G", "T"]:
                        mh.mut_bases[pos] += 1
                        mh.mod_bases[read_bit][pos] += 1
        
        total_snr = 0.0
        total_weight = 0
        for ref_name, mh in mut_histos.items():
            if mh.num_aligned > 0:
                snr = mh.get_signal_to_noise()
                total_snr += snr * mh.num_aligned
                total_weight += mh.num_aligned
                
                metrics["constructs"][ref_name] = {
                    "aligned_reads": mh.num_aligned,
                    "signal_to_noise": snr,
                    "sequence_length": len(mh.ref_seq),
                }
        
        if total_weight > 0:
            metrics["signal_to_noise"] = total_snr / total_weight
        
    except Exception as e:
        log.warning(f"Error generating bit vectors: {e}")
    
    return metrics


def main():
    """Main baseline test function."""
    parser = argparse.ArgumentParser(
        description="Run baseline test with original parameters"
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to reference FASTA file",
    )
    parser.add_argument(
        "--fastq1",
        type=Path,
        required=True,
        help="Path to first FASTQ file",
    )
    parser.add_argument(
        "--fastq2",
        type=Path,
        default=None,
        help="Path to second FASTQ file (optional, for paired-end)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for alignment",
    )
    parser.add_argument(
        "--mapq-cutoff",
        type=int,
        default=20,
        help="MAPQ cutoff for high-quality alignments",
    )
    
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_dir = args.output_dir / "index"
    results_dir = args.output_dir / "results" / "baseline"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("Baseline Test with Original Parameters")
    print("=" * 80)
    print(f"FASTA: {args.fasta}")
    print(f"FASTQ1: {args.fastq1}")
    if args.fastq2:
        print(f"FASTQ2: {args.fastq2}")
    print(f"Output directory: {args.output_dir}")
    print()
    
    ref_seqs = fasta_to_dict(args.fasta)
    num_sequences = len(ref_seqs)
    
    print(f"Reference sequences found: {num_sequences}")
    for name, seq in ref_seqs.items():
        print(f"  - {name}: {len(seq)} bp")
    print()
    
    baseline_params = get_baseline_parameters(num_sequences=num_sequences)
    bt2_args = params_to_bowtie2_args(baseline_params)
    
    print("Baseline parameters:")
    for key, value in baseline_params.items():
        if value:
            print(f"  {key}: {value}")
    print()
    print(f"Bowtie2 command: bowtie2 {' '.join(bt2_args)} ...")
    print()
    
    index = build_bowtie2_index(args.fasta, index_dir)
    
    output_sam = results_dir / "aligned.sam"
    
    print("Running Bowtie2 alignment...")
    stderr_file, elapsed_time = run_bowtie2_alignment(
        index, args.fastq1, args.fastq2, output_sam, bt2_args, args.threads
    )
    
    print(f"Alignment completed in {elapsed_time:.2f} seconds")
    print()
    
    alignment_stats = parse_bowtie2_stats(stderr_file)
    quality_metrics = parse_sam_quality(output_sam, args.mapq_cutoff)
    
    print("Generating bit vectors and analyzing...")
    bit_vector_metrics = generate_bit_vectors_and_analyze(
        output_sam,
        ref_seqs,
        paired=args.fastq2 is not None,
        qscore_cutoff=25,
        mapq_cutoff=args.mapq_cutoff,
    )
    
    alignment_rate = alignment_stats.get("overall_alignment_rate", 0.0) / 100.0
    snr = bit_vector_metrics.get("signal_to_noise", 0.0)
    avg_mapq = quality_metrics.get("avg_mapq", 0.0)
    bv_acceptance_rate = (
        bit_vector_metrics.get("accepted_bit_vectors", 0)
        / bit_vector_metrics.get("total_bit_vectors", 1)
        if bit_vector_metrics.get("total_bit_vectors", 0) > 0
        else 0.0
    )
    
    normalized_snr = min(snr / 10.0, 1.0) if snr > 0 else 0.0
    quality_score = (
        0.40 * normalized_snr
        + 0.30 * alignment_rate
        + 0.20 * (avg_mapq / 60.0)
        + 0.10 * bv_acceptance_rate
    )
    
    results = {
        "parameters": baseline_params,
        "bowtie2_args": " ".join(bt2_args),
        "alignment_stats": alignment_stats,
        "quality_metrics": quality_metrics,
        "bit_vector_metrics": bit_vector_metrics,
        "metrics": {
            "signal_to_noise": snr,
            "alignment_rate": alignment_rate,
            "avg_mapq": avg_mapq,
            "quality_score": quality_score,
            "bit_vector_acceptance_rate": bv_acceptance_rate,
            "accepted_bit_vectors": bit_vector_metrics.get("accepted_bit_vectors", 0),
        },
        "elapsed_time": elapsed_time,
    }
    
    results_file = args.output_dir / "baseline_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print("=" * 80)
    print("Baseline Test Results")
    print("=" * 80)
    print(f"Signal-to-Noise (AC/GU): {snr:.2f}")
    print(f"Alignment Rate: {alignment_rate:.2%}")
    print(f"Avg MAPQ: {avg_mapq:.1f}")
    print(f"Bit Vectors Accepted: {bit_vector_metrics.get('accepted_bit_vectors', 0):,}")
    print(f"Quality Score: {quality_score:.4f}")
    print(f"Elapsed Time: {elapsed_time:.2f} seconds")
    print()
    print(f"Results saved to: {results_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

