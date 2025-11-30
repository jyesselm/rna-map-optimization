#!/usr/bin/env python3
"""Optimize Bowtie2 parameters using Optuna for intelligent hyperparameter search.

This script uses Optuna (Bayesian optimization) to efficiently search the parameter
space and find optimal Bowtie2 parameters that maximize signal-to-noise ratio.
More efficient than exhaustive grid search - intelligently explores parameter space.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

try:
    import optuna
    from optuna.visualization import (
        plot_optimization_history,
        plot_param_importances,
        plot_parallel_coordinate,
    )
except ImportError:
    print("ERROR: Optuna is not installed. Install with: pip install optuna")
    sys.exit(1)

# rna-map package should be installed from main repository:
# pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python
# Or if published to PyPI: pip install rna-map

from rna_map.io.fasta import fasta_to_dict
from rna_map.logger import get_logger
from rna_map.analysis.bit_vector_iterator import BitVectorIterator
from rna_map.analysis.mutation_histogram import MutationHistogram

log = get_logger("BOWTIE2_OPTIMIZER_OPTUNA")


def parse_bowtie2_stats(stderr_file: Path) -> Dict:
    """Parse Bowtie2 alignment statistics from stderr output.
    
    Args:
        stderr_file: Path to Bowtie2 stderr output
        
    Returns:
        Dictionary with alignment statistics
    """
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
    """Parse SAM file to calculate quality metrics.
    
    Args:
        sam_file: Path to SAM file
        mapq_cutoff: Minimum MAPQ score for high-quality alignments
        
    Returns:
        Dictionary with quality metrics
    """
    metrics = {
        "total_alignments": 0,
        "high_quality_alignments": 0,
        "avg_mapq": 0.0,
        "median_mapq": 0.0,
        "mapq_distribution": {},
    }
    
    if not sam_file.exists():
        return metrics
    
    mapq_scores = []
    mapq_counts = {}
    
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
            mapq_counts[mapq] = mapq_counts.get(mapq, 0) + 1
            
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
        metrics["mapq_distribution"] = mapq_counts
    
    return metrics


def build_bowtie2_index(fasta: Path, index_dir: Path) -> Path:
    """Build Bowtie2 index from FASTA file.
    
    Args:
        fasta: Path to FASTA file
        index_dir: Directory to store index
        
    Returns:
        Path to index basename
    """
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
    """Run Bowtie2 alignment with given parameters.
    
    Args:
        index: Path to Bowtie2 index
        fastq1: Path to first FASTQ file
        fastq2: Path to second FASTQ file (optional)
        output_sam: Path to output SAM file
        bt2_args: List of Bowtie2 arguments
        threads: Number of threads
        
    Returns:
        Tuple of (stderr_file, elapsed_time)
    """
    stderr_file = output_sam.parent / "bowtie2_stderr.txt"
    
    # Build command
    cmd = ["bowtie2"] + bt2_args + ["-x", str(index), "-S", str(output_sam)]
    
    if fastq2 and fastq2.exists():
        cmd.extend(["-1", str(fastq1), "-2", str(fastq2)])
    else:
        cmd.extend(["-U", str(fastq1)])
    
    # Add threads if not already specified
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


def generate_bit_vectors_and_analyze(
    sam_file: Path,
    ref_seqs: dict[str, str],
    paired: bool,
    qscore_cutoff: int = 25,
    mapq_cutoff: int = 20,
) -> Dict:
    """Generate bit vectors and extract quality metrics including signal-to-noise.
    
    Args:
        sam_file: Path to SAM file
        ref_seqs: Dictionary of reference sequences
        paired: Whether reads are paired-end
        qscore_cutoff: Quality score cutoff for bit vector generation
        mapq_cutoff: MAPQ cutoff for filtering
        
    Returns:
        Dictionary with bit vector metrics including signal-to-noise ratio
        and per-construct statistics
    """
    from collections import defaultdict
    
    metrics = {
        "total_bit_vectors": 0,
        "accepted_bit_vectors": 0,
        "rejected_low_mapq": 0,
        "mutation_distribution": defaultdict(int),
        "avg_mutations_per_read": 0.0,
        "reads_with_mutations": 0,
        "reads_without_mutations": 0,
        "total_mutations": 0,
        "coverage_positions": 0,
        "avg_coverage": 0.0,
        "signal_to_noise": 0.0,  # AC/GU ratio from mutation histogram
        "constructs": {},  # Per-construct statistics
    }
    
    if not sam_file.exists():
        return metrics
    
    try:
        # Create mutation histograms for each reference
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
            
            # Check if reads pass MAPQ filter
            if not bit_vector.reads:
                continue
            
            # Filter by MAPQ
            if any(read.mapq < mapq_cutoff for read in bit_vector.reads):
                metrics["rejected_low_mapq"] += 1
                continue
            
            if not bit_vector.data:
                continue
            
            metrics["accepted_bit_vectors"] += 1
            
            # Update mutation histogram
            ref_name = bit_vector.reads[0].rname
            if ref_name in mut_histos:
                mh = mut_histos[ref_name]
                mh.num_reads += 1
                mh.num_aligned += 1
                
                # Update mutation histogram with bit vector data
                total_muts = 0
                for pos in mh.get_nuc_coords():
                    if pos not in bit_vector.data:
                        continue
                    read_bit = bit_vector.data[pos]
                    mh.cov_bases[pos] += 1
                    
                    if read_bit in ["A", "C", "G", "T"]:
                        mh.mut_bases[pos] += 1
                        mh.mod_bases[read_bit][pos] += 1
                        total_muts += 1
                    elif read_bit == "1":
                        mh.del_bases[pos] += 1
                    elif read_bit == "?":
                        mh.info_bases[pos] += 1
                
                mh.num_of_mutations[total_muts] += 1
            
            # Count mutations for distribution
            mut_count = sum(
                1 for v in bit_vector.data.values() 
                if v in ["A", "C", "G", "T"]
            )
            metrics["mutation_distribution"][mut_count] += 1
            metrics["total_mutations"] += mut_count
            
            if mut_count > 0:
                metrics["reads_with_mutations"] += 1
            else:
                metrics["reads_without_mutations"] += 1
            
            # Count coverage
            metrics["coverage_positions"] += len(bit_vector.data)
        
        # Calculate signal-to-noise from mutation histograms
        # Average across all references (weighted by number of aligned reads)
        # Also store per-construct statistics
        total_snr = 0.0
        total_weight = 0
        for ref_name, mh in mut_histos.items():
            if mh.num_aligned > 0:
                snr = mh.get_signal_to_noise()
                total_snr += snr * mh.num_aligned
                total_weight += mh.num_aligned
                
                # Store per-construct statistics
                metrics["constructs"][ref_name] = {
                    "aligned_reads": mh.num_aligned,
                    "signal_to_noise": snr,
                    "sequence_length": len(mh.ref_seq),
                    "total_reads": mh.num_reads,
                }
        
        if total_weight > 0:
            metrics["signal_to_noise"] = total_snr / total_weight
        
        # Calculate averages
        if metrics["accepted_bit_vectors"] > 0:
            metrics["avg_mutations_per_read"] = (
                metrics["total_mutations"] / metrics["accepted_bit_vectors"]
            )
            metrics["avg_coverage"] = (
                metrics["coverage_positions"] / metrics["accepted_bit_vectors"]
            )
        
        # Convert defaultdict to regular dict
        metrics["mutation_distribution"] = dict(metrics["mutation_distribution"])
        
    except Exception as e:
        log.warning(f"Error generating bit vectors: {e}")
    
    return metrics


def params_to_bowtie2_args(params: Dict) -> List[str]:
    """Convert parameter dictionary to Bowtie2 argument list.
    
    Args:
        params: Parameter dictionary
        
    Returns:
        List of Bowtie2 arguments
    """
    args = []
    
    # Mapping from param names to Bowtie2 flags
    flag_map = {
        "local": "--local",
        "end_to_end": "--end-to-end",
        "no_unal": "--no-unal",
        "no_discordant": "--no-discordant",
        "no_mixed": "--no-mixed",
        "seed_length": "-L",
        "seed_mismatches": "-N",
        "maxins": "-X",
        "minins": "-I",
        "score_min": "--score-min",
        "mismatch_penalty": "--mp",
        "gap_penalty_read": "--rdg",
        "gap_penalty_ref": "--rfg",
        "seed_interval": "-i",
        "np_penalty": "--np",
        "n_ceil": "--n-ceil",
        "gbar": "--gbar",
        "match_bonus": "--ma",
        "extension_effort": "-D",
        "repetitive_effort": "-R",
        "very_fast_local": "--very-fast-local",
        "fast_local": "--fast-local",
        "sensitive_local": "--sensitive-local",
        "very_sensitive_local": "--very-sensitive-local",
    }
    
    for key, value in params.items():
        if key in flag_map and value:
            if isinstance(value, bool) and value:
                args.append(flag_map[key])
            elif isinstance(value, (int, str)):
                args.extend([flag_map[key], str(value)])
    
    return args


def objective(
    trial: optuna.Trial,
    index: Path,
    fastq1: Path,
    fastq2: Path | None,
    ref_seqs: dict[str, str],
    output_dir: Path,
    mapq_cutoff: int,
    threads: int,
    optimize_threads: bool,
) -> float:
    """Optuna objective function to maximize.
    
    Args:
        trial: Optuna trial object
        index: Path to Bowtie2 index
        fastq1: Path to first FASTQ file
        fastq2: Path to second FASTQ file (optional)
        ref_seqs: Dictionary of reference sequences
        output_dir: Output directory
        mapq_cutoff: MAPQ cutoff
        threads: Number of threads (if not optimizing)
        optimize_threads: Whether to optimize threads
        
    Returns:
        Quality score to maximize
    """
    # Core alignment parameters
    seed_length = trial.suggest_int("seed_length", 10, 22, step=2)
    maxins = trial.suggest_int("maxins", 200, 1200, step=100)
    
    # Seed mismatches (0 is most common, constrain to 0-1)
    seed_mismatches = trial.suggest_int("seed_mismatches", 0, 1)
    
    # Score minimum options (expanded - all three types: L, G, S)
    # L-type (Linear): score = A + B*read_length
    # G-type (Max): score = max(A, B*read_length)  
    # S-type (Sum): score = A + B*read_length (similar to L)
    # Note: L-type with negative values incompatible with mismatch penalty > 0
    # (will be validated later)
    score_min_options = [
        None,
        # L-type (Linear) - positive values (compatible with any mismatch penalty)
        "L,0,0.2",
        "L,0,0.3",
        "L,5,0.1",
        "L,10,0.2",
        "L,15,0.1",
        # L-type with negative values (only compatible when mismatch penalty allows)
        "L,0,-0.6",
        "L,0,-0.4",
        "L,0,-0.2",
        "L,-0.6,-0.6",
        # G-type (Max) - always positive, works with any mismatch penalty
        "G,15,0.1",
        "G,20,8",
        "G,20,15",
        "G,25,10",
        "G,30,12",
        "G,30,15",
        # S-type (Sum) - positive values
        "S,10,0.1",
        "S,15,0.2",
        "S,20,0.15",
        "S,25,0.2",
        "S,30,0.1",
    ]
    
    score_min_choice = trial.suggest_categorical("score_min", score_min_options)
    
    # Mismatch penalty options (expanded)
    mp_choice = trial.suggest_categorical(
        "mismatch_penalty",
        [None, "6,2", "4,2", "2,2"]
    )
    
    # Gap penalties (important for RNA alignments)
    gap_read_choice = trial.suggest_categorical(
        "gap_penalty_read",
        [None, "5,3", "6,4", "8,4"]
    )
    gap_ref_choice = trial.suggest_categorical(
        "gap_penalty_ref",
        [None, "5,3", "6,4", "8,4"]
    )
    
    # Sensitivity mode (expanded)
    sens_mode = trial.suggest_categorical(
        "sensitivity_mode",
        [
            None,
            "very-fast-local",
            "fast-local",
            "sensitive-local",
            "very-sensitive-local",
        ],
    )
    
    # Seed interval (HIGH PRIORITY - controls seed extraction frequency)
    seed_interval_choice = trial.suggest_categorical(
        "seed_interval",
        [None, "S,1,0.5", "S,1,0.75", "S,1,1.15", "S,1,1.5", "S,1,2.0"]
    )
    
    # Non-A/C/G/T penalty (HIGH PRIORITY FOR RNA - modified nucleotides)
    np_penalty_choice = trial.suggest_categorical(
        "np_penalty",
        [None, 0, 1, 2, 3]
    )
    
    # Non-A/C/G/T ceiling (HIGH PRIORITY FOR RNA - max non-standard bases)
    n_ceil_choice = trial.suggest_categorical(
        "n_ceil",
        [None, "L,0,0.1", "L,0,0.15", "L,0,0.2", "L,0,0.3"]
    )
    
    # Gap barrier (MEDIUM PRIORITY - prevents gaps near read ends)
    gbar_choice = trial.suggest_categorical(
        "gbar",
        [None, 0, 2, 4, 6, 8]
    )
    
    # Match bonus (MEDIUM PRIORITY - scoring parameter)
    match_bonus_choice = trial.suggest_categorical(
        "match_bonus",
        [None, 0, 1, 2, 3]
    )
    
    # Extension effort (MEDIUM PRIORITY - how hard to try extending)
    extension_effort_choice = trial.suggest_categorical(
        "extension_effort",
        [None, 5, 10, 15, 20, 25]
    )
    
    # Repetitive seed effort (MEDIUM PRIORITY - repeat handling)
    repetitive_effort_choice = trial.suggest_categorical(
        "repetitive_effort",
        [None, 1, 2, 3, 4]
    )
    
    # Minimum insert size (MEDIUM PRIORITY - for paired-end)
    minins_choice = None
    if fastq2 is not None:
        minins_choice = trial.suggest_categorical(
            "minins",
            [None, 0, 50, 100, 150]
        )
    
    # Threads (optimize if requested, but note: more threads = faster, not better quality)
    trial_threads = threads
    if optimize_threads:
        trial_threads = trial.suggest_int("threads", 2, min(threads * 2, 16))
    
    # Build parameter dictionary
    params = {
        "local": True,
        "no_unal": True,
        "no_discordant": True,
        "no_mixed": True,
        "seed_length": seed_length,
        "maxins": maxins,
    }
    
    # Only add seed_mismatches if > 0 (reduces parameter space)
    if seed_mismatches > 0:
        params["seed_mismatches"] = seed_mismatches
    
    # Validate score_min and mismatch_penalty compatibility
    if score_min_choice and mp_choice:
        # L-type with negative values incompatible with mismatch penalty > 0
        if score_min_choice.startswith("L,"):
            score_parts = score_min_choice.split(",")
            if len(score_parts) >= 3:
                try:
                    # Check if any part is negative
                    if float(score_parts[1]) < 0 or float(score_parts[2]) < 0:
                        mp_parts = mp_choice.split(",")
                        if len(mp_parts) >= 2 and int(mp_parts[1]) > 0:
                            # Incompatible: return 0 to prune
                            trial.set_user_attr("failure_reason", "incompatible_score_mismatch")
                            return 0.0
                except (ValueError, IndexError):
                    pass  # If parsing fails, let Bowtie2 handle it
    
    if score_min_choice:
        params["score_min"] = score_min_choice
    if mp_choice:
        params["mismatch_penalty"] = mp_choice
    if gap_read_choice:
        params["gap_penalty_read"] = gap_read_choice
    if gap_ref_choice:
        params["gap_penalty_ref"] = gap_ref_choice
    if sens_mode:
        params[sens_mode] = True
    if seed_interval_choice:
        params["seed_interval"] = seed_interval_choice
    if np_penalty_choice is not None:
        params["np_penalty"] = np_penalty_choice
    if n_ceil_choice:
        params["n_ceil"] = n_ceil_choice
    if gbar_choice is not None:
        params["gbar"] = gbar_choice
    if match_bonus_choice is not None:
        params["match_bonus"] = match_bonus_choice
    if extension_effort_choice is not None:
        params["extension_effort"] = extension_effort_choice
    if repetitive_effort_choice is not None:
        params["repetitive_effort"] = repetitive_effort_choice
    if minins_choice is not None:
        params["minins"] = minins_choice
    
    # Convert to Bowtie2 arguments
    bt2_args = params_to_bowtie2_args(params)
    
    # Create temporary output directory for this trial
    trial_dir = output_dir / f"trial_{trial.number}"
    trial_dir.mkdir(exist_ok=True)
    output_sam = trial_dir / "aligned.sam"
    
    try:
        # Run alignment
        stderr_file, _ = run_bowtie2_alignment(
            index, fastq1, fastq2, output_sam, bt2_args, trial_threads
        )
        
        # Parse statistics
        alignment_stats = parse_bowtie2_stats(stderr_file)
        quality_metrics = parse_sam_quality(output_sam, mapq_cutoff)
        
        # Early pruning: check if alignment was successful
        if alignment_stats["total_reads"] == 0:
            trial.set_user_attr("failure_reason", "no_reads")
            return 0.0  # Prune trials with no reads
        
        alignment_rate = alignment_stats.get("overall_alignment_rate", 0.0) / 100.0
        
        # Early pruning: low alignment rate
        if alignment_rate < 0.3:  # Lowered threshold slightly
            trial.set_user_attr("failure_reason", "low_alignment_rate")
            trial.set_user_attr("alignment_rate", alignment_rate)
            return 0.0
        
        # Early pruning: check if SAM file has any valid alignments
        if quality_metrics.get("total_alignments", 0) == 0:
            trial.set_user_attr("failure_reason", "no_alignments")
            return 0.0
        
        # Generate bit vectors and analyze
        bit_vector_metrics = generate_bit_vectors_and_analyze(
            output_sam,
            ref_seqs,
            paired=fastq2 is not None,
            qscore_cutoff=25,
            mapq_cutoff=mapq_cutoff,
        )
        
        # Early pruning: no accepted bit vectors
        if bit_vector_metrics.get("accepted_bit_vectors", 0) == 0:
            trial.set_user_attr("failure_reason", "no_accepted_bit_vectors")
            return 0.0
        
        # Calculate quality score
        snr = bit_vector_metrics.get("signal_to_noise", 0.0)
        avg_mapq = quality_metrics.get("avg_mapq", 0.0)
        bv_acceptance_rate = (
            bit_vector_metrics.get("accepted_bit_vectors", 0)
            / bit_vector_metrics.get("total_bit_vectors", 1)
            if bit_vector_metrics.get("total_bit_vectors", 0) > 0
            else 0.0
        )
        
        # Normalize metrics
        normalized_snr = min(snr / 10.0, 1.0) if snr > 0 else 0.0
        
        # Quality score (same weights as main script)
        quality_score = (
            0.40 * normalized_snr  # Maximize signal-to-noise (primary goal)
            + 0.30 * alignment_rate  # Maintain alignment rate
            + 0.20 * (avg_mapq / 60.0)  # Normalize MAPQ to 0-1
            + 0.10 * bv_acceptance_rate  # Bit vector acceptance
        )
        
        # Store intermediate values for analysis
        trial.set_user_attr("signal_to_noise", snr)
        trial.set_user_attr("alignment_rate", alignment_rate)
        trial.set_user_attr("avg_mapq", avg_mapq)
        trial.set_user_attr("bv_accepted", bit_vector_metrics.get("accepted_bit_vectors", 0))
        
        # Store per-construct information
        if bit_vector_metrics.get("constructs"):
            for construct_name, construct_stats in bit_vector_metrics["constructs"].items():
                trial.set_user_attr(f"construct_{construct_name}_aligned", construct_stats["aligned_reads"])
                trial.set_user_attr(f"construct_{construct_name}_snr", construct_stats["signal_to_noise"])
        if optimize_threads:
            trial.set_user_attr("threads_used", trial_threads)
        
        return quality_score
        
    except Exception as e:
        log.warning(f"Trial {trial.number} failed: {e}")
        trial.set_user_attr("failure_reason", f"exception: {str(e)[:50]}")
        return 0.0


def main():
    """Main optimization function using Optuna."""
    parser = argparse.ArgumentParser(
        description="Optimize Bowtie2 parameters using Optuna (Bayesian optimization)"
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
        default=Path("bowtie2_optimization_optuna"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Expected read length (default: 150)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for alignment (or max if optimizing)",
    )
    parser.add_argument(
        "--optimize-threads",
        action="store_true",
        help="Optimize number of threads (default: use fixed threads)",
    )
    parser.add_argument(
        "--mapq-cutoff",
        type=int,
        default=20,
        help="MAPQ cutoff for high-quality alignments (default: 20)",
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=100,
        help="Number of Optuna trials to run (default: 100)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Timeout in seconds (optional)",
    )
    parser.add_argument(
        "--study-name",
        type=str,
        default="bowtie2_optimization",
        help="Name for Optuna study (default: bowtie2_optimization)",
    )
    parser.add_argument(
        "--storage",
        type=str,
        default=None,
        help="Optuna storage URL (e.g., sqlite:///study.db) for resuming studies",
    )
    
    args = parser.parse_args()
    
    # Setup output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_dir = args.output_dir / "index"
    results_dir = args.output_dir / "results"
    results_dir.mkdir(exist_ok=True)
    
    print("=" * 80)
    print("Bowtie2 Parameter Optimization with Optuna")
    print("=" * 80)
    print(f"FASTA: {args.fasta}")
    print(f"FASTQ1: {args.fastq1}")
    if args.fastq2:
        print(f"FASTQ2: {args.fastq2}")
    print(f"Read length: {args.read_length}bp")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of trials: {args.n_trials}")
    print(f"Threads: {args.threads} ({'optimizing' if args.optimize_threads else 'fixed'})")
    print()
    sys.stdout.flush()
    
    # Read reference sequences
    ref_seqs = fasta_to_dict(args.fasta)
    
    # Display constructs found
    print(f"\nReference sequences (constructs) found: {len(ref_seqs)}")
    for i, (name, seq) in enumerate(ref_seqs.items(), 1):
        print(f"  {i}. {name} ({len(seq)} bp)")
    print()
    sys.stdout.flush()
    
    # Build index
    print("Building Bowtie2 index...")
    sys.stdout.flush()
    index = build_bowtie2_index(args.fasta, index_dir)
    
    # Create Optuna study
    study = optuna.create_study(
        direction="maximize",
        study_name=args.study_name,
        storage=args.storage,
        load_if_exists=True,
    )
    
    print("\n" + "=" * 80)
    print("Starting Optuna Optimization")
    print("=" * 80)
    print(f"Study: {args.study_name}")
    print(f"Trials: {args.n_trials}")
    print()
    sys.stdout.flush()
    
    # Run optimization
    study.optimize(
        lambda trial: objective(
            trial,
            index,
            args.fastq1,
            args.fastq2,
            ref_seqs,
            results_dir,
            args.mapq_cutoff,
            args.threads,
            args.optimize_threads,
        ),
        n_trials=args.n_trials,
        timeout=args.timeout,
        show_progress_bar=True,
    )
    
    # Results
    print("\n" + "=" * 80)
    print("Optimization Complete")
    print("=" * 80)
    print(f"Best trial: {study.best_trial.number}")
    print(f"Best value (quality score): {study.best_value:.4f}")
    print(f"Best parameters:")
    for key, value in study.best_params.items():
        print(f"  {key}: {value}")
    print()
    
    # Best trial attributes
    best_snr = study.best_trial.user_attrs.get("signal_to_noise", 0.0)
    best_align_rate = study.best_trial.user_attrs.get("alignment_rate", 0.0)
    best_mapq = study.best_trial.user_attrs.get("avg_mapq", 0.0)
    best_bv = study.best_trial.user_attrs.get("bv_accepted", 0)
    
    print(f"Best trial metrics:")
    print(f"  Signal-to-Noise (AC/GU): {best_snr:.2f}")
    print(f"  Alignment Rate: {best_align_rate:.2%}")
    print(f"  Avg MAPQ: {best_mapq:.1f}")
    print(f"  Bit Vectors Accepted: {best_bv}")
    print()
    
    # Display per-construct statistics
    print("Per-construct statistics:")
    construct_found = False
    for key, value in study.best_trial.user_attrs.items():
        if key.startswith("construct_") and key.endswith("_aligned"):
            construct_name = key.replace("construct_", "").replace("_aligned", "")
            aligned = value
            snr_key = f"construct_{construct_name}_snr"
            snr = study.best_trial.user_attrs.get(snr_key, 0.0)
            print(f"  {construct_name}:")
            print(f"    Aligned reads: {aligned:,}")
            print(f"    Signal-to-Noise (AC/GU): {snr:.2f}")
            construct_found = True
    if not construct_found:
        print("  (No construct-specific data available)")
    print()
    
    # Display per-construct statistics
    print("Per-construct statistics:")
    construct_found = False
    for key, value in study.best_trial.user_attrs.items():
        if key.startswith("construct_") and key.endswith("_aligned"):
            construct_name = key.replace("construct_", "").replace("_aligned", "")
            aligned = value
            snr_key = f"construct_{construct_name}_snr"
            snr = study.best_trial.user_attrs.get(snr_key, 0.0)
            print(f"  {construct_name}:")
            print(f"    Aligned reads: {aligned}")
            print(f"    Signal-to-Noise (AC/GU): {snr:.2f}")
            construct_found = True
    if not construct_found:
        print("  (No construct-specific data available)")
    print()
    
    # Convert best parameters to Bowtie2 arguments
    best_params = {
        "local": True,
        "no_unal": True,
        "no_discordant": True,
        "no_mixed": True,
        "seed_length": study.best_params["seed_length"],
        "maxins": study.best_params["maxins"],
    }
    
    # Add optional parameters if present
    if study.best_params.get("seed_mismatches", 0) > 0:
        best_params["seed_mismatches"] = study.best_params["seed_mismatches"]
    if study.best_params.get("score_min"):
        best_params["score_min"] = study.best_params["score_min"]
    if study.best_params.get("mismatch_penalty"):
        best_params["mismatch_penalty"] = study.best_params["mismatch_penalty"]
    if study.best_params.get("gap_penalty_read"):
        best_params["gap_penalty_read"] = study.best_params["gap_penalty_read"]
    if study.best_params.get("gap_penalty_ref"):
        best_params["gap_penalty_ref"] = study.best_params["gap_penalty_ref"]
    if study.best_params.get("sensitivity_mode"):
        best_params[study.best_params["sensitivity_mode"]] = True
    if study.best_params.get("seed_interval"):
        best_params["seed_interval"] = study.best_params["seed_interval"]
    if "np_penalty" in study.best_params and study.best_params["np_penalty"] is not None:
        best_params["np_penalty"] = study.best_params["np_penalty"]
    if study.best_params.get("n_ceil"):
        best_params["n_ceil"] = study.best_params["n_ceil"]
    if "gbar" in study.best_params and study.best_params["gbar"] is not None:
        best_params["gbar"] = study.best_params["gbar"]
    if "match_bonus" in study.best_params and study.best_params["match_bonus"] is not None:
        best_params["match_bonus"] = study.best_params["match_bonus"]
    if "extension_effort" in study.best_params and study.best_params["extension_effort"] is not None:
        best_params["extension_effort"] = study.best_params["extension_effort"]
    if "repetitive_effort" in study.best_params and study.best_params["repetitive_effort"] is not None:
        best_params["repetitive_effort"] = study.best_params["repetitive_effort"]
    if "minins" in study.best_params and study.best_params["minins"] is not None:
        best_params["minins"] = study.best_params["minins"]
    
    best_bt2_args = params_to_bowtie2_args(best_params)
    print(f"Best Bowtie2 arguments:")
    print(f"  {' '.join(best_bt2_args)}")
    print()
    
    # Save results
    results_file = args.output_dir / "optuna_study.json"
    with open(results_file, "w") as f:
        json.dump({
            "best_trial": {
                "number": study.best_trial.number,
                "value": study.best_value,
                "params": study.best_params,
                "user_attrs": study.best_trial.user_attrs,
            },
            "best_bowtie2_args": " ".join(best_bt2_args),
            "n_trials": len(study.trials),
        }, f, indent=2)
    
    # Save study for later analysis
    study_file = args.output_dir / "optuna_study.pkl"
    import pickle
    with open(study_file, "wb") as f:
        pickle.dump(study, f)
    
    # Create summary DataFrame
    df_data = []
    for trial in study.trials:
        if trial.state == optuna.trial.TrialState.COMPLETE:
            df_data.append({
                "trial": trial.number,
                "quality_score": trial.value,
                "signal_to_noise": trial.user_attrs.get("signal_to_noise", 0.0),
                "alignment_rate": trial.user_attrs.get("alignment_rate", 0.0),
                "avg_mapq": trial.user_attrs.get("avg_mapq", 0.0),
                "bv_accepted": trial.user_attrs.get("bv_accepted", 0),
                **trial.params,
            })
    
    if df_data:
        df = pd.DataFrame(df_data)
        summary_file = args.output_dir / "optuna_summary.csv"
        df.to_csv(summary_file, index=False)
        print(f"Results saved to:")
        print(f"  {results_file}")
        print(f"  {summary_file}")
        print(f"  {study_file}")
        print()
    
    # Generate visualizations if possible
    try:
        print("Generating visualizations...")
        vis_dir = args.output_dir / "visualizations"
        vis_dir.mkdir(exist_ok=True)
        
        # Optimization history
        fig = plot_optimization_history(study)
        fig.write_html(str(vis_dir / "optimization_history.html"))
        
        # Parameter importances
        try:
            fig = plot_param_importances(study)
            fig.write_html(str(vis_dir / "param_importances.html"))
        except Exception:
            pass  # May fail if not enough trials
        
        # Parallel coordinate plot
        try:
            fig = plot_parallel_coordinate(study)
            fig.write_html(str(vis_dir / "parallel_coordinate.html"))
        except Exception:
            pass
        
        print(f"Visualizations saved to: {vis_dir}")
        print()
    except Exception as e:
        log.warning(f"Could not generate visualizations: {e}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

