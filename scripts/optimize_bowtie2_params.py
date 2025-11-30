#!/usr/bin/env python3
"""Optimize Bowtie2 parameters for maximum signal-to-noise ratio.

This script tests different Bowtie2 parameter combinations to find the optimal
balance between alignment quality (signal-to-noise) and number of aligned reads.
Designed for ~150bp sequences typical in RNA structure probing experiments.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import yaml

# rna-map package should be installed from main repository:
# pip install git+https://github.com/jyesselm/rna_map_nextflow.git#subdirectory=python
# Or if published to PyPI: pip install rna-map

from rna_map.io.fasta import fasta_to_dict
from rna_map.logger import get_logger
from rna_map.analysis.bit_vector_iterator import BitVectorIterator
from rna_map.analysis.mutation_histogram import MutationHistogram

log = get_logger("BOWTIE2_OPTIMIZER")


def parse_bowtie2_stats(stderr_file: Path) -> Dict:
    """Parse Bowtie2 alignment statistics from stderr output.
    
    Args:
        stderr_file: Path to Bowtie2 stderr output
        
    Returns:
        Dictionary with alignment statistics including multi-mapping info
    """
    stats = {
        "total_reads": 0,
        "aligned_0_times": 0,
        "aligned_exactly_1_time": 0,
        "aligned_more_than_1_time": 0,
        "overall_alignment_rate": 0.0,
        "concordant_0_times": 0,
        "concordant_exactly_1_time": 0,
        "concordant_more_than_1_time": 0,  # Reads aligning to multiple references
        "is_paired_end": False,
    }
    
    if not stderr_file.exists():
        return stats
    
    with open(stderr_file) as f:
        content = f.read()
        lines = content.split('\n')
        
        for line in lines:
            line = line.strip()
            if "reads; of these:" in line:
                stats["total_reads"] = int(line.split()[0])
            elif "were paired" in line:
                stats["is_paired_end"] = True
            elif "aligned 0 times" in line and "concordantly" not in line:
                stats["aligned_0_times"] = int(line.split()[0])
            elif "aligned exactly 1 time" in line and "concordantly" not in line:
                stats["aligned_exactly_1_time"] = int(line.split()[0])
            elif "aligned >1 times" in line and "concordantly" not in line:
                stats["aligned_more_than_1_time"] = int(line.split()[0])
            elif "aligned concordantly 0 times" in line:
                stats["concordant_0_times"] = int(line.split()[0])
            elif "aligned concordantly exactly 1 time" in line:
                stats["concordant_exactly_1_time"] = int(line.split()[0])
            elif "aligned concordantly >1 times" in line:
                stats["concordant_more_than_1_time"] = int(line.split()[0])
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
    
    log.info(f"Running: {' '.join(cmd)}")
    
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
        total_snr = 0.0
        total_weight = 0
        for ref_name, mh in mut_histos.items():
            if mh.num_aligned > 0:
                snr = mh.get_signal_to_noise()
                total_snr += snr * mh.num_aligned
                total_weight += mh.num_aligned
        
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


def calculate_signal_to_noise(
    alignment_stats: Dict,
    quality_metrics: Dict,
    bit_vector_metrics: Dict | None = None,
    mapq_cutoff: int = 20,
) -> Dict:
    """Calculate signal-to-noise ratio and quality score.
    
    Args:
        alignment_stats: Alignment statistics from Bowtie2
        quality_metrics: Quality metrics from SAM file
        bit_vector_metrics: Bit vector metrics including signal-to-noise
        mapq_cutoff: MAPQ cutoff for high-quality alignments
        
    Returns:
        Dictionary with signal-to-noise metrics
    """
    total_reads = alignment_stats.get("total_reads", 0)
    aligned_reads = alignment_stats.get("aligned_exactly_1_time", 0)
    high_quality = quality_metrics.get("high_quality_alignments", 0)
    avg_mapq = quality_metrics.get("avg_mapq", 0.0)
    
    # Alignment rate
    alignment_rate = (
        alignment_stats.get("overall_alignment_rate", 0.0) / 100.0
    )
    
    # Get signal-to-noise from mutation histogram (AC/GU ratio)
    # This is the actual signal-to-noise metric used in RNA-MAP
    if bit_vector_metrics and "signal_to_noise" in bit_vector_metrics:
        snr = bit_vector_metrics["signal_to_noise"]
    else:
        # Fallback: calculate from alignment stats if no bit vector metrics
        signal = high_quality
        noise = total_reads - high_quality
        snr = signal / noise if noise > 0 else float("inf")
    
    # Include bit vector metrics if available
    bv_metrics = {}
    reads_3plus_mutations = 0
    reads_3plus_rate = 0.0
    
    if bit_vector_metrics:
        bv_accepted = bit_vector_metrics.get("accepted_bit_vectors", 0)
        bv_total = bit_vector_metrics.get("total_bit_vectors", 0)
        bv_acceptance_rate = bv_accepted / bv_total if bv_total > 0 else 0.0
        avg_muts = bit_vector_metrics.get("avg_mutations_per_read", 0.0)
        reads_with_muts = bit_vector_metrics.get("reads_with_mutations", 0)
        reads_total = bv_accepted if bv_accepted > 0 else 1
        mutation_rate = reads_with_muts / reads_total
        
        # Count reads with 3+ mutations (penalty factor)
        mut_dist = bit_vector_metrics.get("mutation_distribution", {})
        reads_3plus_mutations = sum(
            count for mut_count, count in mut_dist.items() 
            if int(mut_count) >= 3
        )
        reads_3plus_rate = reads_3plus_mutations / reads_total if reads_total > 0 else 0.0
        
        bv_metrics = {
            "bit_vector_acceptance_rate": bv_acceptance_rate,
            "avg_mutations_per_read": avg_muts,
            "mutation_detection_rate": mutation_rate,
            "accepted_bit_vectors": bv_accepted,
            "reads_3plus_mutations": reads_3plus_mutations,
            "reads_3plus_rate": reads_3plus_rate,
        }
    
    # Quality score: weighted combination
    # Prioritizes signal-to-noise maximization (AC/GU ratio from mutation histogram)
    # Note: 3+ mutations are filtered in bit vector generation, so not penalized here
    if bit_vector_metrics:
        # Normalize SNR (AC/GU ratio typically ranges 0-10, cap at 10 for normalization)
        normalized_snr = min(snr / 10.0, 1.0) if snr > 0 else 0.0
        
        # Quality score weights:
        # - 40% signal-to-noise (maximize - primary goal) - AC/GU ratio
        # - 30% alignment rate (maintain good coverage)
        # - 20% MAPQ quality
        # - 10% bit vector acceptance rate
        quality_score = (
            0.40 * normalized_snr  # Maximize signal-to-noise (primary goal)
            + 0.30 * alignment_rate  # Maintain alignment rate
            + 0.20 * (avg_mapq / 60.0)  # Normalize MAPQ to 0-1
            + 0.10 * bv_metrics["bit_vector_acceptance_rate"]  # Bit vector acceptance
        )
    else:
        # Fallback scoring without bit vector metrics
        normalized_snr = min(snr / 10.0, 1.0) if snr != float("inf") and snr > 0 else 0.0
        quality_score = (
            0.50 * normalized_snr  # Maximize signal-to-noise (primary goal)
            + 0.30 * alignment_rate  # Maintain alignment rate
            + 0.20 * (avg_mapq / 60.0)  # Normalize MAPQ to 0-1
        )
    
    result = {
        "signal_to_noise": snr,  # AC/GU ratio from mutation histogram
        "alignment_rate": alignment_rate,
        "high_quality_rate": high_quality / total_reads if total_reads > 0 else 0.0,
        "quality_score": quality_score,
    }
    result.update(bv_metrics)
    
    return result


def get_baseline_parameters(num_sequences: int = 1) -> Dict:
    """Get baseline parameter set (default/current settings).
    
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
        "max_alignments": 1,  # -k 1: only report best alignment to reduce multi-mapping
    }
    
    if num_sequences > 1:
        # Stricter parameters for multiple sequences to ensure best match
        base_params.update({
            "seed_length": 15,  # Longer seed for better specificity
            "mismatch_penalty": "6,2",  # Stricter mismatch penalty
            "score_min": "L,10,0.2",  # Higher score threshold to prefer better matches
            "gap_penalty_read": "8,4",  # Stricter gap penalties
            "gap_penalty_ref": "5,3",
        })
    else:
        # More permissive for single sequence (allows mutations)
        base_params.update({
            "seed_length": 12,
        })
    
    return base_params


def generate_parameter_combinations(
    read_length: int = 150, 
    num_sequences: int = 1
) -> List[Dict]:
    """Generate parameter combinations to test.
    
    Args:
        read_length: Expected read length (default: 150bp)
        num_sequences: Number of reference sequences (1 = single, >1 = multiple)
        
    Returns:
        List of parameter dictionaries
    """
    combinations = []
    
    # Base parameters
    base_params = {
        "local": True,
        "no_unal": True,
        "no_discordant": True,
        "no_mixed": True,
        "max_alignments": 1,  # -k 1: only report best alignment to reduce multi-mapping
    }
    
    # Different parameter ranges for single vs multiple sequences
    if num_sequences > 1:
        # Multiple sequences: stricter parameters for better discrimination
        seed_lengths = [15, 18, 20]  # Longer seeds for specificity
        seed_mismatches = [0]  # No seed mismatches for strict matching
        maxins_values = [200, 300, 500] if read_length <= 150 else [500, 1000]
        score_mins = [
            "L,10,0.2",  # Higher threshold to prefer better matches
            "L,10,0.3",
            "G,20,15",
        ]
        mismatch_penalties = [
            "6,2",  # Stricter penalties
            "6,3",
        ]
        gap_penalty_reads = ["8,4", "10,5"]  # Stricter gap penalties
        gap_penalty_refs = ["5,3", "6,4"]
        seed_intervals = [
            "S,1,0.5",  # More seeds for better discrimination
            "S,1,0.75",
            "S,1,1.0",
        ]
        extension_efforts = [20, 25]  # High effort to find best match
        repetitive_efforts = [3, 4]  # High effort for similar sequences
        sensitivity_modes = [
            "sensitive-local",  # More sensitive for discrimination
            "very-sensitive-local",
        ]
    else:
        # Single sequence: more permissive parameters to allow mutations
        seed_lengths = [10, 12, 15]  # Shorter seeds allow mutations
        seed_mismatches = [0, 1]  # Allow seed mismatches for mutations
        maxins_values = [200, 300, 500, 1000] if read_length <= 150 else [500, 1000, 2000]
        score_mins = [
            None,
            "L,0,-0.6",  # Lower threshold to allow mutations
            "L,0,-0.4",
            "G,20,15",
        ]
        mismatch_penalties = [
            None,
            "6,2",
            "4,2",  # Lower penalty for mutations
        ]
        gap_penalty_reads = [None, "8,4", "5,3"]
        gap_penalty_refs = [None, "5,3", "4,2"]
        seed_intervals = [
            None,
            "S,1,1.15",  # Default or fewer seeds for speed
            "S,1,1.5",
            "S,1,2.0",
        ]
        extension_efforts = [10, 15, 20]  # Moderate effort
        repetitive_efforts = [2, 3]  # Moderate effort
        sensitivity_modes = [
            None,
            "very-fast-local",
            "fast-local",
            "sensitive-local",
        ]
    
    # Generate combinations with new parameters
    # Limit combinations to avoid explosion - sample strategically
    for seed_len in seed_lengths:
        for seed_mismatch in seed_mismatches:
            for maxins in maxins_values[:2]:  # Limit maxins combinations
                for score_min in score_mins[:2]:  # Limit score_min combinations
                    for mp in mismatch_penalties:
                        for rdg in gap_penalty_reads[:2]:  # Limit gap penalty combinations
                            for rfg in gap_penalty_refs[:2]:
                                for seed_int in seed_intervals[:2]:  # Limit seed interval combinations
                                    for ext_effort in extension_efforts[:2]:  # Limit extension effort
                                        for rep_effort in repetitive_efforts[:2]:  # Limit repetitive effort
                                            for sens_mode in sensitivity_modes[:2]:  # Limit sensitivity modes
                                                params = base_params.copy()
                                                params["seed_length"] = seed_len
                                                params["seed_mismatches"] = seed_mismatch
                                                params["maxins"] = maxins
                                                if score_min:
                                                    params["score_min"] = score_min
                                                if mp:
                                                    params["mismatch_penalty"] = mp
                                                if rdg:
                                                    params["gap_penalty_read"] = rdg
                                                if rfg:
                                                    params["gap_penalty_ref"] = rfg
                                                if seed_int:
                                                    params["seed_interval"] = seed_int
                                                params["extension_effort"] = ext_effort
                                                params["repetitive_effort"] = rep_effort
                                                if sens_mode:
                                                    params[sens_mode] = True
                                                
                                                combinations.append(params)
    
    # Add targeted combinations based on mode
    if num_sequences > 1:
        # Multiple sequences: strict mode for best discrimination
        targeted = [
            {
                **base_params,
                "seed_length": 18,
                "seed_mismatches": 0,
                "maxins": 300,
                "mismatch_penalty": "6,2",
                "gap_penalty_read": "8,4",
                "gap_penalty_ref": "5,3",
                "score_min": "L,10,0.2",
                "seed_interval": "S,1,0.5",
                "extension_effort": 25,
                "repetitive_effort": 4,
                "sensitive-local": True,
            },
            {
                **base_params,
                "seed_length": 20,
                "seed_mismatches": 0,
                "maxins": 300,
                "mismatch_penalty": "6,3",
                "gap_penalty_read": "10,5",
                "gap_penalty_ref": "6,4",
                "score_min": "L,10,0.3",
                "seed_interval": "S,1,0.75",
                "extension_effort": 25,
                "repetitive_effort": 4,
                "very-sensitive-local": True,
            },
            {
                **base_params,
                "seed_length": 15,
                "seed_mismatches": 0,
                "maxins": 200,
                "mismatch_penalty": "6,2",
                "gap_penalty_read": "8,4",
                "gap_penalty_ref": "5,3",
                "score_min": "L,10,0.2",
                "seed_interval": "S,1,0.75",
                "extension_effort": 20,
                "repetitive_effort": 3,
                "sensitive-local": True,
            },
        ]
    else:
        # Single sequence: permissive mode for mutation detection
        targeted = [
            {
                **base_params,
                "seed_length": 10,
                "seed_mismatches": 0,
                "maxins": 300,
                "mismatch_penalty": "6,2",
                "score_min": "L,0,-0.6",
                "seed_interval": "S,1,1.5",
                "extension_effort": 15,
                "repetitive_effort": 2,
                "fast-local": True,
            },
            {
                **base_params,
                "seed_length": 10,
                "seed_mismatches": 1,
                "maxins": 300,
                "mismatch_penalty": "4,2",
                "score_min": "L,0,-0.6",
                "seed_interval": "S,1,2.0",
                "extension_effort": 15,
                "repetitive_effort": 2,
                "fast-local": True,
            },
            {
                **base_params,
                "seed_length": 12,
                "seed_mismatches": 0,
                "maxins": 300,
                "mismatch_penalty": "6,2",
                "score_min": "L,0,-0.6",
                "extension_effort": 15,
                "repetitive_effort": 2,
            },
            {
                **base_params,
                "seed_length": 12,
                "seed_mismatches": 1,
                "maxins": 300,
                "mismatch_penalty": "6,2",
                "score_min": "L,0,-0.6",
                "extension_effort": 15,
                "repetitive_effort": 2,
            },
        ]
    
    combinations.extend(targeted)
    
    return combinations


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
        "max_alignments": "-k",  # Report up to N alignments per read
        "score_min": "--score-min",
        "mismatch_penalty": "--mp",
        "gap_penalty_read": "--rdg",
        "gap_penalty_ref": "--rfg",
        "seed_interval": "-i",  # Seed interval for better discrimination
        "extension_effort": "-D",  # Extension effort
        "repetitive_effort": "-R",  # Repetitive seed effort
        "best": "--best",  # Report alignments in best-to-worst order
        "strata": "--strata",  # Only report alignments in best alignment stratum
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


def main():
    """Main optimization function."""
    parser = argparse.ArgumentParser(
        description="Optimize Bowtie2 parameters for signal-to-noise ratio"
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
        default=Path("bowtie2_optimization"),
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
        help="Number of threads for alignment",
    )
    parser.add_argument(
        "--mapq-cutoff",
        type=int,
        default=20,
        help="MAPQ cutoff for high-quality alignments (default: 20)",
    )
    parser.add_argument(
        "--max-combinations",
        type=int,
        default=50,
        help="Maximum number of parameter combinations to test (default: 50)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick test with fewer combinations",
    )
    
    args = parser.parse_args()
    
    # Setup output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_dir = args.output_dir / "index"
    results_dir = args.output_dir / "results"
    results_dir.mkdir(exist_ok=True)
    
    print("=" * 80)
    print("Bowtie2 Parameter Optimization")
    print("=" * 80)
    print(f"FASTA: {args.fasta}")
    print(f"FASTQ1: {args.fastq1}")
    if args.fastq2:
        print(f"FASTQ2: {args.fastq2}")
    print(f"Read length: {args.read_length}bp")
    print(f"Output directory: {args.output_dir}")
    sys.stdout.flush()
    
    log.info("=" * 80)
    log.info("Bowtie2 Parameter Optimization")
    log.info("=" * 80)
    log.info(f"FASTA: {args.fasta}")
    log.info(f"FASTQ1: {args.fastq1}")
    if args.fastq2:
        log.info(f"FASTQ2: {args.fastq2}")
    log.info(f"Read length: {args.read_length}bp")
    log.info(f"Output directory: {args.output_dir}")
    
    # Read reference sequences
    ref_seqs = fasta_to_dict(args.fasta)
    
    # Detect number of sequences to determine mode
    num_sequences = len(ref_seqs)
    print(f"\nDetected {num_sequences} reference sequence(s)")
    if num_sequences > 1:
        print("Using MULTIPLE-SEQUENCE mode with stricter parameters for better discrimination")
        log.info(f"Using MULTIPLE-SEQUENCE mode ({num_sequences} sequences)")
    else:
        print("Using SINGLE-SEQUENCE mode with permissive parameters for mutation detection")
        log.info("Using SINGLE-SEQUENCE mode")
    sys.stdout.flush()
    
    # Build index
    print("\n" + "=" * 80)
    print("Building Bowtie2 index")
    print("=" * 80)
    sys.stdout.flush()
    log.info("\n" + "=" * 80)
    log.info("Building Bowtie2 index")
    log.info("=" * 80)
    index = build_bowtie2_index(args.fasta, index_dir)
    
    # Generate parameter combinations
    print("\n" + "=" * 80)
    print("Generating parameter combinations")
    print("=" * 80)
    sys.stdout.flush()
    log.info("\n" + "=" * 80)
    log.info("Generating parameter combinations")
    log.info("=" * 80)
    all_combinations = generate_parameter_combinations(
        args.read_length, 
        num_sequences=num_sequences
    )
    
    if args.quick:
        # Select diverse subset for quick test
        step = max(1, len(all_combinations) // 20)
        combinations = all_combinations[::step][:20]
    else:
        combinations = all_combinations[: args.max_combinations]
    
    print(f"Testing {len(combinations)} parameter combinations")
    sys.stdout.flush()
    log.info(f"Testing {len(combinations)} parameter combinations")
    
    # Test baseline first (num_sequences already detected above)
    baseline_params = get_baseline_parameters(num_sequences=num_sequences)
    baseline_result = None
    
    print("\n" + "=" * 80)
    print("Testing Baseline Configuration")
    print("=" * 80)
    sys.stdout.flush()
    log.info("\n" + "=" * 80)
    log.info("Testing Baseline Configuration")
    log.info("=" * 80)
    
    print("\n[Baseline] Testing baseline configuration:")
    param_str = ", ".join(f"{k}={v}" for k, v in baseline_params.items() if v)
    print(f"  {param_str}")
    sys.stdout.flush()
    
    bt2_args = params_to_bowtie2_args(baseline_params)
    baseline_dir = results_dir / "baseline"
    baseline_dir.mkdir(exist_ok=True)
    baseline_sam = baseline_dir / "aligned.sam"
    
    try:
        stderr_file, elapsed_time = run_bowtie2_alignment(
            index,
            args.fastq1,
            args.fastq2,
            baseline_sam,
            bt2_args,
            threads=args.threads,
        )
        
        alignment_stats = parse_bowtie2_stats(stderr_file)
        quality_metrics = parse_sam_quality(baseline_sam, args.mapq_cutoff)
        
        print(f"  Generating bit vectors...")
        sys.stdout.flush()
        bit_vector_metrics = generate_bit_vectors_and_analyze(
            baseline_sam,
            ref_seqs,
            paired=args.fastq2 is not None,
            qscore_cutoff=25,
            mapq_cutoff=args.mapq_cutoff,
        )
        
        snr_metrics = calculate_signal_to_noise(
            alignment_stats, 
            quality_metrics, 
            bit_vector_metrics,
            args.mapq_cutoff
        )
        
        baseline_result = {
            "combination_id": 0,
            "parameters": baseline_params,
            "bowtie2_args": " ".join(bt2_args),
            "alignment_stats": alignment_stats,
            "quality_metrics": quality_metrics,
            "bit_vector_metrics": bit_vector_metrics,
            "signal_to_noise": snr_metrics,
            "elapsed_time": elapsed_time,
        }
        
        bv_info = ""
        if "accepted_bit_vectors" in snr_metrics:
            reads_3plus = snr_metrics.get("reads_3plus_mutations", 0)
            reads_3plus_pct = snr_metrics.get("reads_3plus_rate", 0.0) * 100
            bv_info = (
                f", BV accepted: {snr_metrics['accepted_bit_vectors']}, "
                f"Avg muts: {snr_metrics.get('avg_mutations_per_read', 0):.2f}, "
                f"3+ muts: {reads_3plus} ({reads_3plus_pct:.1f}%)"
            )
        
        # Report multi-mapping statistics
        multi_map_info = ""
        if alignment_stats.get("is_paired_end"):
            multi_map_count = alignment_stats.get("concordant_more_than_1_time", 0)
            if multi_map_count > 0:
                multi_map_pct = (multi_map_count / alignment_stats["total_reads"]) * 100
                multi_map_info = f", Multi-mapping reads: {multi_map_count} ({multi_map_pct:.2f}%)"
            else:
                multi_map_info = ", Multi-mapping reads: 0 (0.00%)"
        else:
            multi_map_count = alignment_stats.get("aligned_more_than_1_time", 0)
            if multi_map_count > 0:
                multi_map_pct = (multi_map_count / alignment_stats["total_reads"]) * 100
                multi_map_info = f", Multi-mapping reads: {multi_map_count} ({multi_map_pct:.2f}%)"
            else:
                multi_map_info = ", Multi-mapping reads: 0 (0.00%)"
        
        print(
            f"  Alignment rate: {snr_metrics['alignment_rate']:.2%}, "
            f"High-quality: {snr_metrics['high_quality_rate']:.2%}, "
            f"SNR: {snr_metrics['signal_to_noise']:.2f}, "
            f"Quality score: {snr_metrics['quality_score']:.3f}"
            f"{bv_info}"
            f"{multi_map_info}"
        )
        sys.stdout.flush()
        
    except Exception as e:
        log.error(f"  Error testing baseline: {e}")
        print(f"  Error: {e}")
    
    # Test each combination
    results = []
    
    print("\n" + "=" * 80)
    print("Testing Parameter Combinations")
    print("=" * 80)
    sys.stdout.flush()
    log.info("\n" + "=" * 80)
    log.info("Testing Parameter Combinations")
    log.info("=" * 80)
    
    for i, params in enumerate(combinations, 1):
        print(f"\n[{i}/{len(combinations)}] Testing combination:")
        param_str = ", ".join(f"{k}={v}" for k, v in params.items() if v)
        print(f"  {param_str}")
        sys.stdout.flush()
        log.info(f"\n[{i}/{len(combinations)}] Testing combination:")
        log.info(f"  {param_str}")
        
        # Convert to Bowtie2 arguments
        bt2_args = params_to_bowtie2_args(params)
        
        # Create output directory for this combination
        combo_dir = results_dir / f"combo_{i:03d}"
        combo_dir.mkdir(exist_ok=True)
        output_sam = combo_dir / "aligned.sam"
        
        try:
            # Run alignment
            stderr_file, elapsed_time = run_bowtie2_alignment(
                index,
                args.fastq1,
                args.fastq2,
                output_sam,
                bt2_args,
                threads=args.threads,
            )
            
            # Parse statistics
            alignment_stats = parse_bowtie2_stats(stderr_file)
            quality_metrics = parse_sam_quality(output_sam, args.mapq_cutoff)
            
            # Generate bit vectors and analyze
            print(f"  Generating bit vectors...")
            sys.stdout.flush()
            bit_vector_metrics = generate_bit_vectors_and_analyze(
                output_sam,
                ref_seqs,
                paired=args.fastq2 is not None,
                qscore_cutoff=25,  # Default quality score cutoff
                mapq_cutoff=args.mapq_cutoff,
            )
            
            snr_metrics = calculate_signal_to_noise(
                alignment_stats, 
                quality_metrics, 
                bit_vector_metrics,
                args.mapq_cutoff
            )
            
            # Store results
            result = {
                "combination_id": i,
                "parameters": params,
                "bowtie2_args": " ".join(bt2_args),
                "alignment_stats": alignment_stats,
                "quality_metrics": quality_metrics,
                "bit_vector_metrics": bit_vector_metrics,
                "signal_to_noise": snr_metrics,
                "elapsed_time": elapsed_time,
            }
            results.append(result)
            
            bv_info = ""
            if "accepted_bit_vectors" in snr_metrics:
                reads_3plus = snr_metrics.get("reads_3plus_mutations", 0)
                reads_3plus_pct = snr_metrics.get("reads_3plus_rate", 0.0) * 100
                bv_info = (
                    f", BV accepted: {snr_metrics['accepted_bit_vectors']}, "
                    f"Avg muts: {snr_metrics.get('avg_mutations_per_read', 0):.2f}, "
                    f"3+ muts: {reads_3plus} ({reads_3plus_pct:.1f}%)"
                )
            
            # Report multi-mapping statistics
            multi_map_info = ""
            if alignment_stats.get("is_paired_end"):
                multi_map_count = alignment_stats.get("concordant_more_than_1_time", 0)
                if multi_map_count > 0:
                    multi_map_pct = (multi_map_count / alignment_stats["total_reads"]) * 100
                    multi_map_info = f", Multi-mapping: {multi_map_count} ({multi_map_pct:.2f}%)"
                else:
                    multi_map_info = ", Multi-mapping: 0 (0.00%)"
            else:
                multi_map_count = alignment_stats.get("aligned_more_than_1_time", 0)
                if multi_map_count > 0:
                    multi_map_pct = (multi_map_count / alignment_stats["total_reads"]) * 100
                    multi_map_info = f", Multi-mapping: {multi_map_count} ({multi_map_pct:.2f}%)"
                else:
                    multi_map_info = ", Multi-mapping: 0 (0.00%)"
            
            print(
                f"  Alignment rate: {snr_metrics['alignment_rate']:.2%}, "
                f"High-quality: {snr_metrics['high_quality_rate']:.2%}, "
                f"SNR: {snr_metrics['signal_to_noise']:.2f}, "
                f"Quality score: {snr_metrics['quality_score']:.3f}"
                f"{bv_info}"
                f"{multi_map_info}"
            )
            sys.stdout.flush()
            log.info(
                f"  Alignment rate: {snr_metrics['alignment_rate']:.2%}, "
                f"High-quality: {snr_metrics['high_quality_rate']:.2%}, "
                f"SNR: {snr_metrics['signal_to_noise']:.2f}, "
                f"Quality score: {snr_metrics['quality_score']:.3f}"
            )
            
        except Exception as e:
            log.error(f"  Error testing combination: {e}")
            continue
    
    # Analyze results
    print("\n" + "=" * 80)
    print("Results Analysis")
    print("=" * 80)
    sys.stdout.flush()
    log.info("\n" + "=" * 80)
    log.info("Results Analysis")
    log.info("=" * 80)
    
    if not results:
        log.error("No successful alignments!")
        return 1
    
    # Add baseline to results for comparison
    if baseline_result:
        results.insert(0, baseline_result)
    
    # Create DataFrame for analysis
    df_data = []
    baseline_quality_score = None
    baseline_snr = None
    baseline_alignment_rate = None
    
    for r in results:
        row = {
            "combo_id": r["combination_id"],
            "alignment_rate": r["signal_to_noise"]["alignment_rate"],
            "high_quality_rate": r["signal_to_noise"]["high_quality_rate"],
            "signal_to_noise": r["signal_to_noise"]["signal_to_noise"],
            "quality_score": r["signal_to_noise"]["quality_score"],
            "avg_mapq": r["quality_metrics"]["avg_mapq"],
            "total_alignments": r["quality_metrics"]["total_alignments"],
            "high_quality_alignments": r["quality_metrics"]["high_quality_alignments"],
            "elapsed_time": r["elapsed_time"],
        }
        
        # Store baseline metrics for comparison
        if r["combination_id"] == 0:
            baseline_quality_score = r["signal_to_noise"]["quality_score"]
            baseline_snr = r["signal_to_noise"]["signal_to_noise"]
            baseline_alignment_rate = r["signal_to_noise"]["alignment_rate"]
            row["is_baseline"] = True
        else:
            row["is_baseline"] = False
        
        # Add bit vector metrics if available
        if "bit_vector_metrics" in r and r["bit_vector_metrics"]:
            bv = r["bit_vector_metrics"]
            row.update({
                "bv_accepted": bv.get("accepted_bit_vectors", 0),
                "bv_total": bv.get("total_bit_vectors", 0),
                "bv_acceptance_rate": r["signal_to_noise"].get("bit_vector_acceptance_rate", 0.0),
                "avg_mutations_per_read": r["signal_to_noise"].get("avg_mutations_per_read", 0.0),
                "mutation_detection_rate": r["signal_to_noise"].get("mutation_detection_rate", 0.0),
                "reads_3plus_mutations": r["signal_to_noise"].get("reads_3plus_mutations", 0),
                "reads_3plus_rate": r["signal_to_noise"].get("reads_3plus_rate", 0.0),
            })
        
        df_data.append(row)
    
    # Add relative improvements compared to baseline
    if baseline_quality_score is not None:
        for row in df_data:
            if not row.get("is_baseline", False):
                # Calculate relative improvements
                quality_improvement = (
                    (row["quality_score"] - baseline_quality_score) / baseline_quality_score * 100
                    if baseline_quality_score > 0 else 0.0
                )
                
                # For SNR, handle infinity
                if baseline_snr == float("inf") and row["signal_to_noise"] == float("inf"):
                    snr_improvement = 0.0
                elif baseline_snr == float("inf"):
                    snr_improvement = -100.0  # Baseline is inf, this is worse
                elif row["signal_to_noise"] == float("inf"):
                    snr_improvement = 100.0  # This is inf, baseline is not
                else:
                    snr_improvement = (
                        (row["signal_to_noise"] - baseline_snr) / baseline_snr * 100
                        if baseline_snr > 0 else 0.0
                    )
                
                alignment_improvement = (
                    (row["alignment_rate"] - baseline_alignment_rate) / baseline_alignment_rate * 100
                    if baseline_alignment_rate > 0 else 0.0
                )
                
                row["quality_score_vs_baseline"] = quality_improvement
                row["snr_vs_baseline"] = snr_improvement
                row["alignment_rate_vs_baseline"] = alignment_improvement
            else:
                row["quality_score_vs_baseline"] = 0.0
                row["snr_vs_baseline"] = 0.0
                row["alignment_rate_vs_baseline"] = 0.0
    
    df = pd.DataFrame(df_data)
    
    # Find best combination (exclude baseline from best selection)
    df_test = df[df["combo_id"] != 0] if 0 in df["combo_id"].values else df
    best_by_quality = df_test.loc[df_test["quality_score"].idxmax()]
    best_by_snr = df_test.loc[df_test["signal_to_noise"].idxmax()]
    best_by_alignment = df_test.loc[df_test["alignment_rate"].idxmax()]
    
    # Show baseline summary
    if baseline_result:
        baseline_row = df[df["combo_id"] == 0].iloc[0] if 0 in df["combo_id"].values else None
        if baseline_row is not None:
            print("\n" + "=" * 80)
            print("Baseline Configuration")
            print("=" * 80)
            print(f"  Bowtie2 arguments: {baseline_result['bowtie2_args']}")
            print(f"  Quality Score: {baseline_row['quality_score']:.3f}")
            print(f"  Alignment Rate: {baseline_row['alignment_rate']:.2%}")
            print(f"  Signal-to-Noise: {baseline_row['signal_to_noise']:.2f}")
            if "bv_accepted" in baseline_row:
                print(f"  Bit Vectors Accepted: {int(baseline_row['bv_accepted'])}")
            print()
    
    print("=" * 80)
    print("Best by Quality Score (vs Baseline)")
    print("=" * 80)
    print(f"  Combination ID: {int(best_by_quality['combo_id'])}")
    print(f"  Quality Score: {best_by_quality['quality_score']:.3f}", end="")
    if "quality_score_vs_baseline" in best_by_quality:
        improvement = best_by_quality["quality_score_vs_baseline"]
        sign = "+" if improvement >= 0 else ""
        print(f" ({sign}{improvement:.1f}% vs baseline)")
    else:
        print()
    print(f"  Alignment Rate: {best_by_quality['alignment_rate']:.2%}", end="")
    if "alignment_rate_vs_baseline" in best_by_quality:
        improvement = best_by_quality["alignment_rate_vs_baseline"]
        sign = "+" if improvement >= 0 else ""
        print(f" ({sign}{improvement:.1f}% vs baseline)")
    else:
        print()
    print(f"  High-Quality Rate: {best_by_quality['high_quality_rate']:.2%}")
    print(f"  Signal-to-Noise: {best_by_quality['signal_to_noise']:.2f}", end="")
    if "snr_vs_baseline" in best_by_quality:
        improvement = best_by_quality["snr_vs_baseline"]
        sign = "+" if improvement >= 0 else ""
        print(f" ({sign}{improvement:.1f}% vs baseline)")
    else:
        print()
    print(f"  Avg MAPQ: {best_by_quality['avg_mapq']:.1f}")
    
    # Show bit vector metrics if available
    if "bv_accepted" in best_by_quality:
        print(f"  Bit Vectors Accepted: {int(best_by_quality['bv_accepted'])}")
        print(f"  Avg Mutations/Read: {best_by_quality.get('avg_mutations_per_read', 0):.2f}")
        print(f"  Mutation Detection Rate: {best_by_quality.get('mutation_detection_rate', 0):.2%}")
        # Note: 3+ mutations are filtered in bit vector generation, so shown for info only
        reads_3plus = best_by_quality.get('reads_3plus_mutations', 0)
        reads_3plus_rate = best_by_quality.get('reads_3plus_rate', 0.0)
        if reads_3plus > 0:
            print(f"  Reads with 3+ Mutations (before filtering): {reads_3plus} ({reads_3plus_rate:.2%})")
    
    # Get best result (handle baseline at index 0)
    best_idx = int(best_by_quality["combo_id"])
    if best_idx == 0:
        best_result = baseline_result
    else:
        # Find result with matching combo_id (accounting for baseline at index 0)
        best_result = next((r for r in results if r["combination_id"] == best_idx), None)
        if best_result is None:
            # Fallback: use index directly if combo_id matches position
            best_result = results[best_idx - 1] if best_idx <= len(results) else results[0]
    
    if best_result:
        print(f"\n  Bowtie2 arguments:")
        print(f"    {best_result['bowtie2_args']}")
    sys.stdout.flush()
    log.info("\nBest by Quality Score (balanced):")
    log.info(f"  Combination ID: {int(best_by_quality['combo_id'])}")
    log.info(f"  Quality Score: {best_by_quality['quality_score']:.3f}")
    log.info(f"  Alignment Rate: {best_by_quality['alignment_rate']:.2%}")
    log.info(f"  High-Quality Rate: {best_by_quality['high_quality_rate']:.2%}")
    log.info(f"  Signal-to-Noise: {best_by_quality['signal_to_noise']:.2f}")
    log.info(f"  Avg MAPQ: {best_by_quality['avg_mapq']:.1f}")
    
    log.info(f"\n  Bowtie2 arguments:")
    log.info(f"    {best_result['bowtie2_args']}")
    
    # Save results
    results_file = args.output_dir / "optimization_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    
    summary_file = args.output_dir / "optimization_summary.csv"
    df.to_csv(summary_file, index=False)
    
    # Save best parameters
    best_params_file = args.output_dir / "best_parameters.json"
    
    # Helper to get result by combo_id
    def get_result_by_id(combo_id: int):
        if combo_id == 0:
            return baseline_result
        return next((r for r in results if r["combination_id"] == combo_id), None)
    
    best_snr_result = get_result_by_id(int(best_by_snr["combo_id"]))
    best_align_result = get_result_by_id(int(best_by_alignment["combo_id"]))
    
    output_data = {
        "baseline": {
            "parameters": baseline_params if baseline_result else None,
            "bowtie2_args": baseline_result["bowtie2_args"] if baseline_result else None,
            "metrics": baseline_result["signal_to_noise"] if baseline_result else None,
        },
        "best_by_quality_score": {
            "combination_id": int(best_by_quality["combo_id"]),
            "parameters": best_result["parameters"] if best_result else None,
            "bowtie2_args": best_result["bowtie2_args"] if best_result else None,
            "metrics": best_result["signal_to_noise"] if best_result else None,
        },
    }
    
    if best_snr_result:
        output_data["best_by_signal_to_noise"] = {
            "combination_id": int(best_by_snr["combo_id"]),
            "parameters": best_snr_result["parameters"],
            "bowtie2_args": best_snr_result["bowtie2_args"],
        }
    
    if best_align_result:
        output_data["best_by_alignment_rate"] = {
            "combination_id": int(best_by_alignment["combo_id"]),
            "parameters": best_align_result["parameters"],
            "bowtie2_args": best_align_result["bowtie2_args"],
        }
    
    with open(best_params_file, "w") as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to:")
    print(f"  {results_file}")
    print(f"  {summary_file}")
    print(f"  {best_params_file}")
    sys.stdout.flush()
    log.info(f"\nResults saved to:")
    log.info(f"  {results_file}")
    log.info(f"  {summary_file}")
    log.info(f"  {best_params_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

