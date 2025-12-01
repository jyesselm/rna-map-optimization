"""Utility functions for RNA-MAP optimization.

These functions provide compatibility with rna-map package interfaces.
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
from collections import defaultdict


def fasta_to_dict(fasta_path: Path) -> Dict[str, str]:
    """Read FASTA file and return dictionary of sequence names to sequences.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                
                # Start new sequence
                current_name = line[1:].split()[0]  # Get name (first word after >)
                current_seq = []
            else:
                # Add to current sequence
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def extract_per_construct_stats(trial) -> Dict[str, Dict[str, Any]]:
    """Extract per-construct statistics from a trial's user attributes.
    
    Args:
        trial: Optuna trial object
        
    Returns:
        Dictionary mapping construct names to their statistics
    """
    constructs = {}
    user_attrs = trial.user_attrs
    
    # Find all construct-related attributes
    construct_names = set()
    for key in user_attrs.keys():
        if key.startswith("construct_") and key.endswith("_aligned"):
            # Extract construct name: "construct_NAME_aligned" -> "NAME"
            construct_name = key.replace("construct_", "").replace("_aligned", "")
            construct_names.add(construct_name)
    
    # Build construct statistics
    for construct_name in construct_names:
        aligned_key = f"construct_{construct_name}_aligned"
        snr_key = f"construct_{construct_name}_snr"
        
        constructs[construct_name] = {
            "aligned_reads": user_attrs.get(aligned_key, 0),
            "signal_to_noise": user_attrs.get(snr_key, 0.0),
        }
    
    return constructs


def generate_detailed_analysis(study, best_trial) -> Dict[str, Any]:
    """Generate detailed analysis of optimization results.
    
    Args:
        study: Optuna study object
        best_trial: Best trial from the study
        
    Returns:
        Dictionary with detailed analysis metrics
    """
    analysis = {
        "best_trial": {
            "number": best_trial.number,
            "quality_score": best_trial.value,
            "metrics": {
                "signal_to_noise": best_trial.user_attrs.get("signal_to_noise", 0.0),
                "alignment_rate": best_trial.user_attrs.get("alignment_rate", 0.0),
                "avg_mapq": best_trial.user_attrs.get("avg_mapq", 0.0),
                "bv_accepted": best_trial.user_attrs.get("bv_accepted", 0),
            },
            "constructs": extract_per_construct_stats(best_trial),
        },
        "summary_stats": {},
        "parameter_analysis": {},
    }
    
    # Calculate summary statistics across all completed trials
    completed_trials = [t for t in study.trials if t.state.name == "COMPLETE"]
    
    if completed_trials:
        quality_scores = [t.value for t in completed_trials]
        snr_values = [t.user_attrs.get("signal_to_noise", 0.0) for t in completed_trials]
        align_rates = [t.user_attrs.get("alignment_rate", 0.0) for t in completed_trials]
        mapq_values = [t.user_attrs.get("avg_mapq", 0.0) for t in completed_trials]
        bv_values = [t.user_attrs.get("bv_accepted", 0) for t in completed_trials]
        
        analysis["summary_stats"] = {
            "total_trials": len(completed_trials),
            "quality_score": {
                "mean": sum(quality_scores) / len(quality_scores),
                "min": min(quality_scores),
                "max": max(quality_scores),
                "std": (sum((x - sum(quality_scores)/len(quality_scores))**2 for x in quality_scores) / len(quality_scores))**0.5 if len(quality_scores) > 1 else 0.0,
            },
            "signal_to_noise": {
                "mean": sum(snr_values) / len(snr_values) if snr_values else 0.0,
                "min": min(snr_values) if snr_values else 0.0,
                "max": max(snr_values) if snr_values else 0.0,
            },
            "alignment_rate": {
                "mean": sum(align_rates) / len(align_rates) if align_rates else 0.0,
                "min": min(align_rates) if align_rates else 0.0,
                "max": max(align_rates) if align_rates else 0.0,
            },
            "avg_mapq": {
                "mean": sum(mapq_values) / len(mapq_values) if mapq_values else 0.0,
                "min": min(mapq_values) if mapq_values else 0.0,
                "max": max(mapq_values) if mapq_values else 0.0,
            },
            "bv_accepted": {
                "mean": sum(bv_values) / len(bv_values) if bv_values else 0.0,
                "min": min(bv_values) if bv_values else 0,
                "max": max(bv_values) if bv_values else 0,
            },
        }
        
        # Analyze parameter frequency in top 10% of trials
        top_n = max(1, len(completed_trials) // 10)
        top_trials = sorted(completed_trials, key=lambda t: t.value, reverse=True)[:top_n]
        
        param_counts = defaultdict(lambda: defaultdict(int))
        for trial in top_trials:
            for param, value in trial.params.items():
                param_counts[param][value] += 1
        
        # Find most common values for each parameter in top trials
        for param, value_counts in param_counts.items():
            if value_counts:
                most_common = max(value_counts.items(), key=lambda x: x[1])
                analysis["parameter_analysis"][param] = {
                    "most_common_value": most_common[0],
                    "frequency": most_common[1],
                    "frequency_pct": (most_common[1] / top_n) * 100,
                }
    
    return analysis


def format_analysis_report(analysis: Dict[str, Any], best_bt2_args: List[str]) -> str:
    """Format analysis results as a human-readable report.
    
    Args:
        analysis: Analysis dictionary from generate_detailed_analysis
        best_bt2_args: Best Bowtie2 arguments as list
        
    Returns:
        Formatted string report
    """
    lines = []
    lines.append("=" * 80)
    lines.append("Detailed Analysis Report")
    lines.append("=" * 80)
    lines.append("")
    
    # Best trial information
    best = analysis["best_trial"]
    lines.append("BEST TRIAL SUMMARY")
    lines.append("-" * 80)
    lines.append(f"Trial Number: {best['number']}")
    lines.append(f"Quality Score: {best['quality_score']:.4f}")
    lines.append("")
    
    # Best trial metrics
    metrics = best["metrics"]
    lines.append("Best Trial Metrics:")
    lines.append(f"  Signal-to-Noise (AC/GU): {metrics['signal_to_noise']:.2f}")
    lines.append(f"  Alignment Rate: {metrics['alignment_rate']:.2%}")
    lines.append(f"  Average MAPQ: {metrics['avg_mapq']:.1f}")
    lines.append(f"  Bit Vectors Accepted: {metrics['bv_accepted']:,}")
    lines.append("")
    
    # Final bit vector metrics (detailed mutation statistics)
    if "final_bit_vector_metrics" in best:
        bv_metrics = best["final_bit_vector_metrics"]
        lines.append("Final Bit Vector Analysis:")
        lines.append(f"  Total Bit Vectors: {bv_metrics.get('total_bit_vectors', 0):,}")
        lines.append(f"  Accepted Bit Vectors: {bv_metrics.get('accepted_bit_vectors', 0):,}")
        lines.append(f"  Rejected (Low MAPQ): {bv_metrics.get('rejected_low_mapq', 0):,}")
        lines.append(f"  Total Mutations: {bv_metrics.get('total_mutations', 0):,}")
        lines.append(f"  Average Mutations per Read: {bv_metrics.get('avg_mutations_per_read', 0.0):.2f}")
        lines.append(f"  Reads with Mutations: {bv_metrics.get('reads_with_mutations', 0):,}")
        lines.append(f"  Reads without Mutations: {bv_metrics.get('reads_without_mutations', 0):,}")
        lines.append(f"  Average Coverage: {bv_metrics.get('avg_coverage', 0.0):.1f} positions")
        
        # Mutation distribution
        mut_dist = bv_metrics.get('mutation_distribution', {})
        if mut_dist:
            lines.append("  Mutation Distribution:")
            sorted_dist = sorted(mut_dist.items(), key=lambda x: int(x[0]) if str(x[0]).isdigit() else 0)
            for mut_count, read_count in sorted_dist[:10]:  # Show top 10
                pct = (read_count / bv_metrics.get('accepted_bit_vectors', 1)) * 100
                lines.append(f"    {mut_count} mutations: {read_count:,} reads ({pct:.1f}%)")
        lines.append("")
    
    # Per-construct statistics
    if best["constructs"]:
        lines.append("Per-Construct Statistics:")
        for construct_name, stats in best["constructs"].items():
            lines.append(f"  {construct_name}:")
            lines.append(f"    Aligned Reads: {stats['aligned_reads']:,}")
            lines.append(f"    Signal-to-Noise: {stats['signal_to_noise']:.2f}")
        lines.append("")
    else:
        lines.append("Per-Construct Statistics: Not available")
        lines.append("")
    
    # Summary statistics
    if analysis["summary_stats"]:
        summary = analysis["summary_stats"]
        lines.append("SUMMARY STATISTICS (All Trials)")
        lines.append("-" * 80)
        lines.append(f"Total Completed Trials: {summary['total_trials']}")
        lines.append("")
        
        # Quality score stats
        qs = summary["quality_score"]
        lines.append("Quality Score:")
        lines.append(f"  Mean: {qs['mean']:.4f}")
        lines.append(f"  Range: {qs['min']:.4f} - {qs['max']:.4f}")
        if qs['std'] > 0:
            lines.append(f"  Std Dev: {qs['std']:.4f}")
        lines.append("")
        
        # Signal-to-noise stats
        snr = summary["signal_to_noise"]
        lines.append("Signal-to-Noise (AC/GU):")
        lines.append(f"  Mean: {snr['mean']:.2f}")
        lines.append(f"  Range: {snr['min']:.2f} - {snr['max']:.2f}")
        lines.append("")
        
        # Alignment rate stats
        ar = summary["alignment_rate"]
        lines.append("Alignment Rate:")
        lines.append(f"  Mean: {ar['mean']:.2%}")
        lines.append(f"  Range: {ar['min']:.2%} - {ar['max']:.2%}")
        lines.append("")
        
        # MAPQ stats
        mapq = summary["avg_mapq"]
        lines.append("Average MAPQ:")
        lines.append(f"  Mean: {mapq['mean']:.1f}")
        lines.append(f"  Range: {mapq['min']:.1f} - {mapq['max']:.1f}")
        lines.append("")
        
        # Bit vectors stats
        bv = summary["bv_accepted"]
        lines.append("Bit Vectors Accepted:")
        lines.append(f"  Mean: {bv['mean']:,.0f}")
        lines.append(f"  Range: {bv['min']:,} - {bv['max']:,}")
        lines.append("")
    
    # Parameter analysis
    if analysis["parameter_analysis"]:
        lines.append("PARAMETER ANALYSIS (Top 10% of Trials)")
        lines.append("-" * 80)
        for param, info in sorted(analysis["parameter_analysis"].items()):
            lines.append(f"  {param}:")
            lines.append(f"    Most Common: {info['most_common_value']}")
            lines.append(f"    Frequency: {info['frequency']}/{summary['total_trials']} ({info['frequency_pct']:.1f}% of top trials)")
        lines.append("")
    
    # Best Bowtie2 arguments
    lines.append("BEST BOWTIE2 ARGUMENTS")
    lines.append("-" * 80)
    lines.append(" ".join(best_bt2_args))
    lines.append("")
    
    return "\n".join(lines)

