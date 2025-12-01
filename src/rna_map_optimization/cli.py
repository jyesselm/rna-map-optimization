"""Command-line interface for RNA-MAP optimization."""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import click
import numpy as np
import pandas as pd

try:
    import optuna
    from optuna.importance import get_param_importances
except ImportError:
    optuna = None
    get_param_importances = None

from rna_map_optimization.optimizer import run_optimization
from rna_map_optimization.utils import (
    format_analysis_report,
    generate_detailed_analysis,
)
from rna_map_optimization.comprehensive_analysis import run_comprehensive_analysis


def find_fastq_files(case_dir: Path) -> tuple[Optional[Path], Optional[Path]]:
    """Find FASTQ files in a case directory.
    
    Args:
        case_dir: Directory containing test case files
        
    Returns:
        Tuple of (fastq1, fastq2) paths, or (None, None) if not found
    """
    fastq1 = None
    fastq2 = None
    
    patterns_r1 = ["*_R1*.fastq*", "*_1.fastq*", "*mate1*.fastq*", "*R1*.fq*", "*1.fq*"]
    patterns_r2 = ["*_R2*.fastq*", "*_2.fastq*", "*mate2*.fastq*", "*R2*.fq*", "*2.fq*"]
    
    for pattern in patterns_r1:
        matches = list(case_dir.glob(pattern))
        if matches:
            fastq1 = matches[0]
            break
    
    for pattern in patterns_r2:
        matches = list(case_dir.glob(pattern))
        if matches:
            fastq2 = matches[0]
            break
    
    return fastq1, fastq2


def find_fasta_file(case_dir: Path) -> Optional[Path]:
    """Find FASTA file in a case directory.
    
    Args:
        case_dir: Directory containing test case files
        
    Returns:
        Path to FASTA file, or None if not found
    """
    patterns = ["*.fasta", "*.fa", "*.fas"]
    for pattern in patterns:
        matches = list(case_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def process_single_case(
    case_dir: Path,
    output_base: Path,
    n_trials: int = 100,
    mapq_cutoff: int = 20,
    threads: int = 4,
    optimize_threads: bool = False,
    timeout: Optional[int] = None,
    study_name: Optional[str] = None,
    storage: Optional[str] = None,
    cleanup: bool = True,
) -> dict:
    """Process a single test case.
    
    Args:
        case_dir: Directory containing test case files
        output_base: Base output directory
        n_trials: Number of Optuna trials
        mapq_cutoff: MAPQ cutoff
        threads: Number of threads
        optimize_threads: Whether to optimize threads
        timeout: Timeout in seconds
        study_name: Name for Optuna study (defaults to case directory name)
        storage: Optuna storage URL
        
    Returns:
        Dictionary with results
    """
    case_name = case_dir.name
    
    fasta = find_fasta_file(case_dir)
    fastq1, fastq2 = find_fastq_files(case_dir)
    
    if not fasta:
        raise FileNotFoundError(f"No FASTA file found in {case_dir}")
    if not fastq1:
        raise FileNotFoundError(f"No FASTQ file found in {case_dir}")
    
    output_dir = output_base / case_name
    
    if study_name is None:
        study_name = f"bowtie2_optimization_{case_name}"
    
    click.echo(f"\n{'='*80}")
    click.echo(f"Processing case: {case_name}")
    click.echo(f"{'='*80}")
    click.echo(f"FASTA: {fasta}")
    click.echo(f"FASTQ1: {fastq1}")
    if fastq2:
        click.echo(f"FASTQ2: {fastq2}")
    click.echo(f"Output: {output_dir}")
    click.echo()
    
    results = run_optimization(
        fasta=fasta,
        fastq1=fastq1,
        fastq2=fastq2,
        output_dir=output_dir,
        n_trials=n_trials,
        mapq_cutoff=mapq_cutoff,
        threads=threads,
        optimize_threads=optimize_threads,
        timeout=timeout,
        study_name=study_name,
        storage=storage,
        cleanup=cleanup,
    )
    
    click.echo(f"\nCase {case_name} complete:")
    click.echo(f"  Best quality score: {results['best_value']:.4f}")
    click.echo(f"  Best signal-to-noise: {results['best_trial'].user_attrs.get('signal_to_noise', 0.0):.2f}")
    click.echo(f"  Alignment rate: {results['best_trial'].user_attrs.get('alignment_rate', 0.0):.2%}")
    click.echo()
    
    return results


@click.group()
@click.version_option()
def cli():
    """RNA-MAP optimization toolkit."""
    pass


@cli.command()
@click.option(
    '--fasta',
    type=click.Path(exists=True, path_type=Path),
    help='Path to reference FASTA file (single case mode with explicit files)',
)
@click.option(
    '--case-dir',
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help='Directory containing a single test case (auto-detects FASTA/FASTQ files)',
)
@click.option(
    '--cases-dir',
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help='Directory containing multiple test cases (batch mode)',
)
@click.option('--fastq1', type=click.Path(exists=True, path_type=Path), help='Path to first FASTQ file')
@click.option('--fastq2', type=click.Path(exists=True, path_type=Path), help='Path to second FASTQ file')
@click.option(
    '--output-dir',
    type=click.Path(path_type=Path),
    default='bowtie2_optimization_optuna',
    help='Output directory for results',
)
@click.option('--n-trials', type=int, default=100, help='Number of Optuna trials to run (overrides any config file setting)')
@click.option('--mapq-cutoff', type=int, default=20, help='MAPQ cutoff for high-quality alignments')
@click.option('--threads', type=int, default=4, help='Number of threads for alignment')
@click.option('--optimize-threads', is_flag=True, help='Optimize number of threads')
@click.option('--timeout', type=int, help='Timeout in seconds')
@click.option('--study-name', type=str, help='Name for Optuna study')
@click.option('--storage', type=str, help='Optuna storage URL (e.g., sqlite:///study.db)')
@click.option('--keep-intermediates', is_flag=True, help='Keep trial directories and SAM files (default: cleanup to save space)')
def optimize(
    fasta: Optional[Path],
    case_dir: Optional[Path],
    cases_dir: Optional[Path],
    fastq1: Optional[Path],
    fastq2: Optional[Path],
    output_dir: Path,
    n_trials: int,
    mapq_cutoff: int,
    threads: int,
    optimize_threads: bool,
    timeout: Optional[int],
    study_name: Optional[str],
    storage: Optional[str],
    keep_intermediates: bool,
):
    """Run parameter optimization."""
    if not any([fasta, case_dir, cases_dir]):
        click.echo("ERROR: Must specify one of --fasta, --case-dir, or --cases-dir", err=True)
        sys.exit(1)
    
    if case_dir:
        click.echo("=" * 80)
        click.echo("Bowtie2 Parameter Optimization with Optuna (Single Case)")
        click.echo("=" * 80)
        click.echo(f"Case directory: {case_dir}")
        click.echo(f"Output directory: {output_dir}")
        click.echo(f"Number of trials: {n_trials}")
        click.echo(f"Threads: {threads} ({'optimizing' if optimize_threads else 'fixed'})")
        click.echo()
        
        try:
            results = process_single_case(
                case_dir=case_dir,
                output_base=output_dir,
                n_trials=n_trials,
                mapq_cutoff=mapq_cutoff,
                threads=threads,
                optimize_threads=optimize_threads,
                timeout=timeout,
                study_name=study_name,
                storage=storage,
                cleanup=not keep_intermediates,
            )
            
            analysis = generate_detailed_analysis(results['study'], results['best_trial'])
            
            # Add final bit vector metrics if available
            if 'final_bit_vector_metrics' in results:
                analysis['best_trial']['final_bit_vector_metrics'] = results['final_bit_vector_metrics']
            
            report = format_analysis_report(analysis, results['best_bt2_args'])
            click.echo("\n" + report)
            
            output_path = output_dir / case_dir.name
            click.echo(f"Results saved to: {output_path}")
            click.echo()
        except Exception as e:
            click.echo(f"ERROR: Failed to process case: {e}", err=True)
            sys.exit(1)
    
    elif fasta:
        if not fastq1:
            click.echo("ERROR: --fastq1 is required when using --fasta", err=True)
            sys.exit(1)
        
        click.echo("=" * 80)
        click.echo("Bowtie2 Parameter Optimization with Optuna")
        click.echo("=" * 80)
        click.echo(f"FASTA: {fasta}")
        click.echo(f"FASTQ1: {fastq1}")
        if fastq2:
            click.echo(f"FASTQ2: {fastq2}")
        click.echo(f"Output directory: {output_dir}")
        click.echo(f"Number of trials: {n_trials}")
        click.echo(f"Threads: {threads} ({'optimizing' if optimize_threads else 'fixed'})")
        click.echo()
        
        results = run_optimization(
            fasta=fasta,
            fastq1=fastq1,
            fastq2=fastq2,
            output_dir=output_dir,
            n_trials=n_trials,
            mapq_cutoff=mapq_cutoff,
            threads=threads,
            optimize_threads=optimize_threads,
            timeout=timeout,
            study_name=study_name or "bowtie2_optimization",
            storage=storage,
            cleanup=not keep_intermediates,
        )
        
        analysis = generate_detailed_analysis(results['study'], results['best_trial'])
        
        # Add final bit vector metrics if available
        if 'final_bit_vector_metrics' in results:
            analysis['best_trial']['final_bit_vector_metrics'] = results['final_bit_vector_metrics']
        
        report = format_analysis_report(analysis, results['best_bt2_args'])
        click.echo("\n" + report)
        
        click.echo(f"Results saved to: {output_dir}")
        click.echo()
    
    elif cases_dir:
        click.echo("=" * 80)
        click.echo("Batch Bowtie2 Parameter Optimization with Optuna")
        click.echo("=" * 80)
        click.echo(f"Cases directory: {cases_dir}")
        click.echo(f"Output directory: {output_dir}")
        click.echo(f"Number of trials per case: {n_trials}")
        click.echo(f"Threads: {threads} ({'optimizing' if optimize_threads else 'fixed'})")
        click.echo()
        
        case_dirs = [d for d in cases_dir.iterdir() if d.is_dir()]
        
        if not case_dirs:
            click.echo(f"No case directories found in {cases_dir}", err=True)
            sys.exit(1)
        
        click.echo(f"Found {len(case_dirs)} test cases:")
        for case_dir in case_dirs:
            click.echo(f"  - {case_dir.name}")
        click.echo()
        
        all_results = {}
        for i, case_dir in enumerate(case_dirs, 1):
            try:
                click.echo(f"\n[{i}/{len(case_dirs)}] Processing {case_dir.name}...")
                results = process_single_case(
                    case_dir=case_dir,
                    output_base=output_dir,
                    n_trials=n_trials,
                    mapq_cutoff=mapq_cutoff,
                    threads=threads,
                    optimize_threads=optimize_threads,
                    timeout=timeout,
                    study_name=None,
                    storage=storage,
                    cleanup=not keep_intermediates,
                )
                all_results[case_dir.name] = results
            except Exception as e:
                click.echo(f"ERROR: Failed to process {case_dir.name}: {e}", err=True)
                continue
        
        click.echo("\n" + "=" * 80)
        click.echo("Batch Optimization Complete")
        click.echo("=" * 80)
        click.echo(f"Successfully processed {len(all_results)}/{len(case_dirs)} cases")
        click.echo()
        
        if all_results:
            click.echo("Summary of best results:")
            for case_name, results in all_results.items():
                best_snr = results['best_trial'].user_attrs.get("signal_to_noise", 0.0)
                best_align_rate = results['best_trial'].user_attrs.get("alignment_rate", 0.0)
                click.echo(f"  {case_name}:")
                click.echo(f"    Quality score: {results['best_value']:.4f}")
                click.echo(f"    Signal-to-Noise: {best_snr:.2f}")
                click.echo(f"    Alignment Rate: {best_align_rate:.2%}")
            click.echo()
        
        click.echo(f"All results saved to: {output_dir}")
        click.echo()


def fasta_to_dict(fasta_path: Path) -> Dict[str, str]:
    """Read FASTA file and return dictionary of sequence names to sequences."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def longest_common_substring(
    seq1: str,
    seq2: str,
    min_length: int = 20
) -> Tuple[str, int, int] | None:
    """Find longest common substring between two sequences."""
    m, n = len(seq1), len(seq2)
    max_len = 0
    end_pos1 = 0
    end_pos2 = 0
    
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
                if dp[i][j] > max_len:
                    max_len = dp[i][j]
                    end_pos1 = i
                    end_pos2 = j
            else:
                dp[i][j] = 0
    
    if max_len >= min_length:
        start_pos1 = end_pos1 - max_len
        start_pos2 = end_pos2 - max_len
        substring = seq1[start_pos1:end_pos1]
        return (substring, start_pos1, start_pos2)
    
    return None


def find_common_substrings(
    sequences: Dict[str, str],
    min_length: int = 20
) -> List[Tuple[str, int, int, Set[str]]]:
    """Find common substrings across sequences."""
    seq_list = list(sequences.items())
    common_regions = []
    
    if len(seq_list) < 2:
        return common_regions
    
    for i, (name1, seq1) in enumerate(seq_list):
        for j, (name2, seq2) in enumerate(seq_list[i+1:], start=i+1):
            lcs = longest_common_substring(seq1, seq2, min_length)
            if lcs:
                substring, pos1, pos2 = lcs
                common_regions.append((substring, pos1, pos2, {name1, name2}))
    
    merged = []
    seen = set()
    for substring, pos1, pos2, seq_names in common_regions:
        key = (pos1, pos2, tuple(sorted(seq_names)))
        if key in seen:
            continue
        seen.add(key)
        merged.append((substring, pos1, pos2, seq_names))
    
    return merged


def find_unique_regions(
    sequences: Dict[str, str],
    common_regions: List[Tuple[str, int, int, Set[str]]],
    min_unique_length: int = 50
) -> Dict[str, List[Tuple[int, int]]]:
    """Find unique regions for each sequence."""
    unique_regions = {}
    
    for seq_name, seq in sequences.items():
        affected_regions = []
        for substring, pos1, pos2, seq_names in common_regions:
            if seq_name in seq_names:
                if seq_name == list(seq_names)[0]:
                    affected_regions.append((pos1, pos1 + len(substring)))
                else:
                    affected_regions.append((pos2, pos2 + len(substring)))
        
        affected_regions.sort()
        
        unique = []
        start = 0
        
        for common_start, common_end in affected_regions:
            if common_start > start:
                unique_len = common_start - start
                if unique_len >= min_unique_length:
                    unique.append((start, common_start))
            start = max(start, common_end)
        
        if start < len(seq):
            unique_len = len(seq) - start
            if unique_len >= min_unique_length:
                unique.append((start, len(seq)))
        
        unique_regions[seq_name] = unique
    
    return unique_regions


def remove_common_sequences(
    sequences: Dict[str, str],
    min_unique_length: int = 50
) -> Dict[str, str]:
    """Remove common sequences, keeping only unique regions."""
    common_regions = find_common_substrings(sequences, min_length=20)
    unique_regions = find_unique_regions(sequences, common_regions, min_unique_length)
    
    result = {}
    for seq_name, seq in sequences.items():
        if seq_name in unique_regions and unique_regions[seq_name]:
            unique_parts = [seq[start:end] for start, end in unique_regions[seq_name]]
            result[seq_name] = 'N' * 10 + 'N'.join(unique_parts) + 'N' * 10
        else:
            result[seq_name] = seq
    
    return result


def write_fasta(sequences: Dict[str, str], output_path: Path):
    """Write sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


@cli.command()
@click.option('--fasta', type=click.Path(exists=True, path_type=Path), required=True, help='Input FASTA file with multiple sequences')
@click.option('--output', type=click.Path(path_type=Path), help='Output FASTA file with common sequences removed')
@click.option('--min-common-length', type=int, default=20, help='Minimum length of common sequence to detect')
@click.option('--min-unique-length', type=int, default=50, help='Minimum length of unique region to keep')
@click.option('--analyze-only', is_flag=True, help='Only analyze, don\'t create output file')
def analyze_sequences(
    fasta: Path,
    output: Optional[Path],
    min_common_length: int,
    min_unique_length: int,
    analyze_only: bool,
):
    """Analyze common sequences in FASTA file."""
    click.echo(f"Reading sequences from {fasta}")
    sequences = fasta_to_dict(fasta)
    click.echo(f"Found {len(sequences)} sequences")
    
    for name, seq in sequences.items():
        click.echo(f"  - {name}: {len(seq)} bp")
    
    if len(sequences) < 2:
        click.echo("\nWarning: Need at least 2 sequences to find common regions", err=True)
        sys.exit(1)
    
    click.echo(f"\nAnalyzing common sequences (min length: {min_common_length} bp)...")
    common_regions = find_common_substrings(sequences, min_length=min_common_length)
    
    if not common_regions:
        click.echo("No common sequences found!")
        return
    
    click.echo(f"\nFound {len(common_regions)} common region(s):")
    for i, (substring, pos1, pos2, seq_names) in enumerate(common_regions, 1):
        click.echo(f"\n  Common region {i}:")
        click.echo(f"    Length: {len(substring)} bp")
        click.echo(f"    Sequences: {', '.join(sorted(seq_names))}")
        if len(substring) > 50:
            click.echo(f"    Preview: {substring[:50]}...")
        else:
            click.echo(f"    Sequence: {substring}")
    
    total_common = sum(len(substring) for substring, _, _, _ in common_regions)
    total_length = sum(len(seq) for seq in sequences.values())
    avg_length = total_length / len(sequences)
    common_percent = (total_common / avg_length) * 100 if avg_length > 0 else 0
    
    click.echo(f"\nStatistics:")
    click.echo(f"  Total common sequence length: {total_common} bp")
    click.echo(f"  Average sequence length: {avg_length:.1f} bp")
    click.echo(f"  Common sequence percentage: {common_percent:.1f}%")
    
    click.echo(f"\nFinding unique regions (min length: {min_unique_length} bp)...")
    unique_regions = find_unique_regions(sequences, common_regions, min_unique_length)
    
    click.echo("\nUnique regions per sequence:")
    for seq_name, regions in unique_regions.items():
        if regions:
            total_unique = sum(end - start for start, end in regions)
            click.echo(f"  {seq_name}: {len(regions)} region(s), {total_unique} bp total")
            for start, end in regions:
                click.echo(f"    [{start}:{end}] ({end-start} bp)")
        else:
            click.echo(f"  {seq_name}: No unique regions found (all common)")
    
    if not analyze_only and output:
        click.echo(f"\nCreating output file with common sequences removed: {output}")
        unique_sequences = remove_common_sequences(sequences, min_unique_length)
        write_fasta(unique_sequences, output)
        
        click.echo("\nOutput sequences:")
        for name, seq in unique_sequences.items():
            click.echo(f"  - {name}: {len(seq)} bp (original: {len(sequences[name])} bp)")
        
        click.echo(f"\n✓ Saved to {output}")
        click.echo("\nRecommendation: Test alignment with this file and compare:")
        click.echo("  - MAPQ scores (should be higher)")
        click.echo("  - Multi-mapping rate (should be lower)")
        click.echo("  - Alignment rate (may be slightly lower if reads span common regions)")


def load_case_results(case_dir: Path, top_n: int = 100) -> Dict[str, Any] | None:
    """Load top N results from a single case."""
    # Try Optuna format first (new format)
    summary_file = case_dir / "optuna_summary.csv"
    best_params_file = case_dir / "optuna_study.json"
    
    # Fall back to old grid search format
    if not summary_file.exists():
        summary_file = case_dir / "optimization_summary.csv"
        best_params_file = case_dir / "best_parameters.json"
    
    if not summary_file.exists():
        return None
    
    df = pd.read_csv(summary_file)
    
    # Handle Optuna format (has 'trial' column) vs old format (has 'combo_id')
    if "trial" in df.columns:
        # Optuna format - filter out failed trials (quality_score == 0)
        valid_trials = df[df["quality_score"] > 0].copy()
        top_results = valid_trials.nlargest(top_n, "quality_score")
        baseline_quality_score = None  # Optuna doesn't have baseline
    elif "combo_id" in df.columns:
        # Old grid search format
        non_baseline = df[df["combo_id"] != 0].copy()
        top_results = non_baseline.nlargest(top_n, "quality_score")
        baseline_df = df[df.get("is_baseline", pd.Series([False] * len(df))) == True]
        baseline_quality_score = float(baseline_df["quality_score"].iloc[0]) if len(baseline_df) > 0 else None
    else:
        # Unknown format, just take top N
        top_results = df.nlargest(top_n, "quality_score")
        baseline_quality_score = None
    
    best_params = None
    if best_params_file.exists():
        with open(best_params_file, "r") as f:
            data = json.load(f)
            # Extract best parameters from Optuna format
            if "best_trial" in data:
                best_params = {
                    "params": data["best_trial"].get("params", {}),
                    "bowtie2_args": data.get("best_bowtie2_args", ""),
                    "quality_score": data["best_trial"].get("value", 0),
                }
            else:
                best_params = data
    
    return {
        "case_name": case_dir.name,
        "total_combinations": len(df),
        "baseline_quality_score": baseline_quality_score,
        "top_n": top_n,
        "top_results": top_results.to_dict("records"),
        "best_parameters": best_params,
    }


def aggregate_results(all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate results across all cases."""
    if not all_results:
        return {"error": "No results to aggregate"}
    
    all_top_results = []
    for case_result in all_results:
        case_name = case_result["case_name"]
        for result in case_result["top_results"]:
            result["case_name"] = case_name
            all_top_results.append(result)
    
    combined_df = pd.DataFrame(all_top_results)
    
    # Handle missing columns gracefully
    def safe_mean(col_name, default=0.0):
        if col_name in combined_df.columns:
            return float(combined_df[col_name].mean())
        return default
    
    def safe_max(col_name, default=0.0):
        if col_name in combined_df.columns:
            return float(combined_df[col_name].max())
        return default
    
    def safe_min(col_name, default=0.0):
        if col_name in combined_df.columns:
            return float(combined_df[col_name].min())
        return default
    
    # Handle signal_to_noise which might have inf values
    snr_values = combined_df["signal_to_noise"] if "signal_to_noise" in combined_df.columns else pd.Series([0.0])
    snr_clean = snr_values.replace([float("inf"), float("-inf")], None)
    
    stats = {
        "total_cases": len(all_results),
        "total_combinations": len(combined_df),
        "avg_quality_score": safe_mean("quality_score"),
        "max_quality_score": safe_max("quality_score"),
        "min_quality_score": safe_min("quality_score"),
        "avg_signal_to_noise": float(snr_clean.mean()) if len(snr_clean) > 0 else 0.0,
        "avg_alignment_rate": safe_mean("alignment_rate"),
        "avg_high_quality_rate": safe_mean("high_quality_rate"),
    }
    
    best_overall = combined_df.nlargest(1, "quality_score").iloc[0].to_dict()
    
    best_per_case = {}
    for case_result in all_results:
        case_name = case_result["case_name"]
        if case_result["top_results"]:
            best = max(case_result["top_results"], key=lambda x: x["quality_score"])
            # Handle both Optuna format (trial) and old format (combo_id)
            trial_id = best.get("trial", best.get("combo_id", 0))
            best_per_case[case_name] = {
                "trial_id": int(trial_id) if trial_id is not None else 0,
                "combo_id": int(trial_id) if trial_id is not None else 0,  # For backward compatibility
                "quality_score": float(best["quality_score"]),
                "signal_to_noise": float(best["signal_to_noise"]) if best.get("signal_to_noise") != "inf" and best.get("signal_to_noise") is not None else None,
                "alignment_rate": float(best.get("alignment_rate", 0.0)),
            }
    
    return {
        "statistics": stats,
        "best_overall": best_overall,
        "best_per_case": best_per_case,
        "all_results": all_results,
    }


@cli.command()
@click.option('--results-dir', type=click.Path(exists=True, path_type=Path), default='optimization_results', help='Base directory containing optimization results')
@click.option('--top-n', type=int, default=100, help='Number of top results to extract from each case')
@click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory for aggregated results')
@click.option('--min-quality', type=float, help='Minimum quality score threshold')
def collect_results(
    results_dir: Path,
    top_n: int,
    output_dir: Optional[Path],
    min_quality: Optional[float],
):
    """Collect top results from optimization runs."""
    if not results_dir.exists():
        click.echo(f"ERROR: Results directory not found: {results_dir}", err=True)
        sys.exit(1)
    
    if output_dir is None:
        output_dir = results_dir / "aggregated"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    click.echo("=" * 80)
    click.echo("Collecting Top Results from Optimization Runs")
    click.echo("=" * 80)
    click.echo(f"Results directory: {results_dir}")
    click.echo(f"Top N per case: {top_n}")
    click.echo(f"Output directory: {output_dir}")
    click.echo("")
    
    # Check if this is a case directory itself (has optuna_summary.csv or optimization_summary.csv)
    case_summary = results_dir / "optuna_summary.csv"
    if not case_summary.exists():
        case_summary = results_dir / "optimization_summary.csv"
    
    if case_summary.exists():
        # This is a case directory, process it directly
        case_dirs = [results_dir]
    else:
        # This is a parent directory, look for case subdirectories
        case_dirs = [
            d for d in results_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".") and d.name != "aggregated" and d.name != "job_scripts"
            and d.name not in ["index", "results", "visualizations"]  # Skip Optuna output subdirs
        ]
    
    if not case_dirs:
        click.echo(f"ERROR: No case directories found in {results_dir}", err=True)
        click.echo("Expected either:")
        click.echo("  1. A case directory containing optuna_summary.csv")
        click.echo("  2. A parent directory containing case subdirectories")
        sys.exit(1)
    
    click.echo(f"Found {len(case_dirs)} case directory(ies)")
    click.echo("")
    
    all_results = []
    for case_dir in sorted(case_dirs):
        click.echo(f"Loading results from {case_dir.name}...")
        case_result = load_case_results(case_dir, top_n=top_n)
        
        if case_result is None:
            click.echo(f"  ⚠️  Skipping {case_dir.name}: optuna_summary.csv or optimization_summary.csv not found")
            continue
        
        if min_quality is not None:
            filtered = [
                r for r in case_result["top_results"]
                if r["quality_score"] >= min_quality
            ]
            case_result["top_results"] = filtered
            case_result["top_n"] = len(filtered)
        
        all_results.append(case_result)
        click.echo(f"  ✓ Loaded {len(case_result['top_results'])} top results")
        if case_result["baseline_quality_score"] is not None:
            click.echo(f"    Baseline quality score: {case_result['baseline_quality_score']:.4f}")
        click.echo("")
    
    if not all_results:
        click.echo("ERROR: No valid results found", err=True)
        sys.exit(1)
    
    click.echo("Aggregating results...")
    aggregated = aggregate_results(all_results)
    
    click.echo("Saving aggregated results...")
    
    aggregated_json = output_dir / "aggregated_top_results.json"
    with open(aggregated_json, "w") as f:
        json.dump(aggregated, f, indent=2, default=str)
    click.echo(f"  ✓ Saved: {aggregated_json}")
    
    all_top_records = []
    for case_result in all_results:
        for record in case_result["top_results"]:
            all_top_records.append(record)
    
    if all_top_records:
        combined_df = pd.DataFrame(all_top_records)
        combined_df = combined_df.sort_values("quality_score", ascending=False)
        
        combined_csv = output_dir / "combined_top_results.csv"
        combined_df.to_csv(combined_csv, index=False)
        click.echo(f"  ✓ Saved: {combined_csv}")
        
        top_100_overall = combined_df.head(100)
        top_100_csv = output_dir / "top_100_overall.csv"
        top_100_overall.to_csv(top_100_csv, index=False)
        click.echo(f"  ✓ Saved: {top_100_csv}")
    
    summary = {
        "total_cases": aggregated["statistics"]["total_cases"],
        "total_combinations": aggregated["statistics"]["total_combinations"],
        "statistics": aggregated["statistics"],
        "best_per_case": aggregated["best_per_case"],
        "best_overall": {
            "case_name": aggregated["best_overall"].get("case_name"),
            "combo_id": int(aggregated["best_overall"].get("combo_id", 0)),
            "quality_score": float(aggregated["best_overall"].get("quality_score", 0)),
            "signal_to_noise": aggregated["best_overall"].get("signal_to_noise"),
            "alignment_rate": float(aggregated["best_overall"].get("alignment_rate", 0)),
        },
    }
    
    summary_json = output_dir / "summary_statistics.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    click.echo(f"  ✓ Saved: {summary_json}")
    
    click.echo("")
    click.echo("=" * 80)
    click.echo("Summary")
    click.echo("=" * 80)
    click.echo(f"Total cases processed: {aggregated['statistics']['total_cases']}")
    click.echo(f"Total combinations: {aggregated['statistics']['total_combinations']}")
    click.echo(f"Average quality score: {aggregated['statistics']['avg_quality_score']:.4f}")
    click.echo(f"Best quality score: {aggregated['statistics']['max_quality_score']:.4f}")
    click.echo(f"Best case: {aggregated['best_overall'].get('case_name', 'N/A')}")
    # Try to get trial or combo_id from best_overall
    trial_or_combo = aggregated['best_overall'].get('trial', aggregated['best_overall'].get('combo_id', aggregated['best_overall'].get('trial_id', 'N/A')))
    click.echo(f"Best trial/combo_id: {trial_or_combo}")
    click.echo("")
    click.echo("Best per case:")
    for case_name, best in aggregated["best_per_case"].items():
        trial_or_combo = best.get('trial_id', best.get('combo_id', 'N/A'))
        click.echo(f"  {case_name}: quality_score={best['quality_score']:.4f}, "
              f"trial/combo_id={trial_or_combo}")
    click.echo("")
    click.echo(f"Results saved to: {output_dir}")


def analyze_parameter_importance(
    results_dir: Path,
    top_n: int = 100,
    output_file: Optional[Path] = None,
) -> Dict[str, Any]:
    """Analyze parameter importance from optimization results.
    
    Args:
        results_dir: Directory containing optimization results
        top_n: Number of top trials to analyze
        output_file: Optional file to save analysis results
        
    Returns:
        Dictionary with analysis results
    """
    # Load results
    summary_file = results_dir / "optuna_summary.csv"
    study_file = results_dir / "optuna_study.pkl"
    
    if not summary_file.exists():
        click.echo(f"ERROR: {summary_file} not found", err=True)
        return {}
    
    df = pd.read_csv(summary_file)
    
    # Filter to top N trials
    valid_trials = df[df["quality_score"] > 0].copy()
    top_trials = valid_trials.nlargest(top_n, "quality_score")
    
    if len(top_trials) == 0:
        click.echo("ERROR: No valid trials found", err=True)
        return {}
    
    # Get parameter columns (exclude metrics)
    metric_cols = ["trial", "quality_score", "signal_to_noise", "alignment_rate", "avg_mapq", "bv_accepted"]
    param_cols = [col for col in df.columns if col not in metric_cols]
    
    analysis = {
        "total_trials": len(valid_trials),
        "top_n_analyzed": len(top_trials),
        "parameter_importance": {},
        "parameter_correlations": {},
        "parameter_variance": {},
        "top_value_frequency": {},
        "recommendations": {},
    }
    
    # 1. Parameter importance using correlation with quality_score
    click.echo("Calculating parameter correlations with quality score...")
    correlations = {}
    for param in param_cols:
        if param in top_trials.columns:
            # Handle categorical parameters
            if top_trials[param].dtype == 'object':
                # Convert to numeric for correlation (use value counts)
                param_encoded = pd.get_dummies(top_trials[param], prefix=param)
                if len(param_encoded.columns) > 0:
                    corr = param_encoded.corrwith(top_trials["quality_score"]).abs().max()
                    correlations[param] = float(corr) if not pd.isna(corr) else 0.0
            else:
                corr = top_trials[param].corr(top_trials["quality_score"])
                correlations[param] = abs(float(corr)) if not pd.isna(corr) else 0.0
    
    # Sort by importance
    sorted_corrs = sorted(correlations.items(), key=lambda x: x[1], reverse=True)
    analysis["parameter_importance"] = dict(sorted_corrs)
    
    # 2. Parameter variance in top results (high variance = worth exploring more)
    click.echo("Analyzing parameter variance in top results...")
    variance_analysis = {}
    for param in param_cols:
        if param in top_trials.columns:
            param_series = top_trials[param].dropna()
            if len(param_series) == 0:
                continue
                
            if param_series.dtype == 'object' or param_series.dtype.name == 'category':
                # For categorical, calculate entropy-based diversity
                value_counts = param_series.value_counts()
                total_count = len(param_series)
                unique_count = len(value_counts)
                
                # Calculate normalized entropy (diversity measure)
                if unique_count > 1:
                    probs = value_counts / total_count
                    entropy = -sum(p * np.log2(p) for p in probs if p > 0)
                    max_entropy = np.log2(unique_count)
                    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0
                else:
                    normalized_entropy = 0.0
                
                variance_analysis[param] = {
                    "unique_values": unique_count,
                    "total_values": total_count,
                    "diversity": float(normalized_entropy),
                    "entropy": float(entropy) if unique_count > 1 else 0.0,
                }
            else:
                # For numeric, calculate coefficient of variation
                param_values = top_trials[param].dropna()
                if len(param_values) > 0 and param_values.std() > 0:
                    cv = param_values.std() / param_values.mean()
                    variance_analysis[param] = {
                        "mean": float(param_values.mean()),
                        "std": float(param_values.std()),
                        "cv": float(cv),
                        "min": float(param_values.min()),
                        "max": float(param_values.max()),
                    }
                else:
                    variance_analysis[param] = {
                        "mean": float(param_values.mean()) if len(param_values) > 0 else 0.0,
                        "std": 0.0,
                        "cv": 0.0,
                    }
    
    analysis["parameter_variance"] = variance_analysis
    
    # 3. Most common values in top results
    click.echo("Analyzing most common parameter values in top results...")
    top_value_freq = {}
    for param in param_cols:
        if param in top_trials.columns:
            value_counts = top_trials[param].value_counts()
            if len(value_counts) > 0:
                top_value_freq[param] = {
                    "most_common": value_counts.index[0] if len(value_counts) > 0 else None,
                    "frequency": int(value_counts.iloc[0]),
                    "frequency_pct": float((value_counts.iloc[0] / len(top_trials)) * 100),
                    "all_values": value_counts.to_dict(),
                }
    
    analysis["top_value_frequency"] = top_value_freq
    
    # 4. Try Optuna's built-in importance if study file exists
    if study_file.exists() and optuna is not None:
        try:
            with open(study_file, "rb") as f:
                import pickle
                study = pickle.load(f)
            
            click.echo("Calculating Optuna parameter importance...")
            try:
                optuna_importance = get_param_importances(study)
                analysis["optuna_importance"] = {k: float(v) for k, v in optuna_importance.items()}
            except Exception as e:
                click.echo(f"  Warning: Could not calculate Optuna importance: {e}")
        except Exception as e:
            click.echo(f"  Warning: Could not load study file: {e}")
    
    # 5. Generate recommendations
    click.echo("Generating recommendations...")
    recommendations = {
        "high_importance": [],
        "low_importance": [],
        "explore_more": [],
        "explore_less": [],
        "stable_parameters": [],
    }
    
    # High importance (top 25% by correlation)
    n_high = max(1, len(sorted_corrs) // 4)
    recommendations["high_importance"] = [p[0] for p in sorted_corrs[:n_high]]
    
    # Low importance (bottom 25%)
    recommendations["low_importance"] = [p[0] for p in sorted_corrs[-n_high:]]
    
    # Parameters to explore more (high importance + low diversity)
    for param, corr in sorted_corrs[:n_high]:
        var_info = variance_analysis.get(param, {})
        diversity = var_info.get("diversity", 0.0)
        cv = var_info.get("cv", 0.0)
        
        # Low diversity means we should explore more
        if diversity < 0.3 or (isinstance(cv, float) and cv < 0.1):
            recommendations["explore_more"].append(param)
    
    # Parameters to explore less (low importance + high diversity)
    for param, corr in sorted_corrs[-n_high:]:
        var_info = variance_analysis.get(param, {})
        diversity = var_info.get("diversity", 1.0)
        cv = var_info.get("cv", 1.0)
        
        if diversity > 0.7 or (isinstance(cv, float) and cv > 0.5):
            recommendations["explore_less"].append(param)
    
    # Stable parameters (same value in >80% of top results)
    for param, freq_info in top_value_freq.items():
        if freq_info["frequency_pct"] > 80:
            recommendations["stable_parameters"].append({
                "parameter": param,
                "value": freq_info["most_common"],
                "frequency_pct": freq_info["frequency_pct"],
            })
    
    analysis["recommendations"] = recommendations
    
    # Save if output file specified
    if output_file:
        with open(output_file, "w") as f:
            json.dump(analysis, f, indent=2, default=str)
        click.echo(f"Analysis saved to {output_file}")
    
    return analysis


@cli.command()
@click.option('--results-dir', type=click.Path(exists=True, path_type=Path), required=True, help='Directory containing optimization results')
@click.option('--top-n', type=int, default=100, help='Number of top trials to analyze')
@click.option('--output', type=click.Path(path_type=Path), help='Output JSON file for analysis results')
def analyze_parameters(
    results_dir: Path,
    top_n: int,
    output: Optional[Path],
):
    """Analyze parameter importance and generate recommendations."""
    click.echo("=" * 80)
    click.echo("Parameter Importance Analysis")
    click.echo("=" * 80)
    click.echo(f"Results directory: {results_dir}")
    click.echo(f"Analyzing top {top_n} trials")
    click.echo("")
    
    analysis = analyze_parameter_importance(results_dir, top_n, output)
    
    if not analysis:
        click.echo("ERROR: Could not perform analysis", err=True)
        sys.exit(1)
    
    # Display results
    click.echo("=" * 80)
    click.echo("PARAMETER IMPORTANCE (Correlation with Quality Score)")
    click.echo("=" * 80)
    importance = analysis["parameter_importance"]
    sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=True)
    
    for i, (param, score) in enumerate(sorted_importance[:15], 1):
        click.echo(f"{i:2d}. {param:25s} {score:.4f}")
    click.echo("")
    
    # Recommendations
    recs = analysis["recommendations"]
    click.echo("=" * 80)
    click.echo("RECOMMENDATIONS")
    click.echo("=" * 80)
    
    click.echo("\nHigh Importance Parameters (focus optimization here):")
    for param in recs["high_importance"][:10]:
        importance_score = importance.get(param, 0.0)
        click.echo(f"  • {param} (importance: {importance_score:.4f})")
    
    click.echo("\nParameters to Explore More (high importance, low diversity):")
    for param in recs["explore_more"][:10]:
        var_info = analysis["parameter_variance"].get(param, {})
        diversity = var_info.get("diversity", 0.0)
        click.echo(f"  • {param} (diversity: {diversity:.2%})")
    
    click.echo("\nStable Parameters (same value in >80% of top results):")
    for item in recs["stable_parameters"][:10]:
        click.echo(f"  • {item['parameter']} = {item['value']} ({item['frequency_pct']:.1f}% of top trials)")
    
    click.echo("\nLow Importance Parameters (can reduce search space):")
    for param in recs["low_importance"][:10]:
        importance_score = importance.get(param, 0.0)
        click.echo(f"  • {param} (importance: {importance_score:.4f})")
    
    click.echo("")
    
    # Parameter variance summary
    click.echo("=" * 80)
    click.echo("PARAMETER DIVERSITY IN TOP RESULTS")
    click.echo("=" * 80)
    click.echo("(High diversity = many different values work well)")
    click.echo("(Low diversity = specific values are best)")
    click.echo("")
    
    diversity_scores = []
    for param, var_info in analysis["parameter_variance"].items():
        diversity = var_info.get("diversity", 0.0)
        cv = var_info.get("cv", 0.0)
        # Use diversity for categorical, CV for numeric, normalize both to 0-1
        if isinstance(diversity, float) and diversity > 0:
            # Already normalized entropy (0-1)
            score = diversity
        elif isinstance(cv, float) and cv > 0:
            # Normalize CV (typically 0-1, but cap at 1.0)
            score = min(cv, 1.0)
        else:
            score = 0.0
        diversity_scores.append((param, score))
    
    diversity_scores.sort(key=lambda x: x[1], reverse=True)
    
    click.echo("Most Diverse (explore more values):")
    for param, div in diversity_scores[:10]:
        click.echo(f"  • {param:25s} {div:.4f}")
    
    click.echo("\nLeast Diverse (specific values work best):")
    for param, div in diversity_scores[-10:]:
        click.echo(f"  • {param:25s} {div:.4f}")
    
    click.echo("")
    
    if output:
        click.echo(f"Full analysis saved to: {output}")
    click.echo("")


@cli.command()
@click.option('--results-dir', type=click.Path(exists=True, path_type=Path), required=True, help='Directory containing optimization results (with optuna_summary.csv)')
@click.option('--output-dir', type=click.Path(path_type=Path), help='Output directory for analysis (default: results-dir/comprehensive_analysis)')
def analyze_all_runs(
    results_dir: Path,
    output_dir: Optional[Path],
):
    """Comprehensive analysis of all optimization runs.
    
    This command analyzes all runs (not just top N) to identify:
    - Features that correlate with good/bad scores
    - Parameter distributions in top vs bottom runs
    - Correlation heatmaps and scatter plots
    - Comprehensive HTML report
    """
    click.echo("=" * 80)
    click.echo("Comprehensive Analysis of All Runs")
    click.echo("=" * 80)
    click.echo(f"Results directory: {results_dir}")
    click.echo("")
    
    try:
        results = run_comprehensive_analysis(results_dir, output_dir)
        
        if "error" in results:
            click.echo(f"ERROR: {results['error']}", err=True)
            sys.exit(1)
        
        output_path = Path(results["output_dir"])
        click.echo("Analysis complete!")
        click.echo("")
        click.echo("Generated files:")
        click.echo(f"  ✓ {output_path / 'report.html'} - Complete HTML report with all plots embedded")
        click.echo(f"  ✓ {output_path / 'correlations.csv'} - Parameter correlations")
        click.echo(f"  ✓ {output_path / 'comparison.json'} - Good vs bad run comparison")
        click.echo("")
        click.echo(f"Open {output_path / 'report.html'} in your browser to view the complete report.")
        click.echo("")
        
        if results.get("correlations"):
            click.echo("Top 10 Correlated Parameters:")
            for i, corr in enumerate(results["correlations"][:10], 1):
                click.echo(f"  {i:2d}. {corr['parameter']:25s} {corr['correlation']:.4f}")
            click.echo("")
        
    except Exception as e:
        click.echo(f"ERROR: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    """Main entry point for CLI."""
    cli()


if __name__ == "__main__":
    main()
