"""Command-line interface for RNA-MAP optimization."""

import argparse
import sys
from pathlib import Path
from typing import Optional

from rna_map_optimization.optimizer import run_optimization
from rna_map_optimization.utils import generate_detailed_analysis, format_analysis_report


def find_fastq_files(case_dir: Path) -> tuple[Optional[Path], Optional[Path]]:
    """Find FASTQ files in a case directory.
    
    Args:
        case_dir: Directory containing test case files
        
    Returns:
        Tuple of (fastq1, fastq2) paths, or (None, None) if not found
    """
    fastq1 = None
    fastq2 = None
    
    # Common patterns for FASTQ files
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
    
    # Find input files
    fasta = find_fasta_file(case_dir)
    fastq1, fastq2 = find_fastq_files(case_dir)
    
    if not fasta:
        raise FileNotFoundError(f"No FASTA file found in {case_dir}")
    if not fastq1:
        raise FileNotFoundError(f"No FASTQ file found in {case_dir}")
    
    # Create output directory for this case
    output_dir = output_base / case_name
    
    # Use case name as study name if not provided
    if study_name is None:
        study_name = f"bowtie2_optimization_{case_name}"
    
    print(f"\n{'='*80}")
    print(f"Processing case: {case_name}")
    print(f"{'='*80}")
    print(f"FASTA: {fasta}")
    print(f"FASTQ1: {fastq1}")
    if fastq2:
        print(f"FASTQ2: {fastq2}")
    print(f"Output: {output_dir}")
    print()
    sys.stdout.flush()
    
    # Run optimization
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
    )
    
    # Print brief summary (detailed analysis will be shown in main)
    print(f"\nCase {case_name} complete:")
    print(f"  Best quality score: {results['best_value']:.4f}")
    print(f"  Best signal-to-noise: {results['best_trial'].user_attrs.get('signal_to_noise', 0.0):.2f}")
    print(f"  Alignment rate: {results['best_trial'].user_attrs.get('alignment_rate', 0.0):.2%}")
    print()
    
    return results


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Optimize Bowtie2 parameters using Optuna (Bayesian optimization)"
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--fasta",
        type=Path,
        help="Path to reference FASTA file (single case mode with explicit files)",
    )
    input_group.add_argument(
        "--case-dir",
        type=Path,
        help="Directory containing a single test case (auto-detects FASTA/FASTQ files)",
    )
    input_group.add_argument(
        "--cases-dir",
        type=Path,
        help="Directory containing multiple test cases (batch mode)",
    )
    
    parser.add_argument(
        "--fastq1",
        type=Path,
        help="Path to first FASTQ file (required in single case mode)",
    )
    parser.add_argument(
        "--fastq2",
        type=Path,
        default=None,
        help="Path to second FASTQ file (optional, for paired-end)",
    )
    
    # Output options
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("bowtie2_optimization_optuna"),
        help="Output directory for results",
    )
    
    # Optimization parameters
    parser.add_argument(
        "--n-trials",
        type=int,
        default=100,
        help="Number of Optuna trials to run (default: 100)",
    )
    parser.add_argument(
        "--mapq-cutoff",
        type=int,
        default=20,
        help="MAPQ cutoff for high-quality alignments (default: 20)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for alignment (or max if optimizing). Can be any positive integer.",
    )
    parser.add_argument(
        "--optimize-threads",
        action="store_true",
        help="Optimize number of threads (default: use fixed threads)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Timeout in seconds (optional)",
    )
    
    # Optuna options
    parser.add_argument(
        "--study-name",
        type=str,
        default=None,
        help="Name for Optuna study (default: auto-generated)",
    )
    parser.add_argument(
        "--storage",
        type=str,
        default=None,
        help="Optuna storage URL (e.g., sqlite:///study.db) for resuming studies",
    )
    
    args = parser.parse_args()
    
    # Single case directory mode (auto-detect files)
    if args.case_dir:
        if not args.case_dir.exists():
            parser.error(f"Case directory does not exist: {args.case_dir}")
        if not args.case_dir.is_dir():
            parser.error(f"Not a directory: {args.case_dir}")
        
        print("=" * 80)
        print("Bowtie2 Parameter Optimization with Optuna (Single Case)")
        print("=" * 80)
        print(f"Case directory: {args.case_dir}")
        print(f"Output directory: {args.output_dir}")
        print(f"Number of trials: {args.n_trials}")
        print(f"Threads: {args.threads} ({'optimizing' if args.optimize_threads else 'fixed'})")
        print()
        sys.stdout.flush()
        
        try:
            results = process_single_case(
                case_dir=args.case_dir,
                output_base=args.output_dir,
                n_trials=args.n_trials,
                mapq_cutoff=args.mapq_cutoff,
                threads=args.threads,
                optimize_threads=args.optimize_threads,
                timeout=args.timeout,
                study_name=args.study_name,
                storage=args.storage,
            )
            
            # Generate and display detailed analysis
            analysis = generate_detailed_analysis(results['study'], results['best_trial'])
            report = format_analysis_report(analysis, results['best_bt2_args'])
            print("\n" + report)
            
            output_path = args.output_dir / args.case_dir.name
            print(f"Results saved to: {output_path}")
            print()
        except Exception as e:
            print(f"ERROR: Failed to process case: {e}")
            return 1
    
    # Single case mode with explicit files
    elif args.fasta:
        if not args.fastq1:
            parser.error("--fastq1 is required when using --fasta")
        
        print("=" * 80)
        print("Bowtie2 Parameter Optimization with Optuna")
        print("=" * 80)
        print(f"FASTA: {args.fasta}")
        print(f"FASTQ1: {args.fastq1}")
        if args.fastq2:
            print(f"FASTQ2: {args.fastq2}")
        print(f"Output directory: {args.output_dir}")
        print(f"Number of trials: {args.n_trials}")
        print(f"Threads: {args.threads} ({'optimizing' if args.optimize_threads else 'fixed'})")
        print()
        sys.stdout.flush()
        
        results = run_optimization(
            fasta=args.fasta,
            fastq1=args.fastq1,
            fastq2=args.fastq2,
            output_dir=args.output_dir,
            n_trials=args.n_trials,
            mapq_cutoff=args.mapq_cutoff,
            threads=args.threads,
            optimize_threads=args.optimize_threads,
            timeout=args.timeout,
            study_name=args.study_name or "bowtie2_optimization",
            storage=args.storage,
        )
        
        # Generate and display detailed analysis
        analysis = generate_detailed_analysis(results['study'], results['best_trial'])
        report = format_analysis_report(analysis, results['best_bt2_args'])
        print("\n" + report)
        
        print(f"Results saved to: {args.output_dir}")
        print()
    
    # Batch mode (multiple cases)
    elif args.cases_dir:
        if not args.cases_dir.exists():
            parser.error(f"Cases directory does not exist: {args.cases_dir}")
        
        print("=" * 80)
        print("Batch Bowtie2 Parameter Optimization with Optuna")
        print("=" * 80)
        print(f"Cases directory: {args.cases_dir}")
        print(f"Output directory: {args.output_dir}")
        print(f"Number of trials per case: {args.n_trials}")
        print(f"Threads: {args.threads} ({'optimizing' if args.optimize_threads else 'fixed'})")
        print()
        sys.stdout.flush()
        
        # Find all case directories
        case_dirs = [d for d in args.cases_dir.iterdir() if d.is_dir()]
        
        if not case_dirs:
            print(f"No case directories found in {args.cases_dir}")
            return 1
        
        print(f"Found {len(case_dirs)} test cases:")
        for case_dir in case_dirs:
            print(f"  - {case_dir.name}")
        print()
        sys.stdout.flush()
        
        # Process each case
        all_results = {}
        for i, case_dir in enumerate(case_dirs, 1):
            try:
                print(f"\n[{i}/{len(case_dirs)}] Processing {case_dir.name}...")
                results = process_single_case(
                    case_dir=case_dir,
                    output_base=args.output_dir,
                    n_trials=args.n_trials,
                    mapq_cutoff=args.mapq_cutoff,
                    threads=args.threads,
                    optimize_threads=args.optimize_threads,
                    timeout=args.timeout,
                    study_name=None,  # Auto-generate from case name
                    storage=args.storage,
                )
                all_results[case_dir.name] = results
            except Exception as e:
                print(f"ERROR: Failed to process {case_dir.name}: {e}")
                continue
        
        # Summary
        print("\n" + "=" * 80)
        print("Batch Optimization Complete")
        print("=" * 80)
        print(f"Successfully processed {len(all_results)}/{len(case_dirs)} cases")
        print()
        
        if all_results:
            print("Summary of best results:")
            for case_name, results in all_results.items():
                best_snr = results['best_trial'].user_attrs.get("signal_to_noise", 0.0)
                best_align_rate = results['best_trial'].user_attrs.get("alignment_rate", 0.0)
                print(f"  {case_name}:")
                print(f"    Quality score: {results['best_value']:.4f}")
                print(f"    Signal-to-Noise: {best_snr:.2f}")
                print(f"    Alignment Rate: {best_align_rate:.2%}")
            print()
        
        print(f"All results saved to: {args.output_dir}")
        print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

