#!/usr/bin/env python3
"""Collect top N results from all optimization runs and aggregate them.

This script scans the optimization_results directory, extracts the top N
parameter combinations from each case, and creates aggregated summaries.
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd


def load_case_results(case_dir: Path, top_n: int = 100) -> Dict[str, Any] | None:
    """Load top N results from a single case.
    
    Args:
        case_dir: Directory containing optimization results for one case
        top_n: Number of top results to extract
        
    Returns:
        Dictionary with case results or None if files not found
    """
    summary_file = case_dir / "optimization_summary.csv"
    best_params_file = case_dir / "best_parameters.json"
    
    if not summary_file.exists():
        return None
    
    # Load summary CSV
    df = pd.read_csv(summary_file)
    
    # Sort by quality_score (descending) and take top N
    # Exclude baseline (combo_id == 0) from top N if requested
    non_baseline = df[df["combo_id"] != 0].copy()
    top_results = non_baseline.nlargest(top_n, "quality_score")
    
    # Load best parameters JSON if available
    best_params = None
    if best_params_file.exists():
        with open(best_params_file, "r") as f:
            best_params = json.load(f)
    
    return {
        "case_name": case_dir.name,
        "total_combinations": len(df),
        "baseline_quality_score": float(df[df["is_baseline"] == True]["quality_score"].iloc[0]) if len(df[df["is_baseline"] == True]) > 0 else None,
        "top_n": top_n,
        "top_results": top_results.to_dict("records"),
        "best_parameters": best_params,
    }


def aggregate_results(all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate results across all cases.
    
    Args:
        all_results: List of case result dictionaries
        
    Returns:
        Aggregated analysis dictionary
    """
    if not all_results:
        return {"error": "No results to aggregate"}
    
    # Combine all top results
    all_top_results = []
    for case_result in all_results:
        case_name = case_result["case_name"]
        for result in case_result["top_results"]:
            result["case_name"] = case_name
            all_top_results.append(result)
    
    # Create combined DataFrame
    combined_df = pd.DataFrame(all_top_results)
    
    # Calculate statistics
    stats = {
        "total_cases": len(all_results),
        "total_combinations": len(combined_df),
        "avg_quality_score": float(combined_df["quality_score"].mean()),
        "max_quality_score": float(combined_df["quality_score"].max()),
        "min_quality_score": float(combined_df["quality_score"].min()),
        "avg_signal_to_noise": float(combined_df["signal_to_noise"].replace([float("inf"), float("-inf")], None).mean()),
        "avg_alignment_rate": float(combined_df["alignment_rate"].mean()),
        "avg_high_quality_rate": float(combined_df["high_quality_rate"].mean()),
    }
    
    # Find overall best by quality score
    best_overall = combined_df.nlargest(1, "quality_score").iloc[0].to_dict()
    
    # Group by case to show best per case
    best_per_case = {}
    for case_result in all_results:
        case_name = case_result["case_name"]
        if case_result["top_results"]:
            best = max(case_result["top_results"], key=lambda x: x["quality_score"])
            best_per_case[case_name] = {
                "combo_id": int(best["combo_id"]),
                "quality_score": float(best["quality_score"]),
                "signal_to_noise": float(best["signal_to_noise"]) if best["signal_to_noise"] != "inf" else None,
                "alignment_rate": float(best["alignment_rate"]),
            }
    
    return {
        "statistics": stats,
        "best_overall": best_overall,
        "best_per_case": best_per_case,
        "all_results": all_results,
    }


def main() -> int:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Collect top N results from all optimization runs"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("optimization_results"),
        help="Base directory containing optimization results (default: optimization_results)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of top results to extract from each case (default: 100)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for aggregated results (default: results-dir/aggregated)",
    )
    parser.add_argument(
        "--min-quality",
        type=float,
        default=None,
        help="Minimum quality score threshold (optional)",
    )
    
    args = parser.parse_args()
    
    if not args.results_dir.exists():
        print(f"ERROR: Results directory not found: {args.results_dir}")
        return 1
    
    # Determine output directory
    if args.output_dir is None:
        args.output_dir = args.results_dir / "aggregated"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("Collecting Top Results from Optimization Runs")
    print("=" * 80)
    print(f"Results directory: {args.results_dir}")
    print(f"Top N per case: {args.top_n}")
    print(f"Output directory: {args.output_dir}")
    print("")
    
    # Find all case directories
    case_dirs = [
        d for d in args.results_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".") and d.name != "aggregated" and d.name != "job_scripts"
    ]
    
    if not case_dirs:
        print(f"ERROR: No case directories found in {args.results_dir}")
        return 1
    
    print(f"Found {len(case_dirs)} case directory(ies)")
    print("")
    
    # Load results from each case
    all_results = []
    for case_dir in sorted(case_dirs):
        print(f"Loading results from {case_dir.name}...")
        case_result = load_case_results(case_dir, top_n=args.top_n)
        
        if case_result is None:
            print(f"  ⚠️  Skipping {case_dir.name}: optimization_summary.csv not found")
            continue
        
        # Apply quality threshold if specified
        if args.min_quality is not None:
            filtered = [
                r for r in case_result["top_results"]
                if r["quality_score"] >= args.min_quality
            ]
            case_result["top_results"] = filtered
            case_result["top_n"] = len(filtered)
        
        all_results.append(case_result)
        print(f"  ✓ Loaded {len(case_result['top_results'])} top results")
        if case_result["baseline_quality_score"] is not None:
            print(f"    Baseline quality score: {case_result['baseline_quality_score']:.4f}")
        print("")
    
    if not all_results:
        print("ERROR: No valid results found")
        return 1
    
    # Aggregate results
    print("Aggregating results...")
    aggregated = aggregate_results(all_results)
    
    # Save aggregated results
    print("Saving aggregated results...")
    
    # Save full aggregated JSON
    aggregated_json = args.output_dir / "aggregated_top_results.json"
    with open(aggregated_json, "w") as f:
        json.dump(aggregated, f, indent=2, default=str)
    print(f"  ✓ Saved: {aggregated_json}")
    
    # Create combined CSV with top results from all cases
    all_top_records = []
    for case_result in all_results:
        for record in case_result["top_results"]:
            all_top_records.append(record)
    
    if all_top_records:
        combined_df = pd.DataFrame(all_top_records)
        # Sort by quality_score descending
        combined_df = combined_df.sort_values("quality_score", ascending=False)
        
        combined_csv = args.output_dir / "combined_top_results.csv"
        combined_df.to_csv(combined_csv, index=False)
        print(f"  ✓ Saved: {combined_csv}")
        
        # Save top 100 overall
        top_100_overall = combined_df.head(100)
        top_100_csv = args.output_dir / "top_100_overall.csv"
        top_100_overall.to_csv(top_100_csv, index=False)
        print(f"  ✓ Saved: {top_100_csv}")
    
    # Save summary statistics
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
    
    summary_json = args.output_dir / "summary_statistics.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  ✓ Saved: {summary_json}")
    
    # Print summary
    print("")
    print("=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Total cases processed: {aggregated['statistics']['total_cases']}")
    print(f"Total combinations: {aggregated['statistics']['total_combinations']}")
    print(f"Average quality score: {aggregated['statistics']['avg_quality_score']:.4f}")
    print(f"Best quality score: {aggregated['statistics']['max_quality_score']:.4f}")
    print(f"Best case: {aggregated['best_overall'].get('case_name', 'N/A')}")
    print(f"Best combo_id: {aggregated['best_overall'].get('combo_id', 'N/A')}")
    print("")
    print("Best per case:")
    for case_name, best in aggregated["best_per_case"].items():
        print(f"  {case_name}: quality_score={best['quality_score']:.4f}, "
              f"combo_id={best['combo_id']}")
    print("")
    print(f"Results saved to: {args.output_dir}")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

