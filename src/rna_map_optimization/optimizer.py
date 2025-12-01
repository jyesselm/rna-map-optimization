"""Optuna-based Bowtie2 parameter optimization.

This module provides functions for optimizing Bowtie2 parameters using Optuna
Bayesian optimization to maximize signal-to-noise ratio in RNA-MAP alignments.
"""

import json
import pickle
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

try:
    import optuna
    from optuna.visualization import (
        plot_optimization_history,
        plot_param_importances,
        plot_parallel_coordinate,
    )
except ImportError:
    optuna = None
    plot_optimization_history = None
    plot_param_importances = None
    plot_parallel_coordinate = None

from rna_map_optimization.bowtie2 import (
    build_bowtie2_index,
    generate_bit_vectors_and_analyze,
    params_to_bowtie2_args,
    parse_bowtie2_stats,
    parse_sam_quality,
    run_bowtie2_alignment,
)
from rna_map_optimization.config_loader import get_parameter_ranges, load_config
from rna_map_optimization.utils import fasta_to_dict
from rna_map_mini.logger import get_logger

log = get_logger("BOWTIE2_OPTIMIZER")


# All functions moved to bowtie2 module


def create_optuna_objective(
    index: Path,
    fastq1: Path,
    fastq2: Optional[Path],
    ref_seqs: dict[str, str],
    output_dir: Path,
    mapq_cutoff: int,
    threads: int,
    optimize_threads: bool,
    param_ranges: Optional[Dict[str, Any]] = None,
):
    """Create Optuna objective function.

    Args:
        index: Path to Bowtie2 index
        fastq1: Path to first FASTQ file
        fastq2: Path to second FASTQ file (optional)
        ref_seqs: Dictionary of reference sequences
        output_dir: Output directory
        mapq_cutoff: MAPQ cutoff
        threads: Number of threads (if not optimizing)
        optimize_threads: Whether to optimize threads
        param_ranges: Optional parameter ranges from config file

    Returns:
        Objective function for Optuna
    """
    if optuna is None:
        raise ImportError("Optuna is not installed. Install with: pip install optuna")

    # Load parameter ranges from config (required)
    if param_ranges is None:
        config = load_config()
        param_ranges = get_parameter_ranges(config)

    if not param_ranges:
        raise ValueError(
            "No parameters found in config file. "
            "Please ensure config/optimization_config.yml contains parameter definitions."
        )

    def suggest_param(trial: optuna.Trial, param_name: str) -> Any:
        """Suggest parameter value from config.

        Args:
            trial: Optuna trial object
            param_name: Name of parameter

        Returns:
            Suggested parameter value

        Raises:
            KeyError: If parameter not found in config
        """
        if param_name not in param_ranges:
            raise KeyError(
                f"Parameter '{param_name}' not found in config file. "
                f"Please add it to config/optimization_config.yml"
            )

        param_config = param_ranges[param_name]
        param_type = param_config.get("type", "categorical")

        if param_type == "int":
            min_val = param_config["min"]
            max_val = param_config["max"]
            step = param_config.get("step", 1)
            return trial.suggest_int(param_name, min_val, max_val, step=step)
        elif param_type == "float":
            min_val = param_config["min"]
            max_val = param_config["max"]
            step = param_config.get("step", None)
            if step:
                return trial.suggest_float(param_name, min_val, max_val, step=step)
            return trial.suggest_float(param_name, min_val, max_val)
        elif param_type == "categorical":
            options = param_config["options"]
            return trial.suggest_categorical(param_name, options)
        else:
            raise ValueError(f"Unknown parameter type '{param_type}' for '{param_name}'")

    def objective(trial: optuna.Trial) -> float:
        """Optuna objective function to maximize.

        Args:
            trial: Optuna trial object

        Returns:
            Quality score to maximize
        """
        # Load all parameters directly from config
        seed_length = suggest_param(trial, "seed_length")
        maxins = suggest_param(trial, "maxins")
        seed_mismatches = suggest_param(trial, "seed_mismatches")
        score_min_choice = suggest_param(trial, "score_min")
        mp_choice = suggest_param(trial, "mismatch_penalty")
        gap_read_choice = suggest_param(trial, "gap_penalty_read")
        gap_ref_choice = suggest_param(trial, "gap_penalty_ref")
        sens_mode = suggest_param(trial, "sensitivity_mode")
        seed_interval_choice = suggest_param(trial, "seed_interval")
        np_penalty_choice = suggest_param(trial, "np_penalty")
        n_ceil_choice = suggest_param(trial, "n_ceil")
        gbar_choice = suggest_param(trial, "gbar")
        match_bonus_choice = suggest_param(trial, "match_bonus")
        extension_effort_choice = suggest_param(trial, "extension_effort")
        repetitive_effort_choice = suggest_param(trial, "repetitive_effort")

        # Minimum insert size (only for paired-end)
        minins_choice = None
        if fastq2 is not None and "minins" in param_ranges:
            minins_choice = suggest_param(trial, "minins")

        # Threads
        trial_threads = threads
        if optimize_threads:
            # Allow optimizing threads up to 2x the specified threads value
            # No hard upper limit - user can specify as many threads as they want
            trial_threads = trial.suggest_int("threads", 2, threads * 2)

        # Build parameter dictionary
        params = {
            "local": True,
            "no_unal": True,
            "no_discordant": True,
            "no_mixed": True,
            "seed_length": seed_length,
            "maxins": maxins,
        }

        if seed_mismatches > 0:
            params["seed_mismatches"] = seed_mismatches

        # Validate score_min and mismatch_penalty compatibility
        if score_min_choice and mp_choice:
            if score_min_choice.startswith("L,"):
                score_parts = score_min_choice.split(",")
                if len(score_parts) >= 3:
                    try:
                        if float(score_parts[1]) < 0 or float(score_parts[2]) < 0:
                            mp_parts = mp_choice.split(",")
                            if len(mp_parts) >= 2 and int(mp_parts[1]) > 0:
                                trial.set_user_attr("failure_reason", "incompatible_score_mismatch")
                                return 0.0
                    except (ValueError, IndexError):
                        pass

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

            # Early pruning
            if alignment_stats["total_reads"] == 0:
                trial.set_user_attr("failure_reason", "no_reads")
                return 0.0

            alignment_rate = alignment_stats.get("overall_alignment_rate", 0.0) / 100.0

            if alignment_rate < 0.3:
                trial.set_user_attr("failure_reason", "low_alignment_rate")
                trial.set_user_attr("alignment_rate", alignment_rate)
                return 0.0

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

            # Quality score
            quality_score = (
                0.40 * normalized_snr
                + 0.30 * alignment_rate
                + 0.20 * (avg_mapq / 60.0)
                + 0.10 * bv_acceptance_rate
            )

            # Store intermediate values
            trial.set_user_attr("signal_to_noise", snr)
            trial.set_user_attr("alignment_rate", alignment_rate)
            trial.set_user_attr("avg_mapq", avg_mapq)
            trial.set_user_attr("bv_accepted", bit_vector_metrics.get("accepted_bit_vectors", 0))
            trial.set_user_attr("total_reads", alignment_stats.get("total_reads", 0))
            trial.set_user_attr("total_alignments", quality_metrics.get("total_alignments", 0))
            trial.set_user_attr("total_bit_vectors", bit_vector_metrics.get("total_bit_vectors", 0))
            trial.set_user_attr(
                "bv_rejected_low_mapq", bit_vector_metrics.get("rejected_low_mapq", 0)
            )

            # Store per-construct information
            if bit_vector_metrics.get("constructs"):
                for construct_name, construct_stats in bit_vector_metrics["constructs"].items():
                    trial.set_user_attr(
                        f"construct_{construct_name}_aligned", construct_stats["aligned_reads"]
                    )
                    trial.set_user_attr(
                        f"construct_{construct_name}_snr", construct_stats["signal_to_noise"]
                    )
            if optimize_threads:
                trial.set_user_attr("threads_used", trial_threads)

            return quality_score

        except Exception as e:
            log.warning(f"Trial {trial.number} failed: {e}")
            trial.set_user_attr("failure_reason", f"exception: {str(e)[:50]}")
            return 0.0

    return objective


def run_optimization(
    fasta: Path,
    fastq1: Path,
    fastq2: Optional[Path],
    output_dir: Path,
    n_trials: int = 100,
    mapq_cutoff: int = 20,
    threads: int = 4,
    optimize_threads: bool = False,
    timeout: Optional[int] = None,
    study_name: str = "bowtie2_optimization",
    storage: Optional[str] = None,
    cleanup: bool = True,
) -> Dict:
    """Run Optuna optimization for Bowtie2 parameters.

    Args:
        fasta: Path to reference FASTA file
        fastq1: Path to first FASTQ file
        fastq2: Path to second FASTQ file (optional)
        output_dir: Output directory for results
        n_trials: Number of Optuna trials
        mapq_cutoff: MAPQ cutoff for high-quality alignments
        threads: Number of threads for alignment
        optimize_threads: Whether to optimize number of threads
        timeout: Timeout in seconds (optional)
        study_name: Name for Optuna study
        storage: Optuna storage URL (optional)

    Returns:
        Dictionary with optimization results
    """
    if optuna is None:
        raise ImportError("Optuna is not installed. Install with: pip install optuna")

    # Setup output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    index_dir = output_dir / "index"
    results_dir = output_dir / "results"
    results_dir.mkdir(exist_ok=True)

    # Read reference sequences
    ref_seqs = fasta_to_dict(fasta)

    # Build index
    log.info("Building Bowtie2 index...")
    index = build_bowtie2_index(fasta, index_dir)

    # Create Optuna study
    study = optuna.create_study(
        direction="maximize",
        study_name=study_name,
        storage=storage,
        load_if_exists=True,
    )

    # Load config for parameter ranges
    config = load_config()
    param_ranges = get_parameter_ranges(config)
    log.info("Loaded parameter ranges from config/optimization_config.yml")

    # Create objective function
    objective = create_optuna_objective(
        index,
        fastq1,
        fastq2,
        ref_seqs,
        results_dir,
        mapq_cutoff,
        threads,
        optimize_threads,
        param_ranges=param_ranges,
    )

    # Run optimization
    study.optimize(
        objective,
        n_trials=n_trials,
        timeout=timeout,
        show_progress_bar=True,
    )

    # Prepare results
    best_params = {
        "local": True,
        "no_unal": True,
        "no_discordant": True,
        "no_mixed": True,
        "seed_length": study.best_params["seed_length"],
        "maxins": study.best_params["maxins"],
    }

    # Add optional parameters
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
    if (
        "extension_effort" in study.best_params
        and study.best_params["extension_effort"] is not None
    ):
        best_params["extension_effort"] = study.best_params["extension_effort"]
    if (
        "repetitive_effort" in study.best_params
        and study.best_params["repetitive_effort"] is not None
    ):
        best_params["repetitive_effort"] = study.best_params["repetitive_effort"]
    if "minins" in study.best_params and study.best_params["minins"] is not None:
        best_params["minins"] = study.best_params["minins"]

    best_bt2_args = params_to_bowtie2_args(best_params)

    # Save results
    results_file = output_dir / "optuna_study.json"
    with open(results_file, "w") as f:
        json.dump(
            {
                "best_trial": {
                    "number": study.best_trial.number,
                    "value": study.best_value,
                    "params": study.best_params,
                    "user_attrs": study.best_trial.user_attrs,
                },
                "best_bowtie2_args": " ".join(best_bt2_args),
                "n_trials": len(study.trials),
            },
            f,
            indent=2,
        )

    # Save study
    study_file = output_dir / "optuna_study.pkl"
    with open(study_file, "wb") as f:
        pickle.dump(study, f)

    # Create summary DataFrame
    df_data = []
    for trial in study.trials:
        if trial.state == optuna.trial.TrialState.COMPLETE:
            df_data.append(
                {
                    "trial": trial.number,
                    "quality_score": trial.value,
                    "signal_to_noise": trial.user_attrs.get("signal_to_noise", 0.0),
                    "alignment_rate": trial.user_attrs.get("alignment_rate", 0.0),
                    "avg_mapq": trial.user_attrs.get("avg_mapq", 0.0),
                    "total_reads": trial.user_attrs.get("total_reads", 0),
                    "total_alignments": trial.user_attrs.get("total_alignments", 0),
                    "total_bit_vectors": trial.user_attrs.get("total_bit_vectors", 0),
                    "bv_accepted": trial.user_attrs.get("bv_accepted", 0),
                    "bv_rejected_low_mapq": trial.user_attrs.get("bv_rejected_low_mapq", 0),
                    **trial.params,
                }
            )

    if df_data:
        df = pd.DataFrame(df_data)
        summary_file = output_dir / "optuna_summary.csv"
        df.to_csv(summary_file, index=False)

    # Generate visualizations
    try:
        vis_dir = output_dir / "visualizations"
        vis_dir.mkdir(exist_ok=True)

        fig = plot_optimization_history(study)
        fig.write_html(str(vis_dir / "optimization_history.html"))

        try:
            fig = plot_param_importances(study)
            fig.write_html(str(vis_dir / "param_importances.html"))
        except Exception:
            pass

        try:
            fig = plot_parallel_coordinate(study)
            fig.write_html(str(vis_dir / "parallel_coordinate.html"))
        except Exception:
            pass
    except Exception as e:
        log.warning(f"Could not generate visualizations: {e}")

    # Run final detailed analysis on best trial
    best_trial_dir = results_dir / f"trial_{study.best_trial.number}"
    best_sam_file = best_trial_dir / "aligned.sam"

    if best_sam_file.exists():
        log.info("Running final bit vector analysis on best trial...")
        final_bit_vector_metrics = generate_bit_vectors_and_analyze(
            best_sam_file,
            ref_seqs,
            paired=fastq2 is not None,
            qscore_cutoff=25,
            mapq_cutoff=mapq_cutoff,
        )

        # Save detailed final metrics
        final_metrics_file = output_dir / "final_bit_vector_metrics.json"
        with open(final_metrics_file, "w") as f:
            json.dump(
                {
                    "trial_number": study.best_trial.number,
                    "bit_vector_metrics": final_bit_vector_metrics,
                    "best_params": best_params,
                    "best_bowtie2_args": " ".join(best_bt2_args),
                },
                f,
                indent=2,
                default=str,
            )

        log.info(f"Final analysis saved to {final_metrics_file}")

        # Update return dict with final metrics
        final_results = {
            "study": study,
            "best_params": best_params,
            "best_bt2_args": best_bt2_args,
            "best_value": study.best_value,
            "best_trial": study.best_trial,
            "final_bit_vector_metrics": final_bit_vector_metrics,
        }
    else:
        log.warning(f"Best trial SAM file not found: {best_sam_file}. Skipping final analysis.")
        final_results = {
            "study": study,
            "best_params": best_params,
            "best_bt2_args": best_bt2_args,
            "best_value": study.best_value,
            "best_trial": study.best_trial,
        }

    # Cleanup results directory to save space
    # All data is stored in optuna files (optuna_study.json, optuna_summary.csv, etc.)
    if cleanup:
        if best_sam_file.exists() and "final_bit_vector_metrics" in final_results:
            # Final analysis complete, safe to delete results directory
            try:
                shutil.rmtree(results_dir)
                log.info(f"Removed results directory: {results_dir}")
            except Exception as e:
                log.warning(f"Failed to remove results directory {results_dir}: {e}")
        else:
            log.info("Keeping results directory (final analysis not completed)")

    return final_results
