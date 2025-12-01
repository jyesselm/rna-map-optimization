"""Comprehensive analysis of all optimization runs.

This module provides functions to analyze all runs, identify features
that correlate with good/bad scores, and generate comprehensive reports.
"""

from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


def load_all_runs(csv_path: Path) -> pd.DataFrame:
    """Load all runs from optuna_summary.csv.
    
    Args:
        csv_path: Path to optuna_summary.csv
        
    Returns:
        DataFrame with all runs
    """
    df = pd.read_csv(csv_path)
    return df


def identify_metric_columns(df: pd.DataFrame) -> List[str]:
    """Identify metric columns vs parameter columns.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        List of metric column names
    """
    metric_cols = [
        "trial", "quality_score", "signal_to_noise", 
        "alignment_rate", "avg_mapq", "total_reads",
        "total_alignments", "total_bit_vectors", 
        "bv_accepted", "bv_rejected_low_mapq"
    ]
    return [col for col in metric_cols if col in df.columns]


def identify_parameter_columns(df: pd.DataFrame) -> List[str]:
    """Identify parameter columns.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        List of parameter column names
    """
    metric_cols = identify_metric_columns(df)
    param_cols = [col for col in df.columns if col not in metric_cols]
    return param_cols


def calculate_correlations(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate correlations between all features and quality_score.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        DataFrame with correlation values
    """
    valid_df = df[df["quality_score"] > 0].copy()
    
    if len(valid_df) == 0:
        return pd.DataFrame()
    
    param_cols = identify_parameter_columns(valid_df)
    correlations = {}
    
    for param in param_cols:
        if param not in valid_df.columns:
            continue
        
        param_series = valid_df[param].dropna()
        if len(param_series) == 0:
            continue
        
        if param_series.dtype == 'object' or param_series.dtype.name == 'category':
            param_encoded = pd.get_dummies(valid_df[param], prefix=param)
            if len(param_encoded.columns) > 0:
                corr = param_encoded.corrwith(valid_df["quality_score"]).abs()
                correlations[param] = float(corr.max()) if len(corr) > 0 else 0.0
        else:
            corr = valid_df[param].corr(valid_df["quality_score"])
            correlations[param] = abs(float(corr)) if not pd.isna(corr) else 0.0
    
    corr_df = pd.DataFrame(
        list(correlations.items()),
        columns=["parameter", "correlation"]
    ).sort_values("correlation", ascending=False)
    
    return corr_df


def compare_good_vs_bad_runs(
    df: pd.DataFrame, 
    top_percentile: float = 0.1,
    bottom_percentile: float = 0.1
) -> Dict[str, Any]:
    """Compare features in good vs bad runs.
    
    Args:
        df: DataFrame with all runs
        top_percentile: Top percentile to consider "good" (default: 0.1 = top 10%)
        bottom_percentile: Bottom percentile to consider "bad" (default: 0.1 = bottom 10%)
        
    Returns:
        Dictionary with comparison results
    """
    valid_df = df[df["quality_score"] > 0].copy()
    
    if len(valid_df) == 0:
        return {}
    
    top_threshold = valid_df["quality_score"].quantile(1 - top_percentile)
    bottom_threshold = valid_df["quality_score"].quantile(bottom_percentile)
    
    good_runs = valid_df[valid_df["quality_score"] >= top_threshold]
    bad_runs = valid_df[valid_df["quality_score"] <= bottom_threshold]
    
    param_cols = identify_parameter_columns(valid_df)
    comparison = {}
    
    for param in param_cols:
        if param not in valid_df.columns:
            continue
        
        good_vals = good_runs[param].dropna()
        bad_vals = bad_runs[param].dropna()
        
        if len(good_vals) == 0 or len(bad_vals) == 0:
            continue
        
        if good_vals.dtype == 'object' or good_vals.dtype.name == 'category':
            good_counts = good_vals.value_counts(normalize=True)
            bad_counts = bad_vals.value_counts(normalize=True)
            
            comparison[param] = {
                "type": "categorical",
                "good_distribution": good_counts.to_dict(),
                "bad_distribution": bad_counts.to_dict(),
            }
        else:
            comparison[param] = {
                "type": "numeric",
                "good_mean": float(good_vals.mean()),
                "bad_mean": float(bad_vals.mean()),
                "good_std": float(good_vals.std()),
                "bad_std": float(bad_vals.std()),
            }
    
    return {
        "top_threshold": float(top_threshold),
        "bottom_threshold": float(bottom_threshold),
        "n_good": len(good_runs),
        "n_bad": len(bad_runs),
        "comparisons": comparison,
    }


def create_correlation_heatmap(df: pd.DataFrame) -> go.Figure:
    """Create correlation heatmap between parameters and metrics.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        Plotly figure
    """
    valid_df = df[df["quality_score"] > 0].copy()
    
    if len(valid_df) == 0:
        return go.Figure()
    
    metric_cols = ["quality_score", "signal_to_noise", "alignment_rate", "avg_mapq"]
    metric_cols = [col for col in metric_cols if col in valid_df.columns]
    param_cols = identify_parameter_columns(valid_df)
    
    corr_matrix = []
    param_names = []
    
    for param in param_cols[:20]:
        if param not in valid_df.columns:
            continue
        
        param_series = valid_df[param].dropna()
        if len(param_series) == 0:
            continue
        
        row = []
        for metric in metric_cols:
            if param_series.dtype == 'object':
                param_encoded = pd.get_dummies(valid_df[param], prefix=param)
                if len(param_encoded.columns) > 0:
                    corr = param_encoded.corrwith(valid_df[metric]).abs().max()
                    row.append(float(corr) if not pd.isna(corr) else 0.0)
                else:
                    row.append(0.0)
            else:
                corr = valid_df[param].corr(valid_df[metric])
                row.append(abs(float(corr)) if not pd.isna(corr) else 0.0)
        
        corr_matrix.append(row)
        param_names.append(param)
    
    if not corr_matrix:
        return go.Figure()
    
    fig = go.Figure(data=go.Heatmap(
        z=corr_matrix,
        x=metric_cols,
        y=param_names,
        colorscale='Viridis',
        text=[[f"{val:.2f}" for val in row] for row in corr_matrix],
        texttemplate="%{text}",
        textfont={"size": 10},
    ))
    
    fig.update_layout(
        title="Parameter-Metric Correlation Heatmap",
        xaxis_title="Metrics",
        yaxis_title="Parameters",
        height=600,
    )
    
    return fig


def create_feature_distribution_plots(
    df: pd.DataFrame,
    comparison: Dict[str, Any]
) -> List[go.Figure]:
    """Create distribution plots for top parameters.
    
    Args:
        df: DataFrame with all runs
        comparison: Comparison results from compare_good_vs_bad_runs
        
    Returns:
        List of Plotly figures
    """
    valid_df = df[df["quality_score"] > 0].copy()
    if len(valid_df) == 0:
        return []
    
    corr_df = calculate_correlations(valid_df)
    top_params = corr_df.head(10)["parameter"].tolist()
    
    figures = []
    comparisons = comparison.get("comparisons", {})
    
    for param in top_params:
        if param not in valid_df.columns:
            continue
        
        param_series = valid_df[param].dropna()
        if len(param_series) == 0:
            continue
        
        if param_series.dtype == 'object':
            param_comp = comparisons.get(param, {})
            good_dist = param_comp.get("good_distribution", {})
            bad_dist = param_comp.get("bad_distribution", {})
            
            all_values = set(list(good_dist.keys()) + list(bad_dist.keys()))
            
            good_vals = [good_dist.get(v, 0) for v in all_values]
            bad_vals = [bad_dist.get(v, 0) for v in all_values]
            
            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=list(all_values),
                y=good_vals,
                name="Top 10% Runs",
                marker_color='green',
                opacity=0.7,
            ))
            fig.add_trace(go.Bar(
                x=list(all_values),
                y=bad_vals,
                name="Bottom 10% Runs",
                marker_color='red',
                opacity=0.7,
            ))
            
            fig.update_layout(
                title=f"Distribution: {param}",
                xaxis_title=param,
                yaxis_title="Frequency",
                barmode='group',
                height=400,
            )
        else:
            param_comp = comparisons.get(param, {})
            good_mean = param_comp.get("good_mean", 0)
            bad_mean = param_comp.get("bad_mean", 0)
            
            fig = go.Figure()
            fig.add_trace(go.Histogram(
                x=valid_df[param],
                name="All Runs",
                opacity=0.5,
                nbinsx=20,
            ))
            fig.add_vline(
                x=good_mean,
                line_dash="dash",
                line_color="green",
                annotation_text=f"Top 10% mean: {good_mean:.2f}",
            )
            fig.add_vline(
                x=bad_mean,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Bottom 10% mean: {bad_mean:.2f}",
            )
            
            fig.update_layout(
                title=f"Distribution: {param}",
                xaxis_title=param,
                yaxis_title="Count",
                height=400,
            )
        
        figures.append(fig)
    
    return figures


def create_scatter_plots(df: pd.DataFrame) -> List[go.Figure]:
    """Create scatter plots of top parameters vs quality_score.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        List of Plotly figures
    """
    valid_df = df[df["quality_score"] > 0].copy()
    if len(valid_df) == 0:
        return []
    
    corr_df = calculate_correlations(valid_df)
    top_params = corr_df.head(6)["parameter"].tolist()
    
    figures = []
    
    for param in top_params:
        if param not in valid_df.columns:
            continue
        
        param_series = valid_df[param].dropna()
        if len(param_series) == 0:
            continue
        
        if param_series.dtype == 'object':
            param_encoded = pd.get_dummies(valid_df[param], prefix=param)
            if len(param_encoded.columns) > 0:
                for col in param_encoded.columns[:3]:
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=param_encoded[col],
                        y=valid_df["quality_score"],
                        mode='markers',
                        name=col,
                        marker=dict(size=5, opacity=0.6),
                    ))
                    fig.update_layout(
                        title=f"{param} ({col}) vs Quality Score",
                        xaxis_title=col,
                        yaxis_title="Quality Score",
                        height=400,
                    )
                    figures.append(fig)
        else:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=valid_df[param],
                y=valid_df["quality_score"],
                mode='markers',
                marker=dict(size=5, opacity=0.6),
            ))
            fig.update_layout(
                title=f"{param} vs Quality Score",
                xaxis_title=param,
                yaxis_title="Quality Score",
                height=400,
            )
            figures.append(fig)
    
    return figures


def get_best_parameters(df: pd.DataFrame) -> Dict[str, Any]:
    """Extract best parameters from the run with highest quality score.
    
    Args:
        df: DataFrame with all runs
        
    Returns:
        Dictionary with best parameters and metrics
    """
    valid_df = df[df["quality_score"] > 0].copy()
    if len(valid_df) == 0:
        return {}
    
    best_run = valid_df.loc[valid_df["quality_score"].idxmax()]
    param_cols = identify_parameter_columns(valid_df)
    
    best_params = {}
    for param in param_cols:
        if param in best_run and pd.notna(best_run[param]):
            best_params[param] = best_run[param]
    
    return {
        "trial": int(best_run.get("trial", 0)),
        "quality_score": float(best_run["quality_score"]),
        "signal_to_noise": float(best_run.get("signal_to_noise", 0.0)),
        "alignment_rate": float(best_run.get("alignment_rate", 0.0)),
        "avg_mapq": float(best_run.get("avg_mapq", 0.0)),
        "total_reads": int(best_run.get("total_reads", 0)),
        "total_alignments": int(best_run.get("total_alignments", 0)),
        "total_bit_vectors": int(best_run.get("total_bit_vectors", 0)),
        "bv_accepted": int(best_run.get("bv_accepted", 0)),
        "bv_rejected_low_mapq": int(best_run.get("bv_rejected_low_mapq", 0)),
        "parameters": best_params,
    }


def get_parameter_consensus(df: pd.DataFrame, top_percentile: float = 0.1) -> Dict[str, Any]:
    """Get consensus parameter values from top runs.
    
    Args:
        df: DataFrame with all runs
        top_percentile: Top percentile to consider (default: 0.1 = top 10%)
        
    Returns:
        Dictionary with consensus information for each parameter
    """
    valid_df = df[df["quality_score"] > 0].copy()
    if len(valid_df) == 0:
        return {}
    
    threshold = valid_df["quality_score"].quantile(1 - top_percentile)
    top_runs = valid_df[valid_df["quality_score"] >= threshold]
    
    param_cols = identify_parameter_columns(valid_df)
    consensus = {}
    
    for param in param_cols:
        if param not in top_runs.columns:
            continue
        
        param_series = top_runs[param].dropna()
        if len(param_series) == 0:
            continue
        
        if param_series.dtype == 'object' or param_series.dtype.name == 'category':
            value_counts = param_series.value_counts()
            total = len(param_series)
            consensus[param] = {
                "type": "categorical",
                "most_common": value_counts.index[0] if len(value_counts) > 0 else None,
                "frequency": int(value_counts.iloc[0]) if len(value_counts) > 0 else 0,
                "frequency_pct": float((value_counts.iloc[0] / total) * 100) if len(value_counts) > 0 else 0.0,
                "all_values": value_counts.to_dict(),
            }
        else:
            consensus[param] = {
                "type": "numeric",
                "mean": float(param_series.mean()),
                "median": float(param_series.median()),
                "std": float(param_series.std()),
                "min": float(param_series.min()),
                "max": float(param_series.max()),
            }
    
    return {
        "threshold": float(threshold),
        "n_runs": len(top_runs),
        "consensus": consensus,
    }


def get_search_expansion_suggestions(
    df: pd.DataFrame,
    corr_df: pd.DataFrame,
    consensus: Dict[str, Any],
) -> Dict[str, List[str]]:
    """Generate suggestions for where to expand parameter search.
    
    Args:
        df: DataFrame with all runs
        corr_df: DataFrame with parameter correlations
        consensus: Consensus results from get_parameter_consensus
        
    Returns:
        Dictionary with suggestions categorized
    """
    valid_df = df[df["quality_score"] > 0].copy()
    if len(valid_df) == 0:
        return {}
    
    suggestions = {
        "expand_range": [],
        "explore_more": [],
        "reduce_range": [],
        "stable": [],
    }
    
    param_cols = identify_parameter_columns(valid_df)
    consensus_data = consensus.get("consensus", {})
    
    for _, row in corr_df.iterrows():
        param = row["parameter"]
        correlation = row["correlation"]
        
        if param not in consensus_data:
            continue
        
        param_info = consensus_data[param]
        param_series = valid_df[param].dropna()
        
        if len(param_series) == 0:
            continue
        
        if param_info["type"] == "numeric":
            mean_val = param_info["mean"]
            min_val = param_info["min"]
            max_val = param_info["max"]
            std_val = param_info["std"]
            
            # Check if values are at edges of current range
            actual_min = float(param_series.min())
            actual_max = float(param_series.max())
            
            # High correlation + low diversity + at edge = expand range
            if correlation > 0.2 and std_val < mean_val * 0.1:
                if abs(mean_val - actual_min) < std_val or abs(mean_val - actual_max) < std_val:
                    suggestions["expand_range"].append({
                        "parameter": param,
                        "reason": f"High correlation ({correlation:.3f}), low diversity, values near edge",
                        "current_range": f"{actual_min:.2f} - {actual_max:.2f}",
                        "consensus_mean": f"{mean_val:.2f}",
                    })
            
            # High correlation + low diversity = explore more
            if correlation > 0.2 and std_val < mean_val * 0.15:
                suggestions["explore_more"].append({
                    "parameter": param,
                    "reason": f"High correlation ({correlation:.3f}) but low diversity (std={std_val:.2f})",
                    "suggestion": "Try more values around the consensus mean",
                })
            
            # Low correlation + high diversity = can reduce
            if correlation < 0.1 and std_val > mean_val * 0.3:
                suggestions["reduce_range"].append({
                    "parameter": param,
                    "reason": f"Low correlation ({correlation:.3f}) but high diversity",
                    "suggestion": "Consider narrowing search range",
                })
            
            # Very stable (low std) = stable parameter
            if std_val < mean_val * 0.05:
                suggestions["stable"].append({
                    "parameter": param,
                    "value": f"{mean_val:.2f}",
                    "std": f"{std_val:.2f}",
                })
        
        else:
            # Categorical parameter
            freq_pct = param_info["frequency_pct"]
            
            # High correlation + high consensus = explore other values
            if correlation > 0.2 and freq_pct > 80:
                suggestions["explore_more"].append({
                    "parameter": param,
                    "reason": f"High correlation ({correlation:.3f}) but one value dominates ({freq_pct:.1f}%)",
                    "suggestion": "Try other categorical options",
                })
            
            # Very high consensus = stable
            if freq_pct > 90:
                suggestions["stable"].append({
                    "parameter": param,
                    "value": param_info["most_common"],
                    "frequency_pct": f"{freq_pct:.1f}%",
                })
    
    return suggestions


def generate_html_report(
    df: pd.DataFrame,
    output_path: Path,
    comparison: Dict[str, Any],
    heatmap_fig: go.Figure,
    dist_figs: List[go.Figure],
    scatter_figs: List[go.Figure],
) -> None:
    """Generate comprehensive HTML report with all plots embedded.
    
    Args:
        df: DataFrame with all runs
        output_path: Path to save HTML report
        comparison: Comparison results from compare_good_vs_bad_runs
        heatmap_fig: Correlation heatmap figure
        dist_figs: List of distribution plot figures
        scatter_figs: List of scatter plot figures
    """
    valid_df = df[df["quality_score"] > 0].copy()
    
    html_parts = []
    html_parts.append("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Comprehensive Optimization Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1 { color: #2c3e50; }
            h2 { color: #34495e; margin-top: 30px; }
            h3 { color: #7f8c8d; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #3498db; color: white; }
            tr:nth-child(even) { background-color: #f2f2f2; }
            .summary-box { background-color: #ecf0f1; padding: 15px; margin: 20px 0; border-radius: 5px; }
            .good { color: #27ae60; font-weight: bold; }
            .bad { color: #e74c3c; font-weight: bold; }
            .plot-container { margin: 30px 0; padding: 20px; background-color: #f8f9fa; border-radius: 5px; }
            .plot-div { margin: 20px 0; }
        </style>
    </head>
    <body>
    """)
    
    html_parts.append("<h1>Comprehensive Optimization Analysis Report</h1>")
    
    html_parts.append("<div class='summary-box'>")
    html_parts.append(f"<h2>Summary Statistics</h2>")
    html_parts.append(f"<p><strong>Total Runs:</strong> {len(df)}</p>")
    html_parts.append(f"<p><strong>Valid Runs (quality_score > 0):</strong> {len(valid_df)}</p>")
    html_parts.append(f"<p><strong>Failed Runs:</strong> {len(df) - len(valid_df)}</p>")
    
    if len(valid_df) > 0:
        html_parts.append(f"<p><strong>Mean Quality Score:</strong> {valid_df['quality_score'].mean():.4f}</p>")
        html_parts.append(f"<p><strong>Max Quality Score:</strong> {valid_df['quality_score'].max():.4f}</p>")
        html_parts.append(f"<p><strong>Min Quality Score:</strong> {valid_df['quality_score'].min():.4f}</p>")
        html_parts.append(f"<p><strong>Std Dev Quality Score:</strong> {valid_df['quality_score'].std():.4f}</p>")
    html_parts.append("</div>")
    
    if len(valid_df) > 0:
        corr_df = calculate_correlations(valid_df)
        best_params = get_best_parameters(df)
        consensus = get_parameter_consensus(df)
        suggestions = get_search_expansion_suggestions(df, corr_df, consensus)
        
        # Best Parameters Section
        html_parts.append("<h2>Best Parameters</h2>")
        html_parts.append("<div class='summary-box'>")
        html_parts.append(f"<p><strong>Trial:</strong> {best_params.get('trial', 'N/A')}</p>")
        html_parts.append(f"<p><strong>Quality Score:</strong> {best_params.get('quality_score', 0):.4f}</p>")
        html_parts.append(f"<p><strong>Signal-to-Noise:</strong> {best_params.get('signal_to_noise', 0):.2f}</p>")
        html_parts.append(f"<p><strong>Alignment Rate:</strong> {best_params.get('alignment_rate', 0):.2%}</p>")
        html_parts.append(f"<p><strong>Average MAPQ:</strong> {best_params.get('avg_mapq', 0):.1f}</p>")
        
        # Read/Alignment/Bit Vector Pipeline
        total_reads = best_params.get('total_reads', 0)
        total_alignments = best_params.get('total_alignments', 0)
        total_bit_vectors = best_params.get('total_bit_vectors', 0)
        bv_accepted = best_params.get('bv_accepted', 0)
        bv_rejected = best_params.get('bv_rejected_low_mapq', 0)
        
        if total_reads > 0:
            html_parts.append("<h3>Read Pipeline Statistics</h3>")
            html_parts.append(f"<p><strong>Total Reads:</strong> {total_reads:,}</p>")
            html_parts.append(f"<p><strong>Total Alignments:</strong> {total_alignments:,} ({total_alignments/total_reads*100:.1f}% of reads)</p>")
            html_parts.append(f"<p><strong>Total Bit Vectors:</strong> {total_bit_vectors:,} ({total_bit_vectors/total_reads*100:.1f}% of reads)</p>")
            html_parts.append(f"<p><strong>Bit Vectors Accepted:</strong> {bv_accepted:,} ({bv_accepted/total_reads*100:.1f}% of reads, {bv_accepted/total_bit_vectors*100:.1f}% of bit vectors)</p>")
            if bv_rejected > 0:
                html_parts.append(f"<p><strong>Bit Vectors Rejected (Low MAPQ):</strong> {bv_rejected:,} ({bv_rejected/total_bit_vectors*100:.1f}% of bit vectors)</p>")
        
        html_parts.append("</div>")
        
        html_parts.append("<h3>Best Parameter Values</h3>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Parameter</th><th>Value</th></tr>")
        for param, value in sorted(best_params.get("parameters", {}).items()):
            html_parts.append(f"<tr><td>{param}</td><td>{value}</td></tr>")
        html_parts.append("</table>")
        
        # Parameter Consensus Section
        html_parts.append("<h2>Parameter Consensus (Top 10% Runs)</h2>")
        html_parts.append(f"<p><strong>Consensus Threshold:</strong> {consensus.get('threshold', 0):.4f}</p>")
        html_parts.append(f"<p><strong>Number of Runs Analyzed:</strong> {consensus.get('n_runs', 0)}</p>")
        
        html_parts.append("<h3>Consensus Values</h3>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Parameter</th><th>Type</th><th>Consensus Value</th><th>Details</th></tr>")
        
        consensus_data = consensus.get("consensus", {})
        for param, info in sorted(consensus_data.items()):
            if info["type"] == "numeric":
                html_parts.append(
                    f"<tr><td>{param}</td><td>Numeric</td>"
                    f"<td>{info['mean']:.2f} (median: {info['median']:.2f})</td>"
                    f"<td>Range: {info['min']:.2f} - {info['max']:.2f}, Std: {info['std']:.2f}</td></tr>"
                )
            else:
                html_parts.append(
                    f"<tr><td>{param}</td><td>Categorical</td>"
                    f"<td class='good'>{info['most_common']}</td>"
                    f"<td>{info['frequency_pct']:.1f}% of top runs</td></tr>"
                )
        html_parts.append("</table>")
        
        # Search Expansion Suggestions
        html_parts.append("<h2>Search Expansion Suggestions</h2>")
        
        if suggestions.get("expand_range"):
            html_parts.append("<h3>Parameters to Expand Range</h3>")
            html_parts.append("<p>These parameters have high correlation and values near the edge of the search space:</p>")
            html_parts.append("<table>")
            html_parts.append("<tr><th>Parameter</th><th>Reason</th><th>Current Range</th><th>Consensus Mean</th></tr>")
            for item in suggestions["expand_range"]:
                html_parts.append(
                    f"<tr><td>{item['parameter']}</td><td>{item['reason']}</td>"
                    f"<td>{item['current_range']}</td><td>{item['consensus_mean']}</td></tr>"
                )
            html_parts.append("</table>")
        
        if suggestions.get("explore_more"):
            html_parts.append("<h3>Parameters to Explore More</h3>")
            html_parts.append("<p>These parameters show promise but need more exploration:</p>")
            html_parts.append("<table>")
            html_parts.append("<tr><th>Parameter</th><th>Reason</th><th>Suggestion</th></tr>")
            for item in suggestions["explore_more"]:
                html_parts.append(
                    f"<tr><td>{item['parameter']}</td><td>{item['reason']}</td>"
                    f"<td>{item.get('suggestion', 'N/A')}</td></tr>"
                )
            html_parts.append("</table>")
        
        if suggestions.get("reduce_range"):
            html_parts.append("<h3>Parameters to Reduce Range</h3>")
            html_parts.append("<p>These parameters have low correlation - consider narrowing search:</p>")
            html_parts.append("<table>")
            html_parts.append("<tr><th>Parameter</th><th>Reason</th><th>Suggestion</th></tr>")
            for item in suggestions["reduce_range"]:
                html_parts.append(
                    f"<tr><td>{item['parameter']}</td><td>{item['reason']}</td>"
                    f"<td>{item.get('suggestion', 'N/A')}</td></tr>"
                )
            html_parts.append("</table>")
        
        if suggestions.get("stable"):
            html_parts.append("<h3>Stable Parameters</h3>")
            html_parts.append("<p>These parameters show consistent values in top runs - consider fixing them:</p>")
            html_parts.append("<table>")
            html_parts.append("<tr><th>Parameter</th><th>Consensus Value</th><th>Details</th></tr>")
            for item in suggestions["stable"]:
                html_parts.append(
                    f"<tr><td>{item['parameter']}</td><td>{item['value']}</td>"
                    f"<td>{item.get('std', item.get('frequency_pct', 'N/A'))}</td></tr>"
                )
            html_parts.append("</table>")
        
        html_parts.append("<h2>Top Correlated Parameters</h2>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Parameter</th><th>Correlation with Quality Score</th></tr>")
        for _, row in corr_df.head(20).iterrows():
            html_parts.append(f"<tr><td>{row['parameter']}</td><td>{row['correlation']:.4f}</td></tr>")
        html_parts.append("</table>")
        
        html_parts.append("<h2>Good vs Bad Runs Comparison</h2>")
        html_parts.append(f"<p><strong>Top 10% Threshold:</strong> {comparison.get('top_threshold', 0):.4f}</p>")
        html_parts.append(f"<p><strong>Bottom 10% Threshold:</strong> {comparison.get('bottom_threshold', 0):.4f}</p>")
        html_parts.append(f"<p><strong>Number of Good Runs:</strong> {comparison.get('n_good', 0)}</p>")
        html_parts.append(f"<p><strong>Number of Bad Runs:</strong> {comparison.get('n_bad', 0)}</p>")
        
        html_parts.append("<h3>Parameter Differences (Good vs Bad)</h3>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Parameter</th><th>Type</th><th>Good Value</th><th>Bad Value</th></tr>")
        
        comparisons = comparison.get("comparisons", {})
        for param, comp_data in list(comparisons.items())[:20]:
            param_type = comp_data.get("type", "unknown")
            if param_type == "numeric":
                good_val = comp_data.get("good_mean", 0)
                bad_val = comp_data.get("bad_mean", 0)
                html_parts.append(
                    f"<tr><td>{param}</td><td>Numeric</td>"
                    f"<td class='good'>{good_val:.4f}</td><td class='bad'>{bad_val:.4f}</td></tr>"
                )
            else:
                good_dist = comp_data.get("good_distribution", {})
                bad_dist = comp_data.get("bad_distribution", {})
                good_top = max(good_dist.items(), key=lambda x: x[1])[0] if good_dist else "N/A"
                bad_top = max(bad_dist.items(), key=lambda x: x[1])[0] if bad_dist else "N/A"
                html_parts.append(
                    f"<tr><td>{param}</td><td>Categorical</td>"
                    f"<td class='good'>{good_top}</td><td class='bad'>{bad_top}</td></tr>"
                )
        html_parts.append("</table>")
    
    html_parts.append("<h2>Visualizations</h2>")
    
    plotly_js_included = False
    
    if heatmap_fig and heatmap_fig.data:
        html_parts.append("<div class='plot-container'>")
        html_parts.append("<h3>Correlation Heatmap</h3>")
        html_parts.append("<div class='plot-div'>")
        html_str = heatmap_fig.to_html(include_plotlyjs='cdn', div_id="heatmap", full_html=False)
        html_parts.append(html_str)
        html_parts.append("</div>")
        html_parts.append("</div>")
        plotly_js_included = True
    
    if dist_figs:
        html_parts.append("<div class='plot-container'>")
        html_parts.append("<h3>Distribution Plots</h3>")
        html_parts.append("<p>Distribution of top parameters in good vs bad runs:</p>")
        for i, fig in enumerate(dist_figs):
            if fig and fig.data:
                html_parts.append(f"<div class='plot-div'>")
                html_str = fig.to_html(
                    include_plotlyjs='cdn' if not plotly_js_included else False, 
                    div_id=f"dist_{i+1}",
                    full_html=False
                )
                html_parts.append(html_str)
                html_parts.append("</div>")
                if not plotly_js_included:
                    plotly_js_included = True
        html_parts.append("</div>")
    
    if scatter_figs:
        html_parts.append("<div class='plot-container'>")
        html_parts.append("<h3>Scatter Plots</h3>")
        html_parts.append("<p>Scatter plots showing relationship between parameters and quality score:</p>")
        for i, fig in enumerate(scatter_figs):
            if fig and fig.data:
                html_parts.append(f"<div class='plot-div'>")
                html_str = fig.to_html(
                    include_plotlyjs='cdn' if not plotly_js_included else False,
                    div_id=f"scatter_{i+1}",
                    full_html=False
                )
                html_parts.append(html_str)
                html_parts.append("</div>")
                if not plotly_js_included:
                    plotly_js_included = True
        html_parts.append("</div>")
    
    html_parts.append("</body></html>")
    
    with open(output_path, "w") as f:
        f.write("\n".join(html_parts))


def run_comprehensive_analysis(
    results_dir: Path,
    output_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run comprehensive analysis on all runs.
    
    Args:
        results_dir: Directory containing optuna_summary.csv
        output_dir: Output directory for plots and report (default: results_dir/analysis)
        
    Returns:
        Dictionary with analysis results
    """
    if output_dir is None:
        output_dir = results_dir / "comprehensive_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = results_dir / "optuna_summary.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"optuna_summary.csv not found in {results_dir}")
    
    df = load_all_runs(csv_path)
    valid_df = df[df["quality_score"] > 0].copy()
    
    if len(valid_df) == 0:
        return {"error": "No valid runs found"}
    
    corr_df = calculate_correlations(valid_df)
    comparison = compare_good_vs_bad_runs(df)
    
    corr_df.to_csv(output_dir / "correlations.csv", index=False)
    
    with open(output_dir / "comparison.json", "w") as f:
        json.dump(comparison, f, indent=2, default=str)
    
    heatmap_fig = create_correlation_heatmap(df)
    dist_figs = create_feature_distribution_plots(df, comparison)
    scatter_figs = create_scatter_plots(df)
    
    generate_html_report(
        df, 
        output_dir / "report.html", 
        comparison, 
        heatmap_fig,
        dist_figs,
        scatter_figs,
    )
    
    return {
        "correlations": corr_df.to_dict("records"),
        "comparison": comparison,
        "output_dir": str(output_dir),
    }

