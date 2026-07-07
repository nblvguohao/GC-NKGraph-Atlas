"""
GC-NKGraph-Atlas SST Axis Analysis (Phase 14R).

Computes the Serine-Sphingomyelin-Topology transcriptional axis scores
from scRNA-seq data (NK subset).

Gene modules from Zheng et al. Nat Immunol 2023:
  tumor_serine_capacity → nk_sm_balance → nk_protrusion_machinery → cytotoxicity

Usage:
    python src/topology/sst_axis.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# --- Shared SST config (single source of truth) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.sst_config import load_sst_modules, get_sst_genes, SST_MODULES  # noqa: E402

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Aliases for backward compatibility — prefer load_sst_modules() in new code
# ---------------------------------------------------------------------------
MODULES: Dict[str, Dict] = SST_MODULES


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def mean_zscore(df: pd.DataFrame, genes: List[str]) -> pd.Series:
    """Compute mean z-score of a gene set across cells.

    Args:
        df: DataFrame of shape (cells × genes).
        genes: List of gene symbols to average over.

    Returns:
        Series of per-cell mean z-scores. Genes not in df are silently skipped;
        if no genes are found, returns zeros.
    """
    available = [g for g in genes if g in df.columns]
    if not available:
        return pd.Series(0.0, index=df.index, dtype=float)
    z = (df[available] - df[available].mean(0)) / df[available].std(0, ddof=0)
    return z.fillna(0).mean(axis=1)


def compute_sst_scores(
    adata,  # AnnData
    modules: Optional[Dict[str, Dict]] = None,
) -> pd.DataFrame:
    """Compute all SST-axis scores from an AnnData object.

    Args:
        adata: AnnData with cells × genes expression matrix.
        modules: Dict of gene modules (default: loaded from sst_axis_config.yaml).

    Returns:
        DataFrame with per-cell module scores + derived scores.
    """
    if modules is None:
        modules = load_sst_modules()

    expr = adata.to_df()  # cells × genes
    scores = pd.DataFrame(index=expr.index)

    # Per-module scores
    for name, mod in modules.items():
        gene_list = mod["genes"]
        scores[f"{name}_score"] = mean_zscore(expr, gene_list)
        n_found = sum(1 for g in gene_list if g in expr.columns)
        log(f"  {name}: {n_found}/{len(gene_list)} genes")

    # Derived scores
    scores["nk_sm_balance_score"] = (
        scores.get("nk_sm_synthesis_score", 0)
        - scores.get("nk_sm_catabolism_score", 0)
    )
    scores["nk_topology_permissive_score"] = (
        scores.get("nk_sm_balance_score", 0)
        + scores.get("nk_protrusion_machinery_score", 0)
    ) / 2

    # Integrated SST axis score (sign of tumor term calibrated on liver control)
    scores["sst_axis_score"] = (
        scores.get("tumor_serine_capacity_score", 0)
        + scores.get("nk_topology_permissive_score", 0)
        + scores.get("nk_synapse_cytotoxicity_outcome_score", 0)
    )

    # Metadata from AnnData .obs
    scores["sample_id"] = adata.obs.get("sample_id", "unknown")
    scores["tissue"] = adata.obs.get("tissue", "unknown")
    scores["condition"] = adata.obs.get("condition", "unknown")

    return scores


def compute_sample_summary(
    scores: pd.DataFrame,
    group_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Aggregate per-cell SST scores to per-sample summary statistics."""
    if group_cols is None:
        group_cols = ["sample_id", "tissue", "condition"]

    metric_cols = [
        c for c in scores.columns
        if c.endswith("_score") and c not in ("checkpoint_link_score",)
        and c in scores.columns
    ]
    # Ensure all group_cols are present
    actual_groups = [c for c in group_cols if c in scores.columns]

    if not actual_groups:
        return scores[metric_cols].agg(["mean", "std"]).T

    return scores.groupby(actual_groups)[metric_cols].agg(["mean", "std"]).round(4)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    log("=" * 50)
    log("SST Axis Analysis (Phase 14R)")
    log("=" * 50)

    # Config
    nk_path = "data/processed/scrna/gc_nk_subset.h5ad"
    out_dir = "results/tables"
    fig_dir = "results/figures"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # Load NK subset
    log(f"Loading NK subset: {nk_path}")
    import scanpy as sc
    adata = sc.read(nk_path)
    log(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    log(f"  Samples: {adata.obs['sample_id'].value_counts().to_dict()}")

    # Compute SST scores
    log("Computing SST scores...")
    scores = compute_sst_scores(adata)
    log(f"  Scores shape: {scores.shape}")

    # Save per-cell scores
    scores.to_csv(
        os.path.join(out_dir, "sst_axis_scores_single_cell.tsv"),
        sep="\t",
        index_label="cell_id",
    )
    log("  Saved: sst_axis_scores_single_cell.tsv")

    # Per-sample summary
    log("\nPer-sample SST axis summary:")
    summary = compute_sample_summary(scores)
    print(summary.to_string())
    summary.to_csv(os.path.join(out_dir, "sst_axis_sample_summary.tsv"), sep="\t")

    # Compare conditions
    log("\nCondition comparison (NK scores):")
    cond_metrics = [
        "nk_sm_balance_score",
        "nk_protrusion_machinery_score",
        "nk_topology_permissive_score",
        "nk_synapse_cytotoxicity_outcome_score",
        "checkpoint_link_score",
    ]
    available_metrics = [m for m in cond_metrics if m in scores.columns]
    cond_summary = scores.groupby("condition")[available_metrics].agg(["mean", "std"]).round(4)
    print(cond_summary.to_string())
    cond_summary.to_csv(os.path.join(out_dir, "sst_axis_condition_comparison.tsv"), sep="\t")

    # Plot
    log("Generating figure...")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        plot_metrics = [
            m for m in [
                "nk_sm_balance_score",
                "nk_protrusion_machinery_score",
                "nk_topology_permissive_score",
                "nk_synapse_cytotoxicity_outcome_score",
                "checkpoint_link_score",
                "sst_axis_score",
            ]
            if m in scores.columns
        ]
        n_cols = 3
        n_rows = int(np.ceil(len(plot_metrics) / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        flat_axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes

        for ax, metric in zip(flat_axes, plot_metrics):
            for cond in scores["condition"].unique():
                data = scores[scores["condition"] == cond][metric]
                ax.hist(data, bins=50, alpha=0.5, label=str(cond))
            ax.set_title(metric)
            ax.legend(fontsize=8)

        # Hide unused subplots
        for ax in flat_axes[len(plot_metrics):]:
            ax.set_visible(False)

        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "fig9_sst_axis_scores.pdf"), dpi=150)
        plt.close()
        log("  Saved: fig9_sst_axis_scores.pdf")
    except Exception as e:
        log(f"  Plotting failed (non-critical): {e}")

    log("=" * 50)
    log("SST Axis Analysis Complete!")
    log("=" * 50)


if __name__ == "__main__":
    main()
