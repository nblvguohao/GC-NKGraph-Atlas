"""
T12 — Ablation study: quantify contribution of metabolic_crosstalk edge.

Compares full-graph GNN vs no-metabolic-edge GNN on identical 5-fold CV splits.
Reports ΔMCC, ΔAUROC, and paired test results.

This addresses the reviewer question: "If the metabolic edge doesn't improve
classification, what is its value?" The honest answer should be that its value
is in the embedding structure (see T11), not classification accuracy — and this
script provides the evidence for that claim.

RUN ON SERVER:
    python src/baselines/run_ablation.py \
        --full-edges data/processed/graph/edges.tsv \
        --nodes data/processed/graph/nodes.tsv \
        --expression data/processed/bulk/tcga_stad_expression.tsv \
        --labels data/processed/bulk/tcga_stad_labels.tsv \
        --output-dir results/
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))


# ── Helpers ──────────────────────────────────────────────────────────────────

def filter_metabolic_edges(edges_df):
    """Remove rows where edge_type == 'metabolic_crosstalk'."""
    if "edge_type" not in edges_df.columns:
        raise ValueError("edges file must have 'edge_type' column")
    return edges_df[edges_df["edge_type"] != "metabolic_crosstalk"].copy()


def paired_stats(values_a, values_b, metric_name):
    """
    Paired Wilcoxon signed-rank and paired t-test.
    Returns dict with delta, p-values, and 95% CI of delta.
    """
    deltas = np.array(values_a) - np.array(values_b)
    n = len(deltas)
    mean_delta = np.mean(deltas)
    std_delta = np.std(deltas, ddof=1)

    # Paired t-test
    t_stat, t_p = stats.ttest_rel(values_a, values_b)

    # Wilcoxon signed-rank
    try:
        w_stat, w_p = stats.wilcoxon(values_a, values_b)
    except ValueError:
        w_stat, w_p = np.nan, np.nan

    # 95% CI of delta (t-distribution with n-1 df)
    t_crit = stats.t.ppf(0.975, n - 1)
    ci_lo = mean_delta - t_crit * std_delta / np.sqrt(n)
    ci_hi = mean_delta + t_crit * std_delta / np.sqrt(n)

    return {
        "metric": metric_name,
        "n_folds": n,
        "mean_A": np.mean(values_a),
        "mean_B": np.mean(values_b),
        "delta_mean": mean_delta,
        "delta_std": std_delta,
        "ci_95_lower": ci_lo,
        "ci_95_upper": ci_hi,
        "ttest_statistic": t_stat,
        "ttest_p": t_p,
        "wilcoxon_statistic": w_stat,
        "wilcoxon_p": w_p,
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="T12: Graph ablation — full vs no-metabolic-edge GNN"
    )
    parser.add_argument(
        "--full-edges", required=True,
        help="Path to full edges.tsv (all edge types)"
    )
    parser.add_argument(
        "--nodes", required=True,
        help="Path to nodes.tsv"
    )
    parser.add_argument(
        "--expression", required=True,
        help="Path to TCGA-STAD expression matrix"
    )
    parser.add_argument(
        "--labels", required=True,
        help="Path to TCGA-STAD NK state labels"
    )
    parser.add_argument(
        "--output-dir", default="results",
        help="Output directory"
    )
    parser.add_argument(
        "--n-folds", type=int, default=5,
        help="Number of CV folds"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for fold splitting"
    )
    parser.add_argument(
        "--epochs", type=int, default=200,
        help="Max training epochs per fold"
    )
    parser.add_argument(
        "--patience", type=int, default=30,
        help="Early stopping patience"
    )
    parser.add_argument(
        "--lr", type=float, default=1e-3,
        help="Learning rate"
    )
    parser.add_argument(
        "--hidden-dim", type=int, default=128,
        help="GNN hidden dimension"
    )
    parser.add_argument(
        "--embedding-dim", type=int, default=64,
        help="Gene embedding dimension"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "tables"), exist_ok=True)

    log = lambda msg: print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    # ── Load data ────────────────────────────────────────────────────────────
    log("Loading graph edges...")
    edges_full = pd.read_csv(args.full_edges, sep="\t")
    log(f"  Full edges: {len(edges_full)} rows")

    n_metabolic = (edges_full["edge_type"] == "metabolic_crosstalk").sum() \
        if "edge_type" in edges_full.columns else 0
    log(f"  metabolic_crosstalk edges: {n_metabolic}")

    edges_ablated = filter_metabolic_edges(edges_full)
    log(f"  Ablated edges: {len(edges_ablated)} rows "
        f"({len(edges_full) - len(edges_ablated)} removed)")

    log("Loading expression and labels...")
    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    labels = pd.read_csv(args.labels, sep="\t", index_col=0)

    # ── Set up fold splits (manual stratified k-fold, no sklearn dependency) ──
    n_samples = len(expr)
    y = labels.iloc[:, 0].values if labels.ndim > 1 else labels.values
    y = np.array([str(v) for v in y])  # ensure string labels

    # Manual stratified k-fold
    rng = np.random.default_rng(args.seed)
    unique_classes = np.unique(y)
    folds = np.zeros(n_samples, dtype=int)
    for cls in unique_classes:
        cls_idx = np.where(y == cls)[0]
        rng.shuffle(cls_idx)
        # Distribute evenly across folds
        fold_assignments = np.arange(len(cls_idx)) % args.n_folds
        rng.shuffle(fold_assignments)
        folds[cls_idx] = fold_assignments

    log(f"  {n_samples} samples, {args.n_folds} folds, classes={dict(zip(*np.unique(y, return_counts=True)))}")

    # ── Run ablation ─────────────────────────────────────────────────────────
    # This is a skeleton — the actual GNN training call depends on the project's
    # model interface. The pattern shown here matches gc_nkgraph_atlas.py.
    # In practice, we import the model module and call it for each configuration.

    full_results = {"fold": [], "accuracy": [], "balanced_acc": [],
                    "macro_f1": [], "MCC": [], "AUROC": [], "AUPRC": []}
    ablated_results = {"fold": [], "accuracy": [], "balanced_acc": [],
                       "macro_f1": [], "MCC": [], "AUROC": [], "AUPRC": []}

    log("Starting 5-fold CV ablation...")

    for fold_idx in range(args.n_folds):
        test_idx = np.where(folds == fold_idx)[0]
        train_idx = np.where(folds != fold_idx)[0]
        log(f"  Fold {fold_idx + 1}/{args.n_folds} (train={len(train_idx)}, test={len(test_idx)})")

        # --- Full graph ---
        # In production, this would call:
        #   model_full = train_gnn(expr, labels, edges_full, nodes, train_idx, ...)
        #   metrics_full = evaluate(model_full, expr, labels, test_idx)

        # Placeholder: random values for structural testing
        # REMOVE THIS BLOCK and replace with actual GNN calls on the server
        rng = np.random.default_rng(args.seed + fold_idx)
        for metric in ["accuracy", "balanced_acc", "macro_f1", "MCC", "AUROC", "AUPRC"]:
            full_results[metric].append(rng.uniform(0.7, 0.95))
            ablated_results[metric].append(rng.uniform(0.7, 0.95))
        full_results["fold"].append(fold_idx + 1)
        ablated_results["fold"].append(fold_idx + 1)

        log(f"    (placeholder values — replace with actual GNN training)")

    log("Ablation CV complete.")

    # ── NOTE ─────────────────────────────────────────────────────────────────
    print()
    print("!" * 60)
    print("  IMPORTANT: This script currently uses PLACEHOLDER random values.")
    print("  To run the actual ablation, replace the training loop above with:")
    print()
    print("    from src.models.gc_nkgraph_atlas import GC_NKGraph_Atlas")
    print("    model_full = GC_NKGraph_Atlas(...)")
    print("    model_full.fit(expr, labels, edges_full, nodes, train_idx)")
    print("    metrics = model_full.evaluate(test_idx)")
    print()
    print("  Then run model_ablated with edges_ablated instead of edges_full.")
    print("!" * 60)
    print()
    # ──────────────────────────────────────────────────────────────────────────

    # ── Statistical comparison ───────────────────────────────────────────────
    all_stats = []
    for metric in ["MCC", "AUROC"]:
        st = paired_stats(
            full_results[metric], ablated_results[metric], metric
        )
        all_stats.append(st)

    stats_df = pd.DataFrame(all_stats)
    stats_path = os.path.join(args.output_dir, "tables", "ablation_results.tsv")
    stats_df.to_csv(stats_path, sep="\t", index=False, float_format="%.6f")
    log(f"Wrote {stats_path}")

    # ── Per-fold detail ──────────────────────────────────────────────────────
    full_df = pd.DataFrame(full_results)
    full_df["mode"] = "full_graph"
    ablated_df = pd.DataFrame(ablated_results)
    ablated_df["mode"] = "no_metabolic_edge"
    detail_df = pd.concat([full_df, ablated_df], ignore_index=True)
    detail_path = os.path.join(args.output_dir, "tables",
                               "ablation_per_fold.tsv")
    detail_df.to_csv(detail_path, sep="\t", index=False)
    log(f"Wrote {detail_path}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("T12 — ABLATION RESULTS")
    print("=" * 60)
    for st in all_stats:
        print(f"  {st['metric']}:")
        print(f"    Full:  {st['mean_A']:.4f}")
        print(f"    Ablated: {st['mean_B']:.4f}")
        print(f"    Δ:     {st['delta_mean']:+.4f} "
              f"[{st['ci_95_lower']:+.4f}, {st['ci_95_upper']:+.4f}]")
        print(f"    Paired t-test: p={st['ttest_p']:.4f}")
        print(f"    Wilcoxon: p={st['wilcoxon_p']:.4f}")
    print()
    delta_mcc = all_stats[0]["delta_mean"]
    if abs(delta_mcc) < 0.02:
        print("  VERDICT: Removing metabolic_crosstalk edge does not")
        print("  significantly change classification performance. The edge's")
        print("  value lies in the embedding structure (see T11), not in")
        print("  classification accuracy. This supports the 'comparable")
        print("  accuracy, added interpretability' framing in §3.4.")
    else:
        print(f"  VERDICT: ΔMCC = {delta_mcc:+.4f} — the metabolic edge")
        print("  contributes measurably to classification. Report the")
        print("  specific contribution in §3.7.")
    print("=" * 60)

    log("T12 complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
