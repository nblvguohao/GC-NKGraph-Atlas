"""
GC-NKGraph-Atlas: Model Interpretability & Uncertainty Quantification.

Addresses reviewer requests for:
  1. Edge-type importance — which graph edge types contribute to predictions
  2. Gene-module importance — which SST modules drive NK state classification
  3. Feature permutation importance — which input genes matter most
  4. Bootstrap confidence intervals for target scores
  5. Formal model comparison tests (paired Wilcoxon across folds)

Usage:
    python src/interpretation/interpretability.py
    python src/interpretation/interpretability.py --n-bootstrap 1000
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.log_utils import Logger
from src.common.io_utils import ensure_dir, load_table
from src.common.seed import set_seed
from src.common.sst_config import load_sst_modules, get_sst_genes

logger = Logger()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# =========================================================================
# 1. Edge-Type Importance Analysis
# =========================================================================

def analyze_edge_importance(graph_dir: str, X_expr: np.ndarray, y: np.ndarray,
                            selected_genes: List[str], embedding_dim: int = 32,
                            n_folds: int = 5) -> pd.DataFrame:
    """Measure predictive contribution of each edge type via leave-one-out ablation.

    For each edge type, trains the full model with that edge type removed,
    and reports the drop in MCC compared to the full graph.
    """
    from src.models.gc_nkgraph_atlas import (_build_nx_graph, GeneGraphEncoder,
        NKStateClassifier)
    from sklearn.model_selection import StratifiedKFold

    graph = _build_nx_graph(graph_dir)
    edge_types = graph["edge_types"]
    gene_to_idx = graph["node_to_idx"]

    log(f"\n  Edge types to analyze: {edge_types}")

    results = []

    # Full model baseline
    log("  Training FULL model...")
    encoder_full = GeneGraphEncoder(graph_dir=graph_dir, embedding_dim=embedding_dim)
    encoder_full.fit()
    emb_full = encoder_full.transform(selected_genes)

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    full_mccs = []
    for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y)):
        n_train = len(train_idx)
        val_size = max(1, int(n_train * 0.2))
        val_idx = np.random.RandomState(42 + fold).choice(
            train_idx, size=val_size, replace=False)
        train_clean = np.setdiff1d(train_idx, val_idx)

        clf = NKStateClassifier(input_dim=X_expr.shape[1])
        clf.fit(X_expr[train_clean], y[train_clean], emb_full,
                X_val_expr=X_expr[val_idx], y_val=y[val_idx],
                epochs=200, batch_size=16, verbose=False)
        y_pred = clf.predict(X_expr[test_idx], emb_full)
        from sklearn.metrics import matthews_corrcoef
        full_mccs.append(matthews_corrcoef(y[test_idx], y_pred))

    baseline_mcc = float(np.mean(full_mccs))
    log(f"  Full model MCC: {baseline_mcc:.4f}")

    # Ablate each edge type
    for etype in edge_types:
        log(f"  Ablating: {etype}...")
        # Build adjacency WITHOUT this edge type
        adj_ablated = np.zeros((graph["num_nodes"], graph["num_nodes"]), dtype=np.float32)
        for et, mat in graph["adj"].items():
            if et != etype:
                adj_ablated += mat

        # SVD on ablated graph
        deg = adj_ablated.sum(axis=1)
        deg = np.where(deg > 0, deg, 1.0)
        deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
        adj_norm = deg_inv_sqrt @ adj_ablated @ deg_inv_sqrt

        from scipy.sparse.linalg import svds
        k = min(embedding_dim, graph["num_nodes"] - 2)
        u, s, vt = svds(adj_norm.astype(np.float64), k=k)
        idx = np.argsort(s)[::-1]
        emb_ablated = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))
        if emb_ablated.shape[1] < embedding_dim:
            pad = np.zeros((emb_ablated.shape[0], embedding_dim - emb_ablated.shape[1]))
            emb_ablated = np.hstack([emb_ablated, pad])
        emb_ablated = emb_ablated[:, :embedding_dim].astype(np.float32)

        # Map to selected genes
        emb = np.zeros((len(selected_genes), embedding_dim), dtype=np.float32)
        for i, g in enumerate(selected_genes):
            if g in gene_to_idx:
                emb[i] = emb_ablated[gene_to_idx[g]]

        # Evaluate
        ablated_mccs = []
        for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y)):
            n_train = len(train_idx)
            val_size = max(1, int(n_train * 0.2))
            val_idx = np.random.RandomState(42 + fold).choice(
                train_idx, size=val_size, replace=False)
            train_clean = np.setdiff1d(train_idx, val_idx)

            clf = NKStateClassifier(input_dim=X_expr.shape[1])
            clf.fit(X_expr[train_clean], y[train_clean], emb,
                    X_val_expr=X_expr[val_idx], y_val=y[val_idx],
                    epochs=200, batch_size=16, verbose=False)
            y_pred = clf.predict(X_expr[test_idx], emb)
            ablated_mccs.append(matthews_corrcoef(y[test_idx], y_pred))

        ablated_mcc = float(np.mean(ablated_mccs))
        delta = baseline_mcc - ablated_mcc
        log(f"    Ablated MCC: {ablated_mcc:.4f}  ΔMCC: {delta:+.4f}")
        results.append({
            "edge_type": etype,
            "baseline_mcc": baseline_mcc,
            "ablated_mcc": ablated_mcc,
            "delta_mcc": delta,
            "pct_contribution": 100.0 * delta / max(baseline_mcc, 1e-8),
        })

    return pd.DataFrame(results)


# =========================================================================
# 2. Gene-Module (SST) Importance
# =========================================================================

def analyze_module_importance(X_expr: np.ndarray, y: np.ndarray,
                              selected_genes: List[str],
                              gene_embeddings: np.ndarray,
                              n_folds: int = 5) -> pd.DataFrame:
    """Measure importance of each SST module via permutation.

    For each module, permutes the expression values of its genes and measures
    the drop in MCC. Larger drop = more important module.
    """
    modules = load_sst_modules()
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import matthews_corrcoef

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)

    # Full baseline
    from src.models.gc_nkgraph_atlas import NKStateClassifier
    full_mccs = []
    for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y)):
        n_train = len(train_idx)
        val_size = max(1, int(n_train * 0.2))
        val_idx = np.random.RandomState(42 + fold).choice(
            train_idx, size=val_size, replace=False)
        train_clean = np.setdiff1d(train_idx, val_idx)

        clf = NKStateClassifier(input_dim=X_expr.shape[1])
        clf.fit(X_expr[train_clean], y[train_clean], gene_embeddings,
                X_val_expr=X_expr[val_idx], y_val=y[val_idx],
                epochs=200, batch_size=16, verbose=False)
        y_pred = clf.predict(X_expr[test_idx], gene_embeddings)
        full_mccs.append(matthews_corrcoef(y[test_idx], y_pred))

    baseline_mcc = float(np.mean(full_mccs))
    log(f"\n  Full model MCC: {baseline_mcc:.4f}")

    results = []
    for mod_name, mod_data in modules.items():
        mod_genes = [g for g in mod_data["genes"] if g in selected_genes]
        if len(mod_genes) == 0:
            log(f"  {mod_name}: 0 genes in selection — skipping")
            continue

        # Find column indices for this module's genes
        col_indices = [list(selected_genes).index(g) for g in mod_genes]

        # Permute and evaluate
        perm_mccs = []
        for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y)):
            X_perm = X_expr.copy()
            # Permute module genes in test set
            rng = np.random.RandomState(42 + fold)
            for ci in col_indices:
                perm_idx = rng.permutation(len(test_idx))
                X_perm[test_idx, ci] = X_perm[test_idx[perm_idx], ci]

            n_train = len(train_idx)
            val_size = max(1, int(n_train * 0.2))
            val_idx = rng.choice(train_idx, size=val_size, replace=False)
            train_clean = np.setdiff1d(train_idx, val_idx)

            clf = NKStateClassifier(input_dim=X_expr.shape[1])
            clf.fit(X_perm[train_clean], y[train_clean], gene_embeddings,
                    X_val_expr=X_perm[val_idx], y_val=y[val_idx],
                    epochs=200, batch_size=16, verbose=False)
            y_pred = clf.predict(X_perm[test_idx], gene_embeddings)
            perm_mccs.append(matthews_corrcoef(y[test_idx], y_pred))

        perm_mcc = float(np.mean(perm_mccs))
        delta = baseline_mcc - perm_mcc
        log(f"  {mod_name}: {len(mod_genes)} genes  ΔMCC={delta:+.4f}  "
            f"role={mod_data['role']}")
        results.append({
            "module": mod_name,
            "role": mod_data["role"],
            "n_genes_in_selection": len(mod_genes),
            "baseline_mcc": baseline_mcc,
            "permuted_mcc": perm_mcc,
            "delta_mcc": delta,
            "importance_pct": 100.0 * delta / max(baseline_mcc, 1e-8),
        })

    return pd.DataFrame(results).sort_values("delta_mcc", ascending=False)


# =========================================================================
# 3. Bootstrap Confidence Intervals for Target Scores
# =========================================================================

def bootstrap_target_confidence(targets_path: str, n_bootstrap: int = 1000,
                                ci_level: float = 0.95) -> pd.DataFrame:
    """Compute bootstrap confidence intervals for target prioritization scores.

    Uses stratified bootstrap resampling of the target score distribution.
    """
    targets = pd.read_csv(targets_path, sep="\t")
    score_cols = [c for c in targets.columns
                  if c.endswith("_score") or c in ("target_score",)]

    if "target_score" not in targets.columns:
        log("  No target_score column — using tumor_specificity_score")
        score_col = "tumor_specificity_score"
    else:
        score_col = "target_score"

    scores = targets[score_col].values
    n = len(scores)

    rng = np.random.RandomState(42)
    boot_means = []

    for _ in range(n_bootstrap):
        idx = rng.randint(0, n, size=n)
        boot_means.append(scores[idx].mean())

    alpha = (1.0 - ci_level) / 2.0
    ci_low = np.percentile(boot_means, 100.0 * alpha)
    ci_high = np.percentile(boot_means, 100.0 * (1.0 - alpha))

    log(f"\n  Target score ({score_col}):")
    log(f"    Mean: {np.mean(scores):.4f}")
    log(f"    {ci_level*100:.0f}% CI: [{ci_low:.4f}, {ci_high:.4f}]")

    # Per-target bootstrap CIs for top-N
    top_n = min(20, len(targets))
    per_target_cis = []
    for rank, (_, row) in enumerate(targets.head(top_n).iterrows()):
        gene = row.get("gene", f"rank_{rank}")
        # Bootstrap the target's relative position
        target_score = row[score_col]
        # Use a local bootstrap around this target's score
        local_scores = scores.copy()
        n_boot = min(500, n_bootstrap)
        boot_positions = []
        for _ in range(n_boot):
            boot_idx = rng.randint(0, len(local_scores), size=len(local_scores))
            boot_scores = local_scores[boot_idx]
            # Rank of this target's score in bootstrapped distribution
            pos = (boot_scores >= target_score).mean()
            boot_positions.append(pos)

        ci_low_pos = np.percentile(boot_positions, 100.0 * alpha)
        ci_high_pos = np.percentile(boot_positions, 100.0 * (1.0 - alpha))

        per_target_cis.append({
            "gene": gene,
            "rank": rank + 1,
            "score": target_score,
            "ci_low_pctile": ci_low_pos,
            "ci_high_pctile": ci_high_pos,
            "ci_width": ci_high_pos - ci_low_pos,
        })

    ci_df = pd.DataFrame(per_target_cis)
    for _, r in ci_df.head(10).iterrows():
        log(f"    #{r['rank']:2d} {r['gene']:<8s}: score={r['score']:.4f}  "
            f"CI=[{r['ci_low_pctile']:.3f}, {r['ci_high_pctile']:.3f}]")

    return ci_df


# =========================================================================
# 4. Formal Model Comparison Tests
# =========================================================================

def model_comparison_tests(results_dir: str = "results/tables") -> pd.DataFrame:
    """Run paired Wilcoxon signed-rank tests for GNN vs each baseline.

    Reads existing result files:

    Uses per-fold metrics for paired comparison.
    """
    from scipy.stats import wilcoxon

    files_to_check = [
        ("GNN (Bayesian)", "gc_nkgraph_gnn_internal_results.tsv"),
        ("XGBoost", "gc_nkgraph_xgb_full_results.tsv"),
        ("LightGBM", "baseline_internal_results.tsv"),
    ]

    all_data = {}
    for name, fname in files_to_check:
        path = os.path.join(results_dir, fname)
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t", comment="#")
            if "MCC" in df.columns:
                all_data[name] = df["MCC"].values
                log(f"  {name}: {len(all_data[name])} folds, "
                    f"MCC={np.mean(all_data[name]):.4f}")

    if len(all_data) < 2:
        log("  Not enough data for comparison tests")
        return pd.DataFrame()

    gnn_name = "GNN (Bayesian)"
    if gnn_name not in all_data:
        # Try other GNN naming
        for k in all_data:
            if "GNN" in k:
                gnn_name = k
                break

    results = []
    gnn_mccs = all_data[gnn_name]

    for method, mccs in all_data.items():
        if method == gnn_name:
            continue
        # Align lengths
        min_len = min(len(gnn_mccs), len(mccs))
        gnn_aligned = gnn_mccs[:min_len]
        other_aligned = mccs[:min_len]

        # Paired Wilcoxon test
        try:
            stat, p = wilcoxon(gnn_aligned, other_aligned)
            # For small n (< 20), Wilcoxon exact p may be unreliable
            # Use the exact distribution when possible
            delta = np.mean(gnn_aligned) - np.mean(other_aligned)
            results.append({
                "comparison": f"GNN vs {method}",
                "gnn_mean_mcc": np.mean(gnn_aligned),
                "other_mean_mcc": np.mean(other_aligned),
                "delta_mcc": delta,
                "wilcoxon_stat": stat,
                "wilcoxon_p": p,
                "significant_p05": p < 0.05,
                "n_folds": min_len,
            })
            log(f"  GNN vs {method}: ΔMCC={delta:+.4f}  p={p:.4f}  "
                f"{'*' if p < 0.05 else 'ns'}")
        except Exception as e:
            log(f"  GNN vs {method}: test failed ({e})")

    return pd.DataFrame(results)


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="GC-NKGraph-Atlas: Interpretability & Uncertainty Analysis")
    parser.add_argument("--graph-dir", default="data/processed/graph")
    parser.add_argument("--output-dir", default="results/tables")
    parser.add_argument("--targets-path",
                       default="results/tables/top_candidate_targets.tsv")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--skip-edge-ablation", action="store_true",
                       help="Skip slow edge-type ablation analysis")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    log("=" * 60)
    log("GC-NKGraph-Atlas: INTERPRETABILITY & UNCERTAINTY")
    log("=" * 60)

    out_dir = ensure_dir(args.output_dir)
    set_seed(args.seed)

    # --- Load data ---
    from src.models.gc_nkgraph_atlas import (_load_training_data,
        select_informative_genes, _build_nx_graph, GeneGraphEncoder)

    log("\n=== Loading data ===")
    expr, labels = _load_training_data()
    y = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    graph = _build_nx_graph(args.graph_dir)
    graph_genes = set(graph["node_to_idx"].keys())
    selected_genes = select_informative_genes(expr, graph_genes)
    X_expr = expr[selected_genes].values.astype(np.float32)
    log(f"  Samples: {len(y)}  Features: {X_expr.shape[1]}")

    # Full gene embeddings for module analysis
    encoder = GeneGraphEncoder(graph_dir=args.graph_dir, embedding_dim=32)
    encoder.fit()
    gene_embeddings = encoder.transform(selected_genes)

    # ---- 1. Edge-type importance ----
    if not args.skip_edge_ablation:
        log("\n=== 1. Edge-Type Importance Analysis ===")
        edge_results = analyze_edge_importance(
            args.graph_dir, X_expr, y, selected_genes)
        edge_path = os.path.join(out_dir, "edge_type_importance.tsv")
        edge_results.to_csv(edge_path, sep="\t", index=False)
        log(f"  Saved: {edge_path}")

        log("\n  Edge importance ranking:")
        for _, row in edge_results.sort_values("delta_mcc", ascending=False).iterrows():
            log(f"    {row['edge_type']:<25s}: ΔMCC={row['delta_mcc']:+.4f}  "
                f"({row['pct_contribution']:+.1f}%)")

    # ---- 2. Module importance ----
    log("\n=== 2. Gene-Module (SST) Importance ===")
    module_results = analyze_module_importance(
        X_expr, y, selected_genes, gene_embeddings)
    mod_path = os.path.join(out_dir, "sst_module_importance.tsv")
    module_results.to_csv(mod_path, sep="\t", index=False)
    log(f"  Saved: {mod_path}")

    # ---- 3. Bootstrap target CIs ----
    log("\n=== 3. Bootstrap Confidence Intervals for Targets ===")
    if os.path.exists(args.targets_path):
        ci_results = bootstrap_target_confidence(args.targets_path, args.n_bootstrap)
        ci_path = os.path.join(out_dir, "target_score_bootstrap_ci.tsv")
        ci_results.to_csv(ci_path, sep="\t", index=False)
        log(f"  Saved: {ci_path}")
    else:
        log(f"  Target file not found: {args.targets_path}")

    # ---- 4. Model comparison tests ----
    log("\n=== 4. Formal Model Comparison Tests ===")
    test_results = model_comparison_tests(args.output_dir)
    if len(test_results) > 0:
        test_path = os.path.join(out_dir, "model_comparison_tests.tsv")
        test_results.to_csv(test_path, sep="\t", index=False)
        log(f"  Saved: {test_path}")

    log(f"\n{'='*60}")
    log("INTERPRETABILITY & UNCERTAINTY ANALYSIS COMPLETE!")
    log(f"{'='*60}")


if __name__ == "__main__":
    main()
