"""
Self-contained synthetic ablation study.

Runs the full edge-type ablation using LOCAL SYNTHETIC DATA so reviewers
can verify the entire code path end-to-end without real TCGA/GEO data.

Produces:
  results/tables/ablation_edge_types_SYNTHETIC.tsv

Usage:
  python src/baselines/run_ablation_synthetic.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir
from src.common.seed import set_seed


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================================
# Synthetic data setup
# ============================================================================

def build_synthetic_graph(n_genes: int = 500, seed: int = 42) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build a realistic synthetic gene graph with 6 edge types."""
    rng = np.random.RandomState(seed)

    # Node types: mostly genes, a few special nodes
    gene_ids = [f"GENE_{i}" for i in range(n_genes)]
    node_types = ["gene"] * n_genes
    # Add some NK receptor and tumor program nodes
    nk_receptors = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "KLRD1", "NCR1", "NCAM1",
                    "HAVCR2", "TIGIT", "CD96", "KLRC1"]
    tumor_serine = ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2"]
    for name in nk_receptors + tumor_serine:
        if name not in gene_ids:
            gene_ids.append(name)
            node_types.append("gene")
    n_nodes = len(gene_ids)

    nodes = pd.DataFrame({"node_id": gene_ids, "node_type": node_types})

    # Edge types with realistic edge counts
    edge_configs = {
        "ppi": (n_nodes * 3, 0.7),           # PPI: dense, many edges
        "ligand_receptor": (n_nodes // 2, 0.9),
        "tf_target": (n_nodes, 0.8),
        "metabolic_crosstalk": (n_nodes // 3, 0.5),
        "sm_topology_axis": (n_nodes // 2, 0.3),
        "dysfunction_correlation": (n_nodes // 3, 0.6),
    }

    edges_list = []
    gene_idx_map = {g: i for i, g in enumerate(gene_ids)}
    for etype, (n_edges_target, weight) in edge_configs.items():
        n_edges = min(n_edges_target, n_nodes * (n_nodes - 1) // 4)
        existing = set()
        for _ in range(n_edges):
            s, d = int(rng.randint(0, n_nodes)), int(rng.randint(0, n_nodes))
            if s == d or (s, d) in existing or (d, s) in existing:
                continue
            existing.add((s, d))
            edges_list.append({
                "src": gene_ids[s],
                "dst": gene_ids[d],
                "edge_type": etype,
                "weight": round(rng.uniform(0.3, 1.0), 3),
            })

    edges = pd.DataFrame(edges_list)
    return nodes, edges


def build_synthetic_expression(
    nodes: pd.DataFrame,
    n_samples: int = 300,
    n_genes_per_sample: int = 400,
    seed: int = 42,
) -> pd.DataFrame:
    """Build synthetic expression matrix matching the gene graph."""
    rng = np.random.RandomState(seed)
    gene_ids = nodes["node_id"].tolist()
    n_genes = len(gene_ids)

    # Base expression from log-normal
    expr = pd.DataFrame(
        rng.lognormal(mean=2.0, sigma=0.8, size=(n_samples, n_genes)),
        columns=gene_ids,
        index=[f"SYN_SAMPLE_{i:04d}" for i in range(n_samples)],
    )

    # Make some genes more variable (simulate real biology)
    for g in ["GZMB", "PRF1", "NKG7", "IFNG"]:
        if g in expr.columns:
            expr[g] = expr[g] * rng.uniform(0.5, 3.0, n_samples)

    return expr


def build_synthetic_labels(
    expr: pd.DataFrame,
    seed: int = 42,
) -> pd.DataFrame:
    """Build synthetic NK state labels correlated with cytotoxic gene expression."""
    rng = np.random.RandomState(seed)
    n_samples = len(expr)

    # Create a cytotoxicity score from known markers
    cyto_genes = ["GZMB", "PRF1", "NKG7", "GNLY", "IFNG"]
    available = [g for g in cyto_genes if g in expr.columns]
    if available:
        cyto_score = expr[available].mean(axis=1)
        # Top 40% are "NK-hot-cytotoxic"
        threshold = cyto_score.quantile(0.6)
        labels = pd.DataFrame({
            "nk_immune_state": np.where(cyto_score >= threshold, "NK-hot-cytotoxic", "NK-cold/excluded"),
        }, index=expr.index)
    else:
        labels = pd.DataFrame({
            "nk_immune_state": rng.choice(["NK-hot-cytotoxic", "NK-cold/excluded"], n_samples),
        }, index=expr.index)

    return labels


# ============================================================================
# GNN training (adapted to work with synthetic data directly)
# ============================================================================

def train_gnn_synthetic(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    expr: pd.DataFrame,
    labels: pd.DataFrame,
    embedding_dim: int = 64,
    epochs: int = 100,
) -> Dict[str, float]:
    """Train GNN on synthetic data and return 5-fold CV metrics."""
    from src.models.gc_nkgraph_atlas import GeneGraphEncoder, NKStateClassifier, evaluate
    from sklearn.model_selection import StratifiedKFold

    set_seed(42)

    # Build graph files in temp dir
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        nodes.to_csv(os.path.join(tmpdir, "nodes.tsv"), sep="\t", index=False)
        edges.to_csv(os.path.join(tmpdir, "edges.tsv"), sep="\t", index=False)

        # Stage 1: Gene embeddings
        encoder = GeneGraphEncoder(graph_dir=tmpdir, embedding_dim=embedding_dim)
        encoder.fit()

        # Match genes
        graph_genes = set(encoder._gene_to_idx.keys())
        common_genes = [g for g in expr.columns if g in graph_genes]
        if len(common_genes) < 5:
            return {"MCC": float("nan"), "AUROC": float("nan"), "n_common_genes": 0}

        X_expr = expr[common_genes].values.astype(np.float32)
        gene_embeddings = encoder.transform(common_genes)

        y_full = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

        # 5-fold CV
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        results = []
        for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y_full)):
            n_train = len(train_idx)
            val_size = max(1, int(n_train * 0.2))
            rng = np.random.RandomState(42 + fold)
            val_idx = rng.choice(train_idx, size=val_size, replace=False)
            train_clean = np.setdiff1d(train_idx, val_idx)

            clf = NKStateClassifier(
                embedding_dim=embedding_dim, hidden_dims=[128, 64],
                num_classes=2, dropout=0.3, learning_rate=1e-3,
            )
            clf.fit(X_expr[train_clean], y_full[train_clean], gene_embeddings,
                    X_val_expr=X_expr[val_idx], y_val=y_full[val_idx],
                    epochs=epochs, batch_size=32, verbose=False)
            y_pred = clf.predict(X_expr[test_idx], gene_embeddings)
            y_prob = clf.predict_proba(X_expr[test_idx], gene_embeddings)
            m = evaluate(y_full[test_idx], y_pred, y_prob)
            m["fold"] = fold
            results.append(m)

    df = pd.DataFrame(results)
    mean_metrics = df.drop(columns=["fold"]).mean().to_dict()
    mean_metrics["n_common_genes"] = len(common_genes)
    return mean_metrics


# ============================================================================
# Main
# ============================================================================

def main():
    log("=" * 60)
    log("SYNTHETIC ABLATION STUDY (Self-Contained)")
    log("Produces real numeric results from local synthetic data")
    log("=" * 60)

    out_dir = ensure_dir("results/tables")
    set_seed(42)

    # Build synthetic data
    log("\n--- Building synthetic gene graph ---")
    nodes, edges = build_synthetic_graph(n_genes=500)
    log(f"  Nodes: {len(nodes)}")
    log(f"  Edges: {len(edges)}")
    all_edge_types = sorted(edges["edge_type"].unique())
    for et in all_edge_types:
        n = int((edges["edge_type"] == et).sum())
        log(f"    {et}: {n} edges")

    log("\n--- Building synthetic expression + labels ---")
    expr = build_synthetic_expression(nodes, n_samples=300)
    labels = build_synthetic_labels(expr)
    y = (labels["nk_immune_state"] == "NK-hot-cytotoxic").sum()
    log(f"  Samples: {len(expr)}, Genes: {expr.shape[1]}")
    log(f"  NK-hot-cytotoxic: {y}/{len(labels)}")

    # ============================================================
    # EXPERIMENT 1: Full graph baseline
    # ============================================================
    log("\n" + "=" * 60)
    log("EXPERIMENT 1: Full Graph")
    log("=" * 60)
    metrics_full = train_gnn_synthetic(nodes, edges, expr, labels, embedding_dim=64, epochs=100)
    baseline_mcc = metrics_full.get("MCC", 0)
    log(f"  Full graph: MCC={baseline_mcc:.4f}  AUROC={metrics_full.get('AUROC', 0):.4f}  "
        f"Acc={metrics_full.get('Accuracy', 0):.4f}")

    all_results = [{
        "condition": "full_graph",
        "edge_type_removed": "none",
        "n_edge_types": len(all_edge_types),
        "n_edges": len(edges),
        **metrics_full,
    }]

    # ============================================================
    # EXPERIMENT 2: Remove one edge type at a time
    # ============================================================
    log("\n" + "=" * 60)
    log("EXPERIMENT 2: Leave-One-Out Ablation")
    log("=" * 60)

    for etype in all_edge_types:
        log(f"\n  Removing: {etype}")
        ablated = edges[edges["edge_type"] != etype].copy()
        remaining = sorted(ablated["edge_type"].unique())
        try:
            m = train_gnn_synthetic(nodes, ablated, expr, labels, embedding_dim=64, epochs=100)
            delta = m.get("MCC", 0) - baseline_mcc
            direction = "v" if delta < -0.005 else ("^" if delta > 0.005 else "-")
            log(f"    MCC={m.get('MCC', 0):.4f} (delta={delta:+.4f}) {direction}  "
                f"AUROC={m.get('AUROC', 0):.4f}")
        except Exception as e:
            log(f"    FAILED: {e}")
            m = {"MCC": float("nan"), "AUROC": float("nan")}
            delta = float("nan")
        all_results.append({
            "condition": f"remove_{etype}",
            "edge_type_removed": etype,
            "n_edge_types": len(remaining),
            "n_edges": len(ablated),
            "delta_mcc": delta,
            **m,
        })

    # ============================================================
    # EXPERIMENT 3: Single edge type only
    # ============================================================
    log("\n" + "=" * 60)
    log("EXPERIMENT 3: Single-Edge-Type Graphs")
    log("=" * 60)

    for etype in all_edge_types:
        log(f"\n  Only: {etype}")
        single = edges[edges["edge_type"] == etype].copy()
        if len(single) == 0:
            log(f"    No edges — skipping")
            all_results.append({
                "condition": f"only_{etype}", "edge_type_removed": etype,
                "n_edge_types": 0, "n_edges": 0,
                "MCC": float("nan"), "AUROC": float("nan"),
            })
            continue
        try:
            m = train_gnn_synthetic(nodes, single, expr, labels, embedding_dim=64, epochs=100)
            delta = m.get("MCC", 0) - baseline_mcc
            log(f"    MCC={m.get('MCC', 0):.4f} (delta={delta:+.4f})  "
                f"AUROC={m.get('AUROC', 0):.4f}")
        except Exception as e:
            log(f"    FAILED: {e}")
            m = {"MCC": float("nan"), "AUROC": float("nan")}
            delta = float("nan")
        all_results.append({
            "condition": f"only_{etype}",
            "edge_type_removed": etype,
            "n_edge_types": 1,
            "n_edges": len(single),
            "delta_mcc": delta,
            **m,
        })

    # ============================================================
    # Results summary
    # ============================================================
    res_df = pd.DataFrame(all_results)
    res_df["data_source"] = "SYNTHETIC_ABLATION"
    res_df = res_df.sort_values("MCC", ascending=False, na_position="last").reset_index(drop=True)

    out_path = os.path.join(out_dir, "ablation_edge_types_SYNTHETIC.tsv")
    res_df.to_csv(out_path, sep="\t", index=False)

    log("\n" + "=" * 60)
    log("ABLATION RESULTS (Synthetic Data)")
    log("=" * 60)
    display = ["condition", "n_edge_types", "n_edges", "MCC", "AUROC", "delta_mcc"]
    avail = [c for c in display if c in res_df.columns]
    print(res_df[avail].to_string(index=False))

    # ---- Analysis ----
    log("\n--- Analysis ---")
    log(f"  Baseline MCC: {baseline_mcc:.4f}")

    remove_rows = res_df[res_df["condition"].str.startswith("remove_")].copy()
    if len(remove_rows) > 0 and not np.isnan(baseline_mcc):
        remove_rows["mcc_drop"] = baseline_mcc - remove_rows["MCC"].fillna(baseline_mcc)
        remove_rows = remove_rows.sort_values("mcc_drop", ascending=False)
        log("\n  Edge types ranked by MCC DROP when REMOVED:")
        for _, r in remove_rows.iterrows():
            drop = r["mcc_drop"]
            bar = "#" * max(0, int(drop * 200))
            log(f"    {r['edge_type_removed']:<25}  drop={drop:+.4f}  {bar}")

    single_rows = res_df[res_df["condition"].str.startswith("only_")].copy()
    single_rows = single_rows.sort_values("MCC", ascending=False, na_position="last")
    log("\n  Edge types ranked by MCC when used ALONE:")
    for _, r in single_rows.iterrows():
        log(f"    {r['edge_type_removed']:<25}  MCC={r['MCC']:.4f}")

    log(f"\n  Results saved: {out_path}")
    log(f"\n  NOTE: Results are from SYNTHETIC data for code-path verification.")
    log(f"  To run with REAL data: python src/baselines/run_ablation.py")
    log("=" * 60)
    log("SYNTHETIC ABLATION COMPLETE")
    log("=" * 60)


if __name__ == "__main__":
    main()
