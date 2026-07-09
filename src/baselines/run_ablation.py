"""
Edge-Type Ablation Study for GC-NKGraph-Atlas Heterogeneous Graph.

Systematically removes each edge type from the heterogeneous gene graph and
measures the GNN model's performance drop. Answers the question:
    "Which edge types contribute to model performance, and by how much?"

Ablation conditions:
  1. Full graph (all 6 edge types) — baseline
  2. Remove PPI edges
  3. Remove ligand_receptor edges
  4. Remove tf_target edges
  5. Remove metabolic_crosstalk edges (the mechanism-specific edge)
  6. Remove sm_topology_axis edges
  7. Remove dysfunction_correlation edges
  8. Each edge type alone (single-edge-type graphs)

Usage:
    python src/baselines/run_ablation.py
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
from src.common.log_utils import Logger  # noqa: E402
from src.common.io_utils import ensure_dir, load_table, load_config  # noqa: E402
from src.common.seed import set_seed  # noqa: E402


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Graph loading & manipulation (no dependency on model internals)
# ---------------------------------------------------------------------------

def load_graph(graph_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load graph nodes and edges."""
    nodes = pd.read_csv(os.path.join(graph_dir, "nodes.tsv"), sep="\t")
    edges = pd.read_csv(os.path.join(graph_dir, "edges.tsv"), sep="\t")
    return nodes, edges


def get_edge_types(edges: pd.DataFrame) -> List[str]:
    """Return sorted list of edge types in the graph."""
    return sorted(edges["edge_type"].unique().tolist())


def ablate_edge_type(
    edges: pd.DataFrame,
    edge_type_to_remove: str,
) -> pd.DataFrame:
    """Return a copy of edges with one edge type removed."""
    return edges[edges["edge_type"] != edge_type_to_remove].copy()


def keep_only_edge_type(
    edges: pd.DataFrame,
    edge_type_to_keep: str,
) -> pd.DataFrame:
    """Return edges containing only one edge type."""
    return edges[edges["edge_type"] == edge_type_to_keep].copy()


def save_graph(nodes: pd.DataFrame, edges: pd.DataFrame, out_dir: str) -> None:
    """Save nodes and edges TSV files to a directory."""
    os.makedirs(out_dir, exist_ok=True)
    nodes.to_csv(os.path.join(out_dir, "nodes.tsv"), sep="\t", index=False)
    edges.to_csv(os.path.join(out_dir, "edges.tsv"), sep="\t", index=False)


def count_edges_by_type(edges: pd.DataFrame) -> Dict[str, int]:
    """Return edge count per type."""
    return edges["edge_type"].value_counts().to_dict()


# ---------------------------------------------------------------------------
# GNN training wrapper (uses existing model code)
# ---------------------------------------------------------------------------

def train_gnn_with_graph(
    graph_dir: str,
    embedding_dim: int = 128,
    epochs: int = 150,
    output_dir: Optional[str] = None,
) -> Dict[str, float]:
    """Train the GNN model using a specific graph directory and return metrics.

    This function mirrors the logic in gc_nkgraph_atlas.py but operates on
    an arbitrary graph directory without command-line argument parsing.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from src.models.gc_nkgraph_atlas import (
        GeneGraphEncoder,
        NKStateClassifier,
        evaluate,
    )
    from sklearn.model_selection import StratifiedKFold

    set_seed(42)

    # Load graph
    if not os.path.exists(os.path.join(graph_dir, "nodes.tsv")):
        raise FileNotFoundError(f"Graph not found at {graph_dir}")

    # Stage 1: Gene embeddings
    encoder = GeneGraphEncoder(
        graph_dir=graph_dir,
        embedding_dim=embedding_dim,
    )
    encoder.fit()

    # Load training data
    config = load_config("configs/data_config.yaml")
    expr = None
    for ds in config.get("bulk_datasets", []):
        if ds["role"] == "train_primary":
            expr = load_table(ds["expression_path"])
            break
    if expr is None:
        raise ValueError("No train_primary dataset found")

    labels = load_table("results/tables/nk_state_labels.tsv")
    common = expr.index.intersection(labels.index)
    expr = expr.loc[common]
    labels = labels.loc[common]
    y_full = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    # Match genes
    graph_genes = set(encoder._gene_to_idx.keys())
    common_genes = [g for g in expr.columns if g in graph_genes]
    if len(common_genes) < 5:
        log(f"    WARNING: Only {len(common_genes)} common genes — returning NaN metrics")
        return {"Accuracy": float("nan"), "MacroF1": float("nan"), "MCC": float("nan"),
                "AUROC": float("nan"), "AUPRC": float("nan"), "n_common_genes": len(common_genes)}

    X_expr = expr[common_genes].values.astype(np.float32)
    gene_embeddings = encoder.transform(common_genes)

    # 5-fold CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    all_results: List[Dict] = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y_full)):
        n_train = len(train_idx)
        val_size = max(1, int(n_train * 0.2))
        rng = np.random.RandomState(42 + fold)
        val_idx = rng.choice(train_idx, size=val_size, replace=False)
        train_idx_clean = np.setdiff1d(train_idx, val_idx)

        classifier = NKStateClassifier(
            embedding_dim=embedding_dim,
            hidden_dims=[256, 128],
            num_classes=2,
            dropout=0.3,
            learning_rate=1e-3,
        )
        classifier.fit(
            X_expr[train_idx_clean], y_full[train_idx_clean],
            gene_embeddings,
            X_val_expr=X_expr[val_idx], y_val=y_full[val_idx],
            epochs=epochs,
            batch_size=32,
            verbose=False,
        )

        y_pred = classifier.predict(X_expr[test_idx], gene_embeddings)
        y_prob = classifier.predict_proba(X_expr[test_idx], gene_embeddings)
        metrics = evaluate(y_full[test_idx], y_pred, y_prob)
        metrics["fold"] = fold
        all_results.append(metrics)

    # Aggregate
    res_df = pd.DataFrame(all_results)
    mean_metrics = res_df.drop(columns=["fold"]).mean().to_dict()
    mean_metrics["n_common_genes"] = len(common_genes)
    return mean_metrics


# ---------------------------------------------------------------------------
# Main ablation pipeline
# ---------------------------------------------------------------------------

def run_ablation(
    graph_dir: str = "data/processed/graph",
    out_dir: str = "results/tables",
    embedding_dim: int = 128,
    epochs: int = 150,
) -> str:
    """Run the full edge-type ablation study.

    Args:
        graph_dir: Path to the full heterogeneous graph.
        out_dir: Output directory for results tables.
        embedding_dim: Gene embedding dimension.
        epochs: Max training epochs per condition.

    Returns:
        Path to the ablation results table.
    """
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(os.path.join(graph_dir, "nodes.tsv")):
        log(f"Graph not found at {graph_dir}. Skipping ablation.")
        log("Run build_heterograph.py first to construct the graph.")
        return ""

    nodes, edges = load_graph(graph_dir)
    all_edge_types = get_edge_types(edges)
    log(f"Full graph: {len(nodes)} nodes, {len(edges)} edges")
    log(f"Edge types ({len(all_edge_types)}): {all_edge_types}")
    for etype in all_edge_types:
        n = count_edges_by_type(edges).get(etype, 0)
        log(f"  {etype}: {n} edges")

    # Temp directory for ablated graphs
    tmp_base = os.path.join(out_dir, "..", "ablation_graphs")
    os.makedirs(tmp_base, exist_ok=True)

    # ---- Experiment 1: Full graph (baseline) ----
    log("\n" + "=" * 60)
    log("EXPERIMENT 1: Full Graph (Baseline)")
    log("=" * 60)
    full_graph_dir = os.path.join(tmp_base, "full")
    save_graph(nodes, edges, full_graph_dir)
    try:
        full_metrics = train_gnn_with_graph(full_graph_dir, embedding_dim, epochs, out_dir)
        log(f"  Full graph: MCC={full_metrics.get('MCC', '?'):.4f}  AUROC={full_metrics.get('AUROC', '?'):.4f}")
    except Exception as e:
        log(f"  Full graph training FAILED: {e}")
        full_metrics = {}

    # ---- Experiment 2: Remove each edge type one at a time ----
    log("\n" + "=" * 60)
    log("EXPERIMENT 2: Leave-One-Edge-Type-Out Ablation")
    log("=" * 60)

    ablation_results: List[Dict[str, Any]] = []

    # Add baseline
    ablation_results.append({
        "condition": "full_graph",
        "edge_type_removed": "none",
        "n_edge_types": len(all_edge_types),
        "n_edges_remaining": len(edges),
        **full_metrics,
    })

    for etype in all_edge_types:
        log(f"\n  Removing: {etype}")
        ablated_edges = ablate_edge_type(edges, etype)
        remaining_types = get_edge_types(ablated_edges)
        graph_dir_abl = os.path.join(tmp_base, f"remove_{etype}")
        save_graph(nodes, ablated_edges, graph_dir_abl)

        try:
            metrics = train_gnn_with_graph(graph_dir_abl, embedding_dim, epochs, out_dir)
            delta_mcc = metrics.get("MCC", float("nan")) - full_metrics.get("MCC", float("nan"))
            log(f"    MCC={metrics.get('MCC', '?'):.4f} (Δ={delta_mcc:+.4f})  "
                f"AUROC={metrics.get('AUROC', '?'):.4f}")
        except Exception as e:
            log(f"    FAILED: {e}")
            metrics = {"Accuracy": float("nan"), "MCC": float("nan"), "AUROC": float("nan")}
            delta_mcc = float("nan")

        ablation_results.append({
            "condition": f"remove_{etype}",
            "edge_type_removed": etype,
            "n_edge_types": len(remaining_types),
            "n_edges_remaining": len(ablated_edges),
            "delta_mcc": delta_mcc,
            **metrics,
        })

    # ---- Experiment 3: Each edge type alone ----
    log("\n" + "=" * 60)
    log("EXPERIMENT 3: Single-Edge-Type Graphs")
    log("=" * 60)

    for etype in all_edge_types:
        log(f"\n  Only: {etype}")
        single_edges = keep_only_edge_type(edges, etype)
        graph_dir_single = os.path.join(tmp_base, f"only_{etype}")
        save_graph(nodes, single_edges, graph_dir_single)

        # Only run if we have nodes AND at least some edges
        if len(single_edges) == 0:
            log(f"    No edges for {etype} — skipping")
            ablation_results.append({
                "condition": f"only_{etype}",
                "edge_type_removed": "none",
                "n_edge_types": 0,
                "n_edges_remaining": 0,
                "delta_mcc": float("nan"),
                "Accuracy": float("nan"), "MCC": float("nan"), "AUROC": float("nan"),
            })
            continue

        try:
            metrics = train_gnn_with_graph(graph_dir_single, embedding_dim, epochs, out_dir)
            delta_mcc = metrics.get("MCC", float("nan")) - full_metrics.get("MCC", float("nan"))
            log(f"    MCC={metrics.get('MCC', '?'):.4f} (Δ={delta_mcc:+.4f})  "
                f"AUROC={metrics.get('AUROC', '?'):.4f}")
        except Exception as e:
            log(f"    FAILED: {e}")
            metrics = {"Accuracy": float("nan"), "MCC": float("nan"), "AUROC": float("nan")}
            delta_mcc = float("nan")

        ablation_results.append({
            "condition": f"only_{etype}",
            "edge_type_removed": etype,
            "n_edge_types": 1,
            "n_edges_remaining": len(single_edges),
            "delta_mcc": delta_mcc,
            **metrics,
        })

    # ---- Save results ----
    res_df = pd.DataFrame(ablation_results)

    # Sort: full_graph first, then by MCC descending
    res_df["_sort_mcc"] = res_df["MCC"].fillna(-999)
    res_df = res_df.sort_values("_sort_mcc", ascending=False).drop(columns=["_sort_mcc"])
    res_df = res_df.reset_index(drop=True)

    out_path = os.path.join(out_dir, "ablation_edge_types.tsv")
    res_df.to_csv(out_path, sep="\t", index=False)

    # ---- Print summary ----
    log("\n" + "=" * 60)
    log("ABLATION RESULTS")
    log("=" * 60)

    display_cols = ["condition", "n_edge_types", "n_edges_remaining", "MCC", "AUROC"]
    available_cols = [c for c in display_cols if c in res_df.columns]
    print(res_df[available_cols].to_string(index=False))

    # ---- Analysis ----
    log("\n--- Ablation Analysis ---")
    baseline_mcc = full_metrics.get("MCC", 0)

    # Which edge type causes the biggest drop when removed?
    remove_rows = res_df[res_df["condition"].str.startswith("remove_")]
    if len(remove_rows) > 0 and not pd.isna(baseline_mcc):
        remove_rows = remove_rows.copy()
        remove_rows["mcc_drop"] = baseline_mcc - remove_rows["MCC"]
        remove_rows = remove_rows.sort_values("mcc_drop", ascending=False)
        log("\n  Edge types ranked by MCC drop when REMOVED (largest = most important):")
        for _, r in remove_rows.iterrows():
            direction = "▼" if r["mcc_drop"] > 0.005 else ("▲" if r["mcc_drop"] < -0.005 else "—")
            log(f"    {r['edge_type_removed']:<30} ΔMCC={r['mcc_drop']:+.4f} {direction}")

    # Which single edge type gives best performance?
    single_rows = res_df[res_df["condition"].str.startswith("only_")]
    if len(single_rows) > 0:
        single_rows = single_rows.sort_values("MCC", ascending=False)
        log("\n  Edge types ranked by MCC when used ALONE:")
        for _, r in single_rows.iterrows():
            log(f"    {r['edge_type_removed']:<30} MCC={r['MCC']:.4f}")

    log(f"\n  Results saved to: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Standalone mode (when real graph data not available)
# ---------------------------------------------------------------------------

def run_ablation_synthetic(out_dir: str = "results/tables") -> str:
    """Run a synthetic ablation analysis when real graph data is unavailable.

    Creates a small toy graph, runs the ablation logic, and saves results
    with clear SYNTHETIC_ABLATION markers so they are never confused with
    real experimental results.
    """
    os.makedirs(out_dir, exist_ok=True)

    log("=" * 60)
    log("SYNTHETIC ABLATION (demo mode)")
    log("No real graph data available — using toy graph for code-path validation.")
    log("Results marked SYNTHETIC_ABLATION and excluded from scientific claims.")
    log("=" * 60)

    # Build a small toy graph
    rng = np.random.RandomState(42)
    n_nodes = 200

    nodes = pd.DataFrame({
        "node_id": [f"GENE_{i}" for i in range(n_nodes)],
        "node_type": rng.choice(["gene", "nk_receptor", "tumor_serine_program"], n_nodes),
    })

    edge_types = ["ppi", "ligand_receptor", "tf_target", "metabolic_crosstalk",
                  "sm_topology_axis", "dysfunction_correlation"]
    edges_list = []
    for etype in edge_types:
        n_edges = rng.randint(50, 500)
        for _ in range(n_edges):
            s, d = rng.randint(0, n_nodes - 1, 2)
            if s != d:
                edges_list.append({
                    "src": f"GENE_{s}",
                    "dst": f"GENE_{d}",
                    "edge_type": etype,
                    "weight": round(rng.uniform(0.3, 1.0), 3),
                })
    edges = pd.DataFrame(edges_list).drop_duplicates(subset=["src", "dst", "edge_type"])

    tmp_base = os.path.join(out_dir, "..", "ablation_graphs_SYNTHETIC")
    full_graph_dir = os.path.join(tmp_base, "full")
    save_graph(nodes, edges, full_graph_dir)

    log(f"  Synthetic graph: {len(nodes)} nodes, {len(edges)} edges")
    for etype in edge_types:
        n = len(edges[edges["edge_type"] == etype])
        log(f"    {etype}: {n} edges")

    # Run the ablation on synthetic data
    try:
        result_path = run_ablation(
            graph_dir=full_graph_dir,
            out_dir=out_dir,
            embedding_dim=64,
            epochs=50,
        )
        # Mark as synthetic
        if result_path:
            df = pd.read_csv(result_path, sep="\t")
            df["data_source"] = "SYNTHETIC_ABLATION_DEMO"
            df.to_csv(result_path, sep="\t", index=False)
            log(f"\n  Synthetic ablation results: {result_path}")
        return result_path
    except Exception as e:
        log(f"  Synthetic ablation also failed (expected without full data): {e}")
        # Still output a mock table showing the expected format
        mock_results = pd.DataFrame({
            "condition": ["full_graph"] + [f"remove_{et}" for et in edge_types] + [f"only_{et}" for et in edge_types],
            "edge_type_removed": ["none"] + edge_types + edge_types,
            "n_edge_types": [len(edge_types)] + [len(edge_types) - 1] * len(edge_types) + [1] * len(edge_types),
            "n_edges_remaining": [len(edges)] + [0] * len(edge_types) + [0] * len(edge_types),
            "MCC": [float("nan")] * (1 + 2 * len(edge_types)),
            "AUROC": [float("nan")] * (1 + 2 * len(edge_types)),
            "data_source": ["SYNTHETIC_ABLATION_DEMO"] * (1 + 2 * len(edge_types)),
        })
        mock_path = os.path.join(out_dir, "ablation_edge_types_SYNTHETIC.tsv")
        mock_results.to_csv(mock_path, sep="\t", index=False)
        log(f"  Mock results format: {mock_path}")
        return mock_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    log("=" * 60)
    log("EDGE-TYPE ABLATION STUDY")
    log("=" * 60)

    out_dir = ensure_dir("results/tables")
    graph_dir = "data/processed/graph"

    # Check for real graph data
    if os.path.exists(os.path.join(graph_dir, "nodes.tsv")):
        log(f"Real graph found at {graph_dir}")
        run_ablation(graph_dir=graph_dir, out_dir=out_dir)
    else:
        log(f"No real graph at {graph_dir}")
        log("Running synthetic ablation demo (results marked SYNTHETIC_ABLATION)...")
        run_ablation_synthetic(out_dir)

    log("\nABLATION STUDY COMPLETE!")


if __name__ == "__main__":
    main()
