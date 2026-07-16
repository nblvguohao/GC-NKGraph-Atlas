"""
T19 -- Multi-view edge-type fusion with a learnable attention weight
(exploratory, GRAFT/TREE-inspired).

Context: T17 and T18 both assess an edge type's value by hard ablation --
build two graph variants (with/without the edge type), embed each
separately, train two classifiers, bootstrap the MCC difference. That
harness has two weaknesses this script is designed to avoid:
  1. Merging multiple edge types into one adjacency matrix forces an
     arbitrary tie-break when two edge types connect the same node pair
     (see the max-weight fix applied to build_adj_matrix in T17 after this
     was found to change results). Multi-view fusion sidesteps this
     entirely -- each edge type keeps its own adjacency matrix and never
     competes for the same matrix cell.
  2. A single ablation run answers "does dropping this edge type change
     held-out MCC" with only a bootstrap over the eval set as uncertainty
     quantification -- no visibility into the edge type's *relative*
     contribution when several are present together.

Design (lightweight adaptation of GRAFT's multi-view graph encoding, Eq 3-5
of Cho & Cho, Brief Bioinform 2026, bbaf706 -- see manuscript S1.4/S4.4):
  - Each edge type (ppi, mechanism [metabolic_crosstalk + sm_topology_axis],
    ligand_receptor, go_prior, msigdb_prior) gets its own adjacency matrix
    and its own spectral (SVD) embedding -- same encoder as every other
    script in this repo, NOT a full GCN per view (this graph has 104 nodes;
    GRAFT/TREE's multi-layer-GCN-per-view design targets 10,000-20,000-node
    genome-scale networks and would almost certainly overfit here).
  - A single learnable softmax weight vector over views (analogous to
    GRAFT's alpha, Eq 4) fuses the per-view graph projections before
    concatenation with raw expression, trained end-to-end with the
    classifier rather than fixed post-hoc.
  - After training, the learned view weights ARE the edge-type-importance
    answer for that run -- no separate multi-variant ablation needed.

This is still evaluated with the same STAD-train / LIHC-test cross-cohort
design as T17/T18, over more seeds (5, vs T18's 3) since GRAFT/TREE's
10-fold-CV convention is not directly portable to a single train/test
cohort split.

Output: results/tables/t19_multiview_fusion_weights.tsv (learned view
weights, one row per seed) and results/tables/t19_multiview_fusion_perf.tsv
(cross-cohort MCC/AUROC per seed, for comparison against T18's variants).

Run: python src/graph_construction/build_heterograph.py \
        --out-dir data/processed/graph_ablation/multiview \
        --enable-go-prior --enable-msigdb-prior
     python src/a100_recompute/run_t19_multiview_fusion.py
"""
import sys, time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir
from src.common.seed import set_seed
from src.a100_recompute.run_t17_edge_external_value import (
    compute_svd_embeddings, load_cohort_expression,
)

GRAPH_DIR = "data/processed/graph_ablation/multiview"
SEEDS = [1234, 2345, 3456, 4567, 5678]
RESULTS = "results/tables"
OUT_WEIGHTS = f"{RESULTS}/t19_multiview_fusion_weights.tsv"
OUT_PERF = f"{RESULTS}/t19_multiview_fusion_perf.tsv"

VIEW_EDGE_TYPES: Dict[str, List[str]] = {
    "ppi": ["ppi"],
    "mechanism": ["metabolic_crosstalk", "sm_topology_axis"],
    "ligand_receptor": ["ligand_receptor"],
    "go_prior": ["go_prior"],
    "msigdb_prior": ["msigdb_prior"],
}
VIEW_NAMES = list(VIEW_EDGE_TYPES.keys())


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_view_adjacency(edges_df: pd.DataFrame, nodes_df: pd.DataFrame,
                          edge_types: List[str]) -> Tuple[np.ndarray, Dict[str, int]]:
    node_ids = nodes_df["node_id"].astype(str).tolist()
    n = len(node_ids)
    node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    adj = np.zeros((n, n), dtype=np.float32)
    sub = edges_df[edges_df["edge_type"].isin(edge_types)]
    for _, row in sub.iterrows():
        src, dst = str(row["src"]), str(row["dst"])
        if src in node_to_idx and dst in node_to_idx:
            i, j = node_to_idx[src], node_to_idx[dst]
            w = float(row.get("weight", 1.0))
            if w > adj[i, j]:
                adj[i, j] = w
                adj[j, i] = w
    return adj, node_to_idx


class MultiViewFusionClassifier(nn.Module):
    """GRAFT-style learnable softmax fusion over K per-edge-type graph views.

    view_weight (analogous to GRAFT's `a` in softmax(Z_i a)) is a single
    learnable vector over views, shared across all samples/genes -- a
    deliberately simplified, global version of GRAFT's per-node attention,
    appropriate for a 104-node graph and a sample-level (not gene-level)
    classification task.
    """

    def __init__(self, expr_dim: int, emb_dim: int, n_views: int,
                 hidden_dims=(128, 64), n_classes: int = 2):
        super().__init__()
        self.view_weight = nn.Parameter(torch.zeros(n_views))  # softmax(0)=uniform init
        dims = [expr_dim + emb_dim] + list(hidden_dims) + [n_classes]
        layers = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < len(dims) - 2:
                layers.append(nn.BatchNorm1d(dims[i + 1]))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(0.5))
        self.head = nn.Sequential(*layers)

    def fused_projection(self, view_projs: torch.Tensor) -> torch.Tensor:
        """view_projs: (K, N, D) -> fused (N, D) via softmax(view_weight)."""
        w = torch.softmax(self.view_weight, dim=0)  # (K,)
        return torch.einsum("k,knd->nd", w, view_projs)

    def forward(self, x_expr: torch.Tensor, view_projs: torch.Tensor) -> torch.Tensor:
        fused = self.fused_projection(view_projs)
        x = torch.cat([x_expr, fused], dim=1)
        return self.head(x)


def raw_graph_projection(expr_df: pd.DataFrame, emb: np.ndarray, node_to_idx: Dict[str, int]) -> np.ndarray:
    genes_in_graph = [g for g in expr_df.columns if g in node_to_idx]
    X_expr = expr_df[genes_in_graph].values.astype(np.float32)
    gene_indices = [node_to_idx[g] for g in genes_in_graph]
    E = emb[gene_indices]
    proj = X_expr @ E
    proj = proj / max(np.std(proj), 1e-8)
    return proj.astype(np.float32)


def main() -> None:
    log("=" * 70)
    log("T19 -- Multi-view edge-type fusion with learnable attention (exploratory)")
    log("=" * 70)

    nodes = pd.read_csv(f"{GRAPH_DIR}/nodes.tsv", sep="\t")
    edges = pd.read_csv(f"{GRAPH_DIR}/edges.tsv", sep="\t")
    log(f"Graph: {len(nodes)} nodes, {len(edges)} edges across {edges['edge_type'].nunique()} types")

    labels = pd.read_csv("results/tables/nk_state_labels.tsv", sep="\t", index_col=0, comment="#")
    labels = labels[labels.index.notnull()]

    expr_train, y_train = load_cohort_expression(labels, "TCGA-STAD")
    expr_test, y_test = load_cohort_expression(labels, "TCGA-LIHC")
    log(f"STAD train: {expr_train.shape[0]} samples, {y_train.sum()} pos")
    log(f"LIHC test:  {expr_test.shape[0]} samples, {y_test.sum()} pos")

    log("\nComputing per-view spectral embeddings...")
    embeddings: Dict[str, np.ndarray] = {}
    node_to_idx = None
    for view, etypes in VIEW_EDGE_TYPES.items():
        adj, n2i = build_view_adjacency(edges, nodes, etypes)
        node_to_idx = n2i  # identical node set/order across views
        n_edges = int((adj > 0).sum() / 2)
        emb = compute_svd_embeddings(adj, embedding_dim=64)
        embeddings[view] = emb
        log(f"  {view:16s} ({','.join(etypes)}): {n_edges} edges")

    genes_in_graph = [g for g in expr_train.columns if g in node_to_idx]
    log(f"\n{len(genes_in_graph)} genes overlap between graph and expression data")

    # Per-view raw (unstandardized) projections, shared across seeds since they
    # don't depend on classifier weights.
    train_projs_raw = np.stack([raw_graph_projection(expr_train, embeddings[v], node_to_idx)
                                 for v in VIEW_NAMES], axis=0)  # (K, N_train, D)
    test_projs_raw = np.stack([raw_graph_projection(expr_test, embeddings[v], node_to_idx)
                                for v in VIEW_NAMES], axis=0)  # (K, N_test, D)

    X_train_expr = expr_train[genes_in_graph].values.astype(np.float32)
    X_test_expr = expr_test[genes_in_graph].values.astype(np.float32)

    weight_rows, perf_rows = [], []

    for seed in SEEDS:
        set_seed(seed)
        # Standardize raw expression using train statistics only.
        mean_e, std_e = X_train_expr.mean(axis=0), X_train_expr.std(axis=0) + 1e-8
        Xtr_e = (X_train_expr - mean_e) / std_e
        Xte_e = (X_test_expr - mean_e) / std_e

        model = MultiViewFusionClassifier(
            expr_dim=Xtr_e.shape[1], emb_dim=64, n_views=len(VIEW_NAMES),
        )
        n_pos, n_neg = int(y_train.sum()), int((1 - y_train).sum())
        total = n_pos + n_neg
        class_weights = torch.tensor([total / (2.0 * n_neg), total / (2.0 * n_pos)], dtype=torch.float32)
        criterion = nn.CrossEntropyLoss(weight=class_weights)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

        Xtr_e_t = torch.from_numpy(Xtr_e).float()
        train_projs_t = torch.from_numpy(train_projs_raw).float()
        y_train_t = torch.from_numpy(y_train.astype(np.int64))

        model.train()
        n_samples = Xtr_e_t.shape[0]
        batch_size = 32
        for epoch in range(200):
            perm = torch.randperm(n_samples)
            for start in range(0, n_samples, batch_size):
                idx = perm[start:start + batch_size]
                optimizer.zero_grad()
                logits = model(Xtr_e_t[idx], train_projs_t[:, idx, :])
                loss = criterion(logits, y_train_t[idx])
                loss.backward()
                optimizer.step()

        model.eval()
        with torch.no_grad():
            train_logits = model(Xtr_e_t, train_projs_t)
            train_preds = train_logits.argmax(axis=1).numpy()
            train_mcc = matthews_corrcoef(y_train, train_preds)

            Xte_e_t = torch.from_numpy(Xte_e).float()
            test_projs_t = torch.from_numpy(test_projs_raw).float()
            test_logits = model(Xte_e_t, test_projs_t)
            test_probs = torch.softmax(test_logits, dim=1).numpy()[:, 1]
            test_preds = test_logits.argmax(axis=1).numpy()
            test_mcc = matthews_corrcoef(y_test, test_preds)
            test_auroc = roc_auc_score(y_test, test_probs)

            view_w = torch.softmax(model.view_weight, dim=0).numpy()

        log(f"  seed={seed}  Train MCC={train_mcc:.4f}  Test MCC={test_mcc:.4f}  "
            f"AUROC={test_auroc:.4f}  weights={dict(zip(VIEW_NAMES, np.round(view_w, 3)))}")

        perf_rows.append({"seed": seed, "train_mcc": round(train_mcc, 4),
                           "test_mcc": round(test_mcc, 4), "test_auroc": round(test_auroc, 4)})
        weight_rows.append({"seed": seed, **{v: round(float(w), 4) for v, w in zip(VIEW_NAMES, view_w)}})

    ensure_dir(RESULTS)
    perf_df = pd.DataFrame(perf_rows)
    weights_df = pd.DataFrame(weight_rows)
    perf_df.to_csv(OUT_PERF, sep="\t", index=False)
    weights_df.to_csv(OUT_WEIGHTS, sep="\t", index=False)

    log(f"\nWritten {len(perf_df)} rows to {OUT_PERF}")
    log(f"Written {len(weights_df)} rows to {OUT_WEIGHTS}")

    log("\n" + "=" * 70)
    log("Aggregate across seeds:")
    log(f"  Test MCC:   {perf_df['test_mcc'].mean():.4f} +/- {perf_df['test_mcc'].std():.4f}")
    log(f"  Test AUROC: {perf_df['test_auroc'].mean():.4f} +/- {perf_df['test_auroc'].std():.4f}")
    log("  Mean learned view weight (softmax, higher = more relied upon):")
    for v in VIEW_NAMES:
        log(f"    {v:16s} {weights_df[v].mean():.4f} +/- {weights_df[v].std():.4f}")
    log("=" * 70)
    log("T19 COMPLETE (exploratory -- not part of the reported manuscript graph)")


if __name__ == "__main__":
    main()
