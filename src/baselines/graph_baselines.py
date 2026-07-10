"""
GC-NKGraph-Atlas: Graph-Method Baselines (Reviewer-requested).

Provides comparison baselines that isolate the contribution of:
  (a) the heterogeneous graph architecture vs simpler graph methods
  (b) mechanism-grounded edges vs data-driven co-expression graph
  (c) graph embeddings vs pure expression features

Baselines:
  1. Node2Vec     — random-walk embeddings on the heterogeneous graph → MLP
  2. GCN           — 2-layer Graph Convolutional Network on the graph
  3. SVD-Graph     — SVD embeddings from the full heterogeneous graph
  4. CoExpression — SVD embeddings from co-expression network (no mechanism edges)
  5. MLP-Expr      — Same MLP architecture, expression features only (no graph)

Usage:
    python src/baselines/graph_baselines.py
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
from src.common.io_utils import ensure_dir
from src.common.seed import set_seed

logger = Logger()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# =========================================================================
# Shared utilities
# =========================================================================

def _load_data():
    """Load training data used by the main GNN model."""
    from src.models.gc_nkgraph_atlas import (_load_training_data,
        select_informative_genes, _build_nx_graph)

    expr, labels = _load_training_data()
    y = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    graph_dir = "data/processed/graph"
    graph = _build_nx_graph(graph_dir)
    graph_genes = set(graph["node_to_idx"].keys())
    selected_genes = select_informative_genes(expr, graph_genes)

    X_expr = expr[selected_genes].values.astype(np.float32)
    return X_expr, y, selected_genes, graph


def _build_full_adj(graph) -> np.ndarray:
    """Build normalized adjacency using all available graph edges."""
    n = graph["num_nodes"]
    adj = np.zeros((n, n), dtype=np.float32)
    for mat in graph["adj"].values():
        adj += mat
    deg = adj.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)
    deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
    adj_norm = deg_inv_sqrt @ adj @ deg_inv_sqrt
    return adj_norm.astype(np.float64)


def _build_coexpression_adj(X_expr, selected_genes, graph,
                            threshold: float = 0.7) -> np.ndarray:
    """Build co-expression adjacency (data-driven baseline, no mechanism edges)."""
    n = graph["num_nodes"]
    gene_to_idx = graph["node_to_idx"]

    adj = np.zeros((n, n), dtype=np.float32)
    graph_genes_in_expr = [g for g in gene_to_idx if g in selected_genes]
    if len(graph_genes_in_expr) < 5:
        return adj

    idx_in_expr = [list(selected_genes).index(g) for g in graph_genes_in_expr]
    expr_subset = X_expr[:, idx_in_expr]
    corr = np.corrcoef(expr_subset.T)
    corr = np.abs(corr)

    for i, gi in enumerate(graph_genes_in_expr):
        for j, gj in enumerate(graph_genes_in_expr):
            if i < j and corr[i, j] > threshold:
                ni, nj = gene_to_idx[gi], gene_to_idx[gj]
                adj[ni, nj] = corr[i, j]
                adj[nj, ni] = corr[i, j]

    deg = adj.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)
    deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
    adj_norm = deg_inv_sqrt @ adj @ deg_inv_sqrt
    return adj_norm.astype(np.float64)


def _evaluate(y_true, y_pred, y_prob=None) -> Dict[str, float]:
    from sklearn.metrics import (
        accuracy_score, balanced_accuracy_score, f1_score,
        matthews_corrcoef, roc_auc_score, average_precision_score)
    m = {
        "Accuracy": float(accuracy_score(y_true, y_pred)),
        "BalancedAccuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "MacroF1": float(f1_score(y_true, y_pred, average="macro")),
        "MCC": float(matthews_corrcoef(y_true, y_pred)),
    }
    if y_prob is not None and y_prob.shape[1] >= 2:
        m["AUROC"] = float(roc_auc_score(y_true, y_prob[:, 1]))
        m["AUPRC"] = float(average_precision_score(y_true, y_prob[:, 1]))
    return m


# =========================================================================
# SVD embedding helper
# =========================================================================

def svd_embed(adj_norm: np.ndarray, embedding_dim: int,
              n_nodes: int) -> np.ndarray:
    """SVD spectral embedding."""
    from scipy.sparse.linalg import svds
    k = min(embedding_dim, n_nodes - 2)
    if k < 2:
        rng = np.random.RandomState(42)
        return rng.randn(n_nodes, embedding_dim).astype(np.float32)

    u, s, vt = svds(adj_norm, k=k)
    idx = np.argsort(s)[::-1]
    emb = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))

    if emb.shape[1] < embedding_dim:
        pad = np.zeros((n_nodes, embedding_dim - emb.shape[1]))
        emb = np.hstack([emb, pad])
    else:
        emb = emb[:, :embedding_dim]
    return emb.astype(np.float32)


# =========================================================================
# Baseline 1: Node2Vec random-walk embeddings
# =========================================================================

class Node2VecBaseline:
    """Node2Vec random-walk embeddings → MLP classifier."""

    def __init__(self, embedding_dim: int = 64, walk_length: int = 30,
                 num_walks: int = 200):
        self.embedding_dim = embedding_dim
        self.walk_length = walk_length
        self.num_walks = num_walks
        self._embeddings: Optional[np.ndarray] = None
        self._gene_to_idx: Dict[str, int] = {}

    def fit(self, graph) -> Node2VecBaseline:
        adj_norm = _build_full_adj(graph)
        self._gene_to_idx = graph["node_to_idx"]
        n_nodes = graph["num_nodes"]

        # Build adjacency list
        adj_list = {i: [] for i in range(n_nodes)}
        rows, cols = np.where(adj_norm > 0)
        for i, j in zip(rows, cols):
            if i != j:
                adj_list[i].append(j)

        # Random walks
        walks = []
        nodes = list(range(n_nodes))
        rng = np.random.RandomState(42)
        for _ in range(self.num_walks):
            rng.shuffle(nodes)
            for start in nodes:
                walk = [start]
                for _ in range(self.walk_length - 1):
                    cur = walk[-1]
                    neighbors = adj_list.get(cur, [])
                    if not neighbors:
                        break
                    walk.append(neighbors[rng.randint(len(neighbors))])
                if len(walk) > 1:
                    walks.append(walk)

        if len(walks) < 10:
            rng = np.random.RandomState(42)
            self._embeddings = rng.randn(n_nodes, self.embedding_dim).astype(np.float32)
            return self

        # Co-occurrence → SVD (DeepWalk approximation)
        window = 5
        cooc = np.zeros((n_nodes, n_nodes), dtype=np.float32)
        for walk in walks:
            for i, node in enumerate(walk):
                for j in range(max(0, i - window), min(len(walk), i + window + 1)):
                    if i != j:
                        cooc[node, walk[j]] += 1.0

        from scipy.sparse.linalg import svds
        k = min(self.embedding_dim, n_nodes - 2)
        u, s, vt = svds(cooc.astype(np.float64), k=k)
        idx = np.argsort(s)[::-1]
        emb = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))
        if emb.shape[1] < self.embedding_dim:
            pad = np.zeros((n_nodes, self.embedding_dim - emb.shape[1]))
            emb = np.hstack([emb, pad])
        else:
            emb = emb[:, :self.embedding_dim]
        self._embeddings = emb.astype(np.float32)
        log(f"  Node2Vec embeddings: {self._embeddings.shape}")
        return self

    def transform(self, gene_list: List[str]) -> np.ndarray:
        X = np.zeros((len(gene_list), self.embedding_dim), dtype=np.float32)
        for i, g in enumerate(gene_list):
            if g in self._gene_to_idx:
                X[i] = self._embeddings[self._gene_to_idx[g]]
        return X


# =========================================================================
# Baseline 2: GCN
# =========================================================================

class SimpleGCNBaseline:
    """2-layer GCN on graph → MLP classifier."""

    def __init__(self, embedding_dim: int = 64):
        self.embedding_dim = embedding_dim
        self._gene_embeddings: Optional[np.ndarray] = None
        self._gene_to_idx: Dict[str, int] = {}

    def fit_embeddings(self, graph, epochs: int = 200) -> SimpleGCNBaseline:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F

        adj_norm = _build_full_adj(graph)
        n_nodes = graph["num_nodes"]
        self._gene_to_idx = graph["node_to_idx"]

        # Node features: degree + log-degree
        adj_combined = np.zeros((n_nodes, n_nodes), dtype=np.float32)
        for mat in graph["adj"].values():
            adj_combined += mat
        degrees = adj_combined.sum(axis=1)
        node_features = np.column_stack([
            degrees, np.log1p(degrees)
        ]).astype(np.float32)
        feat_dim = node_features.shape[1]

        adj_t = torch.from_numpy(adj_norm).float()
        feat_t = torch.from_numpy(node_features)

        class GCN(nn.Module):
            def __init__(self, in_dim, hidden_dim, out_dim):
                super().__init__()
                self.conv1 = nn.Linear(in_dim, hidden_dim)
                self.conv2 = nn.Linear(hidden_dim, out_dim)
            def forward(self, x, adj):
                x = F.relu(self.conv1(adj @ x))
                return self.conv2(adj @ x)

        gcn = GCN(feat_dim, 128, self.embedding_dim)
        optimizer = torch.optim.Adam(gcn.parameters(), lr=1e-3)

        for epoch in range(epochs):
            gcn.train()
            emb = gcn(feat_t, adj_t)
            recon = emb @ emb.T
            loss = F.mse_loss(recon, adj_t)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        gcn.eval()
        with torch.no_grad():
            self._gene_embeddings = gcn(feat_t, adj_t).numpy().astype(np.float32)
        log(f"  GCN embeddings: {self._gene_embeddings.shape}")
        return self

    def transform(self, gene_list: List[str]) -> np.ndarray:
        X = np.zeros((len(gene_list), self.embedding_dim), dtype=np.float32)
        for i, g in enumerate(gene_list):
            if g in self._gene_to_idx:
                X[i] = self._gene_embeddings[self._gene_to_idx[g]]
        return X


# =========================================================================
# Baseline 3: SVD-Graph (full graph SVD)
# =========================================================================

class SVDGraphBaseline:
    """SVD on full heterogeneous graph."""

    def __init__(self, embedding_dim: int = 64):
        self.embedding_dim = embedding_dim
        self._embeddings: Optional[np.ndarray] = None
        self._gene_to_idx: Dict[str, int] = {}

    def fit(self, graph) -> SVDGraphBaseline:
        adj_norm = _build_full_adj(graph)
        self._gene_to_idx = graph["node_to_idx"]
        self._embeddings = svd_embed(adj_norm, self.embedding_dim, graph["num_nodes"])
        log(f"  SVD-Graph embeddings: {self._embeddings.shape}")
        return self

    def transform(self, gene_list: List[str]) -> np.ndarray:
        X = np.zeros((len(gene_list), self.embedding_dim), dtype=np.float32)
        for i, g in enumerate(gene_list):
            if g in self._gene_to_idx:
                X[i] = self._embeddings[self._gene_to_idx[g]]
        return X


# =========================================================================
# Baseline 4: CoExpression graph
# =========================================================================

class CoExpressionBaseline:
    """SVD on co-expression graph (no mechanism edges)."""

    def __init__(self, embedding_dim: int = 64, threshold: float = 0.7):
        self.embedding_dim = embedding_dim
        self.threshold = threshold
        self._embeddings: Optional[np.ndarray] = None
        self._gene_to_idx: Dict[str, int] = {}

    def fit(self, X_expr, selected_genes, graph) -> CoExpressionBaseline:
        adj_norm = _build_coexpression_adj(X_expr, selected_genes, graph, self.threshold)
        self._gene_to_idx = graph["node_to_idx"]
        n_edges = int((adj_norm > 0).sum() / 2)
        self._embeddings = svd_embed(adj_norm, self.embedding_dim, graph["num_nodes"])
        log(f"  CoExpression ({n_edges} edges): embeddings {self._embeddings.shape}")
        return self

    def transform(self, gene_list: List[str]) -> np.ndarray:
        X = np.zeros((len(gene_list), self.embedding_dim), dtype=np.float32)
        for i, g in enumerate(gene_list):
            if g in self._gene_to_idx:
                X[i] = self._embeddings[self._gene_to_idx[g]]
        return X


# =========================================================================
# Shared MLP classifier
# =========================================================================

def _train_mlp_classifier(X_expr, y, gene_embeddings, n_folds=5,
                          hidden_dims=None, dropout=0.6, lr=1.7e-3,
                          wd=5.6e-6, epochs=300, batch_size=16,
                          verbose=True) -> pd.DataFrame:
    """Train MLP classifier with graph-projected features."""
    import torch
    import torch.nn as nn
    from sklearn.model_selection import StratifiedKFold

    hidden_dims = hidden_dims or [256, 128]

    # Prepare features: concat expression + graph projection
    graph_proj = X_expr @ gene_embeddings
    graph_proj = graph_proj / max(np.std(graph_proj), 1e-8)
    X = np.hstack([X_expr, graph_proj]).astype(np.float32)

    # Standardize
    mean = X.mean(axis=0)
    std = X.std(axis=0) + 1e-8
    X = (X - mean) / std

    input_dim = X.shape[1]
    dims = [input_dim] + hidden_dims + [2]

    def _build_fresh_model():
        layers = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < len(dims) - 2:
                layers.append(nn.BatchNorm1d(dims[i + 1]))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(dropout))
        return nn.Sequential(*layers)

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    results = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        rng = np.random.RandomState(42 + fold)
        val_size = max(1, int(len(train_idx) * 0.2))
        val_idx = rng.choice(train_idx, size=val_size, replace=False)
        train_clean = np.setdiff1d(train_idx, val_idx)

        X_t = torch.from_numpy(X[train_clean])
        y_t = torch.from_numpy(y[train_clean].astype(np.int64))
        X_val = torch.from_numpy(X[val_idx])
        y_val = torch.from_numpy(y[val_idx].astype(np.int64))

        n_pos = max((y[train_clean] == 1).sum(), 1)
        n_neg = max((y[train_clean] == 0).sum(), 1)
        total = n_pos + n_neg
        class_weights = torch.tensor(
            [total / (2.0 * n_neg), total / (2.0 * n_pos)], dtype=torch.float32)

        fold_model = _build_fresh_model()
        optimizer = torch.optim.AdamW(fold_model.parameters(), lr=lr, weight_decay=wd)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=20)
        criterion = nn.CrossEntropyLoss(weight=class_weights)

        best_val = 0.0
        best_state = None
        patience = 50
        wait = 0

        for epoch in range(epochs):
            fold_model.train()
            perm = torch.randperm(len(X_t))
            for i in range(0, len(X_t), batch_size):
                idx = perm[i:i + batch_size]
                optimizer.zero_grad()
                loss = criterion(fold_model(X_t[idx]), y_t[idx])
                loss.backward()
                torch.nn.utils.clip_grad_norm_(fold_model.parameters(), 1.0)
                optimizer.step()

            fold_model.eval()
            with torch.no_grad():
                val_acc = (fold_model(X_val).argmax(1) == y_val).float().mean().item()
            if val_acc > best_val:
                best_val = val_acc
                best_state = {k: v.clone() for k, v in fold_model.state_dict().items()}
                wait = 0
            else:
                wait += 1
            scheduler.step(val_acc)
            if wait >= patience:
                break

        if best_state:
            fold_model.load_state_dict(best_state)

        fold_model.eval()
        with torch.no_grad():
            logits = fold_model(torch.from_numpy(X[test_idx]))
            y_pred = logits.argmax(1).numpy()
            y_prob = torch.softmax(logits, 1).numpy()

        metrics = _evaluate(y[test_idx], y_pred, y_prob)
        metrics["fold"] = fold
        results.append(metrics)

    return pd.DataFrame(results)


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Graph-Method Baselines")
    parser.add_argument("--graph-dir", default="data/processed/graph")
    parser.add_argument("--output-dir", default="results/tables")
    parser.add_argument("--embedding-dim", type=int, default=64)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    log("=" * 60)
    log("GC-NKGraph-Atlas: GRAPH-METHOD BASELINES")
    log("=" * 60)

    out_dir = ensure_dir(args.output_dir)
    set_seed(args.seed)

    # --- Load data ---
    log("\n=== Loading data ===")
    X_expr, y, selected_genes, graph = _load_data()
    n_pos, n_neg = (y == 1).sum(), (y == 0).sum()
    log(f"  Samples: {len(y)} (pos={n_pos}, neg={n_neg})")
    log(f"  Features: {X_expr.shape[1]} genes")

    all_results = []

    # --- Baseline 1: Node2Vec ---
    log("\n=== Baseline 1: Node2Vec (random walks) ===")
    n2v = Node2VecBaseline(embedding_dim=args.embedding_dim)
    n2v.fit(graph)
    n2v_emb = n2v.transform(selected_genes)
    res_n2v = _train_mlp_classifier(X_expr, y, n2v_emb, n_folds=args.n_folds,
                                     epochs=args.epochs)
    for col in res_n2v.columns:
        if col != "fold":
            log(f"  {col}: {res_n2v[col].mean():.4f} ± {res_n2v[col].std():.4f}")
    res_n2v["method"] = "Node2Vec"
    all_results.append(res_n2v)

    # --- Baseline 2: GCN ---
    log("\n=== Baseline 2: GCN (2-layer) ===")
    gcn = SimpleGCNBaseline(embedding_dim=args.embedding_dim)
    gcn.fit_embeddings(graph, epochs=200)
    gcn_emb = gcn.transform(selected_genes)
    res_gcn = _train_mlp_classifier(X_expr, y, gcn_emb, n_folds=args.n_folds,
                                     epochs=args.epochs)
    for col in res_gcn.columns:
        if col != "fold":
            log(f"  {col}: {res_gcn[col].mean():.4f} ± {res_gcn[col].std():.4f}")
    res_gcn["method"] = "GCN"
    all_results.append(res_gcn)

    # --- Baseline 3: SVD-Graph ---
    log("\n=== Baseline 3: SVD-Graph (full graph SVD) ===")
    svd_graph = SVDGraphBaseline(embedding_dim=args.embedding_dim)
    svd_graph.fit(graph)
    svd_emb = svd_graph.transform(selected_genes)
    res_svd = _train_mlp_classifier(X_expr, y, svd_emb, n_folds=args.n_folds,
                                     epochs=args.epochs)
    for col in res_svd.columns:
        if col != "fold":
            log(f"  {col}: {res_svd[col].mean():.4f} ± {res_svd[col].std():.4f}")
    res_svd["method"] = "SVD-Graph"
    all_results.append(res_svd)

    # --- Baseline 4: CoExpression ---
    log("\n=== Baseline 4: CoExpression (data-driven, no mechanism edges) ===")
    coexpr = CoExpressionBaseline(embedding_dim=args.embedding_dim)
    coexpr.fit(X_expr, selected_genes, graph)
    coexpr_emb = coexpr.transform(selected_genes)
    res_coexpr = _train_mlp_classifier(X_expr, y, coexpr_emb, n_folds=args.n_folds,
                                        epochs=args.epochs)
    for col in res_coexpr.columns:
        if col != "fold":
            log(f"  {col}: {res_coexpr[col].mean():.4f} ± {res_coexpr[col].std():.4f}")
    res_coexpr["method"] = "CoExpression"
    all_results.append(res_coexpr)

    # --- Baseline 5: MLP-Expr (no graph) ---
    log("\n=== Baseline 5: MLP-Expr (expression only, no graph) ===")
    dummy_emb = np.zeros((len(selected_genes), args.embedding_dim), dtype=np.float32)
    res_mlp = _train_mlp_classifier(X_expr, y, dummy_emb, n_folds=args.n_folds,
                                     epochs=args.epochs)
    for col in res_mlp.columns:
        if col != "fold":
            log(f"  {col}: {res_mlp[col].mean():.4f} ± {res_mlp[col].std():.4f}")
    res_mlp["method"] = "MLP-Expr"
    all_results.append(res_mlp)

    # --- Summary ---
    log(f"\n{'='*60}")
    log("GRAPH BASELINES SUMMARY")
    log(f"{'='*60}")

    summary = pd.concat(all_results, ignore_index=True)
    for method_name in summary["method"].unique():
        subset = summary[summary["method"] == method_name]
        log(f"\n  {method_name}:")
        for col in ["Accuracy", "MCC", "AUROC", "MacroF1"]:
            if col in subset.columns:
                log(f"    {col}: {subset[col].mean():.4f} ± {subset[col].std():.4f}")

    # Save
    out_path = os.path.join(out_dir, "graph_baselines_results.tsv")
    summary.to_csv(out_path, sep="\t", index=False)
    log(f"\nSaved: {out_path}")

    compact = summary.groupby("method")[["Accuracy", "MCC", "AUROC", "MacroF1"]].agg(
        ["mean", "std"]).round(4)
    compact_path = os.path.join(out_dir, "graph_baselines_summary.tsv")
    compact.to_csv(compact_path, sep="\t")
    log(f"Saved: {compact_path}")

    log("\nGRAPH BASELINES COMPLETE!")


if __name__ == "__main__":
    main()
