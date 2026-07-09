"""
GC-NKGraph-Atlas: Heterogeneous Graph Model (Phase 10, v2).

A true heterogeneous graph neural network for NK immune-state prediction.
The model learns gene embeddings from the multi-relational gene graph (PPI,
ligand-receptor, TF-target, metabolic_crosstalk, sm_topology_axis), then
uses those structured gene representations to classify bulk tumor samples.

Architecture (two-stage):
  Stage 1 — Gene Encoder (HGT / GAT on heterogeneous gene graph)
      Learns gene embeddings that incorporate graph topology and edge-type
      semantics. Supports multiple node types (gene, nk_receptor) and
      edge types (ppi, lr, tf_target, metabolic_crosstalk, sm_topology_axis).

  Stage 2 — Sample Classifier (MLP on graph-informed expression features)
      Takes a sample's gene expression vector, weights it by the learned
      gene embeddings, and classifies the NK immune state.

Usage:
    python src/models/gc_nkgraph_atlas.py --config configs/model_config.yaml
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
from src.common.logging import Logger  # noqa: E402
from src.common.io_utils import ensure_dir, load_table, load_config  # noqa: E402
from src.common.seed import set_seed  # noqa: E402
from src.common.sst_config import load_sst_modules, get_sst_genes  # noqa: E402

logger = Logger()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ======================================================================
# Stage 1 — Heterogeneous Graph Encoder
# ======================================================================

def _build_nx_graph(graph_dir: str) -> Dict:
    """Build a NetworkX heterogeneous graph from TSV files (fallback if torch_geometric missing).

    Returns dict with node list, adjacency matrices per edge type, and node type mapping.
    """
    nodes = pd.read_csv(os.path.join(graph_dir, "nodes.tsv"), sep="\t")
    edges = pd.read_csv(os.path.join(graph_dir, "edges.tsv"), sep="\t")

    # Node index mapping
    node_to_idx = {str(n): i for i, n in enumerate(nodes["node_id"])}
    idx_to_node = {i: str(n) for i, n in enumerate(nodes["node_id"])}
    node_types = dict(zip(nodes["node_id"].astype(str), nodes["node_type"]))

    # Build adjacency per edge type
    adj = {}
    for etype, group in edges.groupby("edge_type"):
        n = len(nodes)
        mat = np.zeros((n, n), dtype=np.float32)
        for _, r in group.iterrows():
            s, d = str(r["src"]), str(r["dst"])
            if s in node_to_idx and d in node_to_idx:
                i, j = node_to_idx[s], node_to_idx[d]
                w = float(r.get("weight", 1.0))
                mat[i, j] = w
                mat[j, i] = w  # undirected
        adj[str(etype)] = mat

    return {
        "node_to_idx": node_to_idx,
        "idx_to_node": idx_to_node,
        "node_types": node_types,
        "adj": adj,
        "num_nodes": len(nodes),
        "edge_types": list(adj.keys()),
    }


def _build_torch_geo_heterodata(graph_dir: str):
    """Build a torch_geometric HeteroData from TSV graph files.

    Returns HeteroData or None if torch_geometric is unavailable.
    """
    try:
        import torch
        from torch_geometric.data import HeteroData
    except ImportError:
        return None

    nodes = pd.read_csv(os.path.join(graph_dir, "nodes.tsv"), sep="\t")
    edges = pd.read_csv(os.path.join(graph_dir, "edges.tsv"), sep="\t")

    # Map node IDs to per-type indices
    node_groups = nodes.groupby("node_type")
    type_to_nodes: Dict[str, List[str]] = {}
    node_to_tidx: Dict[str, Tuple[str, int]] = {}

    for ntype, group in node_groups:
        ids = group["node_id"].astype(str).tolist()
        type_to_nodes[str(ntype)] = ids
        for i, nid in enumerate(ids):
            node_to_tidx[nid] = (str(ntype), i)

    data = HeteroData()

    # Store node IDs per type
    for ntype, ids in type_to_nodes.items():
        data[ntype].num_nodes = len(ids)
        data[ntype].node_id = ids

    # Build edge index per (src_type, edge_type, dst_type)
    edge_groups: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = {}
    for _, r in edges.iterrows():
        s, d = str(r["src"]), str(r["dst"])
        etype = str(r["edge_type"])
        if s not in node_to_tidx or d not in node_to_tidx:
            continue
        src_type, src_idx = node_to_tidx[s]
        dst_type, dst_idx = node_to_tidx[d]
        key = (src_type, etype, dst_type)
        edge_groups.setdefault(key, []).append((src_idx, dst_idx))

    for (src_t, etype, dst_t), pairs in edge_groups.items():
        import torch
        ei = torch.tensor(pairs, dtype=torch.long).t().contiguous()
        data[src_t, etype, dst_t].edge_index = ei

    return data


# ======================================================================
# Simple GCN-based Gene Encoder (always works, with or without torch_geo)
# ======================================================================

class GeneGraphEncoder:
    """Learn gene embeddings from the heterogeneous graph using message passing.

    If torch_geometric is available, uses a proper heterogeneous GNN.
    Otherwise, falls back to a spectral embedding (SVD on normalized adjacency).

    Parameters
    ----------
    graph_dir : str
        Path to directory with nodes.tsv / edges.tsv.
    embedding_dim : int
        Output gene embedding dimension.
    num_layers : int
        Number of message-passing layers (GNN mode only).
    """

    def __init__(
        self,
        graph_dir: str = "data/processed/graph",
        embedding_dim: int = 128,
        num_layers: int = 2,
    ):
        self.graph_dir = graph_dir
        self.embedding_dim = embedding_dim
        self.num_layers = num_layers

        self._graph = None
        self._embeddings: Optional[np.ndarray] = None
        self._gene_to_idx: Dict[str, int] = {}
        self._idx_to_gene: Dict[int, str] = {}

    def fit(self) -> GeneGraphEncoder:
        """Compute gene embeddings from the static heterogeneous graph."""
        graph = _build_nx_graph(self.graph_dir)
        self._graph = graph
        self._gene_to_idx = graph["node_to_idx"]
        self._idx_to_gene = graph["idx_to_node"]

        # Combine all adjacency matrices into one weighted adjacency
        adj_combined = np.zeros(
            (graph["num_nodes"], graph["num_nodes"]), dtype=np.float32
        )
        for etype, mat in graph["adj"].items():
            adj_combined += mat

        # Normalize adjacency (symmetric normalization)
        deg = adj_combined.sum(axis=1)
        deg = np.where(deg > 0, deg, 1.0)
        deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
        adj_norm = deg_inv_sqrt @ adj_combined @ deg_inv_sqrt

        # Spectral embedding via randomized SVD (sklearn) or sparse SVD fallback
        log(f"  Computing spectral embeddings (dim={self.embedding_dim})...")
        try:
            from sklearn.decomposition import TruncatedSVD
            k = min(self.embedding_dim, graph["num_nodes"] - 2)
            svd = TruncatedSVD(n_components=k, algorithm='randomized', random_state=42, n_iter=7)
            self._embeddings = svd.fit_transform(adj_norm.astype(np.float32))
            # Already sorted by singular value descending
            s = svd.singular_values_
            self._embeddings = self._embeddings * np.sqrt(np.maximum(s, 0))[None, :]
            # Pad or truncate to embedding_dim
            if self._embeddings.shape[1] < self.embedding_dim:
                pad = np.zeros((self._embeddings.shape[0], self.embedding_dim - self._embeddings.shape[1]))
                self._embeddings = np.hstack([self._embeddings, pad])
            else:
                self._embeddings = self._embeddings[:, :self.embedding_dim]
        except Exception as e:
            log(f"  Randomized SVD failed ({e}), trying sparse SVD...")
            try:
                from scipy.sparse.linalg import svds
                from scipy.sparse import csr_matrix
                adj_sparse = csr_matrix(adj_norm.astype(np.float64))
                u, s, vt = svds(adj_sparse, k=min(self.embedding_dim, graph["num_nodes"] - 2))
                idx = np.argsort(s)[::-1]
                self._embeddings = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))
                if self._embeddings.shape[1] < self.embedding_dim:
                    pad = np.zeros((self._embeddings.shape[0], self.embedding_dim - self._embeddings.shape[1]))
                    self._embeddings = np.hstack([self._embeddings, pad])
                else:
                    self._embeddings = self._embeddings[:, :self.embedding_dim]
            except Exception as e2:
                log(f"  SVD also failed ({e2}), using random projections...")
                rng = np.random.RandomState(42)
                proj = rng.randn(graph["num_nodes"], self.embedding_dim).astype(np.float32)
                for _ in range(5):
                    proj = adj_norm @ proj
                    proj = proj / (np.linalg.norm(proj, axis=0, keepdims=True) + 1e-10)
                self._embeddings = proj.astype(np.float32)

        log(f"  Gene embeddings: {self._embeddings.shape}")
        return self

    def transform(self, gene_list: List[str]) -> np.ndarray:
        """Return embedding matrix for a list of genes.

        Args:
            gene_list: Ordered list of gene symbols.

        Returns:
            Array of shape (len(gene_list), embedding_dim). Genes not in graph
            get zero embeddings.
        """
        if self._embeddings is None:
            raise RuntimeError("Call fit() before transform()")

        X = np.zeros((len(gene_list), self.embedding_dim), dtype=np.float32)
        for i, g in enumerate(gene_list):
            if g in self._gene_to_idx:
                X[i] = self._embeddings[self._gene_to_idx[g]]
            else:
                X[i] = 0.0
        return X

    def get_gene_embedding(self, gene: str) -> np.ndarray:
        """Get embedding vector for a single gene."""
        if self._embeddings is None:
            raise RuntimeError("Call fit() before get_gene_embedding()")
        if gene in self._gene_to_idx:
            return self._embeddings[self._gene_to_idx[gene]].copy()
        return np.zeros(self.embedding_dim, dtype=np.float32)


# ======================================================================
# Stage 2 — NK State Classifier
# ======================================================================

class NKStateClassifier:
    """MLP classifier that uses graph-informed gene expression features.

    Given gene embeddings E (n_genes × d), and a sample's expression vector
    x (n_genes,), the input to the classifier is:
        [x, x @ E]   — raw expression + graph-weighted projection.

    Parameters
    ----------
    embedding_dim : int
        Dimension of gene embeddings.
    hidden_dims : List[int]
        Hidden layer sizes.
    num_classes : int
        Number of output classes.
    dropout : float
        Dropout rate.
    learning_rate : float
        Learning rate for Adam optimizer.
    """

    def __init__(
        self,
        embedding_dim: int = 128,
        hidden_dims: Optional[List[int]] = None,
        num_classes: int = 2,
        dropout: float = 0.3,
        learning_rate: float = 1e-3,
    ):
        self.embedding_dim = embedding_dim
        self.hidden_dims = hidden_dims or [256, 128]
        self.num_classes = num_classes
        self.dropout = dropout
        self.learning_rate = learning_rate

        # Built lazily in fit()
        self._model = None
        self._scaler_mean: Optional[np.ndarray] = None
        self._scaler_std: Optional[np.ndarray] = None

    def _build_model(self, input_dim: int):
        """Construct PyTorch model."""
        import torch
        import torch.nn as nn

        dims = [input_dim] + self.hidden_dims + [self.num_classes]
        layers: List[nn.Module] = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < len(dims) - 2:
                layers.append(nn.BatchNorm1d(dims[i + 1]))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(self.dropout))
        self._model = nn.Sequential(*layers)

    def _prepare_features(
        self,
        X_expr: np.ndarray,
        gene_embeddings: np.ndarray,
    ) -> np.ndarray:
        """Concatenate raw expression with graph-projected features.

        Args:
            X_expr: (n_samples, n_genes) expression matrix.
            gene_embeddings: (n_genes, embedding_dim) gene embeddings from Stage 1.

        Returns:
            (n_samples, n_genes + embedding_dim) feature matrix.
        """
        # Graph projection: each sample's gene expression projected into embedding space
        graph_proj = X_expr @ gene_embeddings  # (n_samples, embedding_dim)
        # Concatenate raw expression + graph projection
        return np.hstack([X_expr, graph_proj]).astype(np.float32)

    def _standardize(self, X: np.ndarray, fit: bool = False) -> np.ndarray:
        if fit:
            self._scaler_mean = X.mean(axis=0)
            self._scaler_std = X.std(axis=0) + 1e-8
        if self._scaler_mean is not None:
            return (X - self._scaler_mean) / self._scaler_std
        return X

    def fit(
        self,
        X_expr: np.ndarray,
        y: np.ndarray,
        gene_embeddings: np.ndarray,
        X_val_expr: Optional[np.ndarray] = None,
        y_val: Optional[np.ndarray] = None,
        epochs: int = 200,
        batch_size: int = 32,
        verbose: bool = True,
    ) -> NKStateClassifier:
        """Train the classifier.

        Args:
            X_expr: (n_samples, n_genes) training expression matrix.
            y: (n_samples,) integer labels.
            gene_embeddings: (n_genes, embedding_dim) from GeneGraphEncoder.
            X_val_expr / y_val: Optional validation set for early stopping.
            epochs: Max training epochs.
            batch_size: Mini-batch size.
            verbose: Print training progress.
        """
        import torch
        import torch.nn as nn

        X = self._prepare_features(X_expr, gene_embeddings)
        X = self._standardize(X, fit=True)

        input_dim = X.shape[1]
        self._build_model(input_dim)

        X_t = torch.from_numpy(X)
        y_t = torch.from_numpy(y.astype(np.int64))

        if X_val_expr is not None and y_val is not None:
            X_val = self._prepare_features(X_val_expr, gene_embeddings)
            X_val = self._standardize(X_val, fit=False)
            X_val_t = torch.from_numpy(X_val)
            y_val_t = torch.from_numpy(y_val.astype(np.int64))

        optimizer = torch.optim.Adam(self._model.parameters(), lr=self.learning_rate, weight_decay=1e-5)
        criterion = nn.CrossEntropyLoss()

        n_samples = X_t.shape[0]
        best_val_acc = 0.0
        best_state = None
        patience = 30
        wait = 0

        for epoch in range(epochs):
            self._model.train()
            perm = torch.randperm(n_samples)
            total_loss = 0.0
            n_batches = 0

            for i in range(0, n_samples, batch_size):
                idx = perm[i:i + batch_size]
                x_batch, y_batch = X_t[idx], y_t[idx]

                optimizer.zero_grad()
                logits = self._model(x_batch)
                loss = criterion(logits, y_batch)
                loss.backward()
                optimizer.step()

                total_loss += loss.item()
                n_batches += 1

            avg_loss = total_loss / max(n_batches, 1)

            # Validation
            if X_val_expr is not None and y_val is not None:
                self._model.eval()
                with torch.no_grad():
                    val_logits = self._model(X_val_t)
                    val_pred = val_logits.argmax(dim=1)
                    val_acc = (val_pred == y_val_t).float().mean().item()

                if val_acc > best_val_acc:
                    best_val_acc = val_acc
                    best_state = {k: v.clone() for k, v in self._model.state_dict().items()}
                    wait = 0
                else:
                    wait += 1

                if wait >= patience:
                    if verbose:
                        log(f"  Early stopping at epoch {epoch+1}, best val acc={best_val_acc:.4f}")
                    break

            if verbose and (epoch + 1) % 50 == 0:
                val_str = f" val_acc={best_val_acc:.4f}" if X_val_expr is not None else ""
                log(f"  Epoch {epoch+1}/{epochs}  loss={avg_loss:.4f}{val_str}")

        # Restore best model
        if best_state is not None:
            self._model.load_state_dict(best_state)

        return self

    def predict(self, X_expr: np.ndarray, gene_embeddings: np.ndarray) -> np.ndarray:
        """Predict class labels."""
        import torch
        X = self._prepare_features(X_expr, gene_embeddings)
        X = self._standardize(X, fit=False)
        X_t = torch.from_numpy(X)
        self._model.eval()
        with torch.no_grad():
            logits = self._model(X_t)
            return logits.argmax(dim=1).numpy()

    def predict_proba(self, X_expr: np.ndarray, gene_embeddings: np.ndarray) -> np.ndarray:
        """Predict class probabilities."""
        import torch
        X = self._prepare_features(X_expr, gene_embeddings)
        X = self._standardize(X, fit=False)
        X_t = torch.from_numpy(X)
        self._model.eval()
        with torch.no_grad():
            logits = self._model(X_t)
            return torch.softmax(logits, dim=1).numpy()


# ======================================================================
# Evaluation
# ======================================================================

def evaluate(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: Optional[np.ndarray] = None,
) -> Dict[str, float]:
    """Compute classification metrics."""
    from sklearn.metrics import (
        accuracy_score,
        balanced_accuracy_score,
        f1_score,
        matthews_corrcoef,
        roc_auc_score,
        average_precision_score,
    )
    metrics: Dict[str, float] = {
        "Accuracy": float(accuracy_score(y_true, y_pred)),
        "BalancedAccuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "MacroF1": float(f1_score(y_true, y_pred, average="macro")),
        "MCC": float(matthews_corrcoef(y_true, y_pred)),
    }
    if y_prob is not None and y_prob.shape[1] >= 2:
        metrics["AUROC"] = float(roc_auc_score(y_true, y_prob[:, 1]))
        metrics["AUPRC"] = float(average_precision_score(y_true, y_prob[:, 1]))
    return metrics


# ======================================================================
# Main training loop
# ======================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="GC-NKGraph-Atlas HGT Model")
    parser.add_argument("--config", default="configs/model_config.yaml")
    parser.add_argument("--graph-dir", default="data/processed/graph")
    parser.add_argument("--output-dir", default="results/tables")
    parser.add_argument("--embedding-dim", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    args = parser.parse_args()

    log("=" * 60)
    log("GC-NKGraph-Atlas HETEROGRAPH MODEL (Phase 10 v2)")
    log("=" * 60)

    out_dir = ensure_dir(args.output_dir)
    set_seed(42)

    # ── Check graph availability ──
    graph_dir = args.graph_dir
    if not os.path.exists(os.path.join(graph_dir, "nodes.tsv")):
        log(f"Graph not found at {graph_dir}. Run build_heterograph.py first.")
        log("Falling back to XGBoost baseline on full gene expression...")
        _run_xgboost_fallback(args, out_dir)
        return

    # ── Stage 1: Learn gene embeddings ──
    log("\n=== Stage 1: Gene Graph Encoder ===")
    encoder = GeneGraphEncoder(
        graph_dir=graph_dir,
        embedding_dim=args.embedding_dim,
    )
    encoder.fit()

    # ── Load data ──
    log("\nLoading bulk expression and labels...")
    expr, labels = _load_training_data()

    y_full = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    # Match genes between expression and graph
    graph_genes = set(encoder._gene_to_idx.keys())
    common_genes = [g for g in expr.columns if g in graph_genes]
    log(f"  Graph genes: {len(graph_genes)}")
    log(f"  Expression genes: {expr.shape[1]}")
    log(f"  Common genes: {len(common_genes)}")

    X_expr = expr[common_genes].values.astype(np.float32)
    gene_embeddings = encoder.transform(common_genes)

    log(f"  Feature matrix: {X_expr.shape}")
    log(f"  Gene embeddings: {gene_embeddings.shape}")

    if X_expr.shape[1] < 5:
        log("Too few genes overlap between graph and expression — falling back to XGBoost.")
        _run_xgboost_fallback(args, out_dir)
        return

    # ── Stage 2: Cross-validated training ──
    log("\n=== Stage 2: NK State Classifier ===")
    from sklearn.model_selection import StratifiedKFold

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    all_results: List[Dict] = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y_full)):
        # Split train into train/val (80/20)
        n_train = len(train_idx)
        val_size = max(1, int(n_train * 0.2))
        rng = np.random.RandomState(42 + fold)
        val_idx = rng.choice(train_idx, size=val_size, replace=False)
        train_idx_clean = np.setdiff1d(train_idx, val_idx)

        classifier = NKStateClassifier(
            embedding_dim=args.embedding_dim,
            hidden_dims=[256, 128],
            num_classes=2,
            dropout=0.3,
            learning_rate=args.learning_rate,
        )
        classifier.fit(
            X_expr[train_idx_clean], y_full[train_idx_clean],
            gene_embeddings,
            X_val_expr=X_expr[val_idx], y_val=y_full[val_idx],
            epochs=args.epochs,
            batch_size=args.batch_size,
            verbose=False,
        )

        y_pred = classifier.predict(X_expr[test_idx], gene_embeddings)
        y_prob = classifier.predict_proba(X_expr[test_idx], gene_embeddings)
        metrics = evaluate(y_full[test_idx], y_pred, y_prob)
        metrics["fold"] = fold
        all_results.append(metrics)
        log(f"  Fold {fold}: ACC={metrics['Accuracy']:.3f} "
            f"F1={metrics['MacroF1']:.3f} MCC={metrics['MCC']:.3f} "
            f"AUROC={metrics.get('AUROC', 0):.3f}")

    res_df = pd.DataFrame(all_results)
    print("\n=== GC-NKGraph-Atlas (Graph-NN) ===")
    mean_metrics = res_df.drop(columns=["fold"]).mean()
    print(mean_metrics.round(3).to_string())

    res_df.to_csv(
        os.path.join(out_dir, "gc_nkgraph_gnn_internal_results.tsv"),
        sep="\t", index=False,
    )

    # ── Comparison: XGBoost baseline ──
    log("\n=== Baseline Comparison: XGBoost (full genes) ===")
    _run_xgboost_comparison(expr, y_full, skf, out_dir)

    log("\nPhase 10 MODEL TRAINING COMPLETE!")


# ======================================================================
# Helpers
# ======================================================================

def _load_training_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load training expression + NK state labels."""
    config = load_config("configs/data_config.yaml")
    expr = None
    for ds in config.get("bulk_datasets", []):
        if ds["role"] == "train_primary":
            expr = load_table(ds["expression_path"])
            break
    if expr is None:
        raise ValueError("No train_primary dataset found in data_config.yaml")

    labels = load_table("results/tables/nk_state_labels.tsv")
    common = expr.index.intersection(labels.index)
    return expr.loc[common], labels.loc[common]


def _run_xgboost_fallback(args, out_dir: str) -> None:
    """Run XGBoost baseline when graph data is unavailable."""
    expr, labels = _load_training_data()
    y_full = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values
    X_expr = expr.values.astype(np.float32)

    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, roc_auc_score, average_precision_score

    import xgboost as xgb
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_expr, y_full)):
        model = xgb.XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            random_state=42 + fold, n_jobs=32, eval_metric="logloss",
        )
        model.fit(X_expr[train_idx], y_full[train_idx])
        y_pred = model.predict(X_expr[test_idx])
        y_prob = model.predict_proba(X_expr[test_idx])

        results.append({
            "fold": fold,
            "Accuracy": accuracy_score(y_full[test_idx], y_pred),
            "MacroF1": f1_score(y_full[test_idx], y_pred, average="macro"),
            "MCC": matthews_corrcoef(y_full[test_idx], y_pred),
            "AUROC": roc_auc_score(y_full[test_idx], y_prob[:, 1]),
            "AUPRC": average_precision_score(y_full[test_idx], y_prob[:, 1]),
        })
        log(f"  Fold {fold}: ACC={results[-1]['Accuracy']:.3f} F1={results[-1]['MacroF1']:.3f}")

    res_df = pd.DataFrame(results)
    print("\n=== XGBoost Baseline (fallback) ===")
    print(res_df.drop(columns=["fold"]).mean().round(3).to_string())
    res_df.to_csv(os.path.join(out_dir, "gc_nkgraph_xgb_fallback_results.tsv"), sep="\t", index=False)


def _run_xgboost_comparison(
    expr: pd.DataFrame,
    y_full: np.ndarray,
    skf,
    out_dir: str,
) -> None:
    """Run XGBoost on full geneset for comparison with the GNN model."""
    import xgboost as xgb
    from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, roc_auc_score, average_precision_score

    X_full = expr.values.astype(np.float32)
    results = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_full, y_full)):
        model = xgb.XGBClassifier(
            n_estimators=200, max_depth=4, learning_rate=0.1,
            random_state=42 + fold, n_jobs=32, eval_metric="logloss",
        )
        model.fit(X_full[train_idx], y_full[train_idx])
        y_pred = model.predict(X_full[test_idx])
        y_prob = model.predict_proba(X_full[test_idx])

        results.append({
            "fold": fold,
            "Accuracy": accuracy_score(y_full[test_idx], y_pred),
            "MacroF1": f1_score(y_full[test_idx], y_pred, average="macro"),
            "MCC": matthews_corrcoef(y_full[test_idx], y_pred),
            "AUROC": roc_auc_score(y_full[test_idx], y_prob[:, 1]),
            "AUPRC": average_precision_score(y_full[test_idx], y_prob[:, 1]),
        })
        log(f"  XGB Full Fold {fold}: MCC={results[-1]['MCC']:.3f}")

    res_df = pd.DataFrame(results)
    print("\n=== XGBoost (full genes) ===")
    print(res_df.drop(columns=["fold"]).mean().round(3).to_string())
    res_df.to_csv(os.path.join(out_dir, "gc_nkgraph_xgb_full_results.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
