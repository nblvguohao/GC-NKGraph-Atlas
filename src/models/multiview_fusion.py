"""Auditable lightweight multi-view graph fusion for GC-NKGraph-Atlas.

Each biological edge family is encoded independently with a deterministic
spectral embedding.  The sample-level classifier can either average those
views or learn one global softmax weight per view.  It intentionally does not
implement a Transformer: the graph contains roughly one hundred genes, and the
submission uses the graph as a mechanism-structured diagnostic probe.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from src.immune_scoring.nk_scores import (
    CAF_ECM_TGFB_GENES,
    NK_CYTOTOXICITY_GENES,
    NK_DYSFUNCTION_GENES,
    NK_MARKERS,
)


DEFAULT_VIEW_SPEC: Dict[str, Tuple[str, ...]] = {
    "generic_interaction": ("ppi", "coexpression"),
    "ligand_receptor": ("ligand_receptor",),
    "metabolic_crosstalk": ("metabolic_crosstalk",),
    "sm_topology_axis": ("sm_topology_axis",),
    "go_prior": ("go_prior",),
    "msigdb_prior": ("msigdb_prior",),
}


@dataclass(frozen=True)
class GraphView:
    """One independently encoded graph view over a shared node ordering."""

    name: str
    edge_types: Tuple[str, ...]
    adjacency: np.ndarray
    node_ids: Tuple[str, ...]
    n_edges: int


def build_graph_views(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    view_spec: Mapping[str, Sequence[str]] = DEFAULT_VIEW_SPEC,
    *,
    strict: bool = True,
    require_nonempty: bool = True,
) -> Dict[str, GraphView]:
    """Partition supported edges into independent weighted adjacency views."""

    required_edge_columns = {"src", "dst", "edge_type", "weight"}
    missing_columns = required_edge_columns - set(edges.columns)
    if missing_columns:
        raise ValueError(f"Edges are missing required columns: {sorted(missing_columns)}")
    if "node_id" not in nodes.columns:
        raise ValueError("Nodes are missing required column: node_id")

    assignments: Dict[str, str] = {}
    for view_name, edge_types in view_spec.items():
        for edge_type in edge_types:
            if edge_type in assignments:
                raise ValueError(
                    f"Edge type {edge_type!r} is assigned to both "
                    f"{assignments[edge_type]!r} and {view_name!r}"
                )
            assignments[str(edge_type)] = str(view_name)

    observed = set(edges["edge_type"].astype(str))
    unassigned = observed - set(assignments)
    if strict and unassigned:
        raise ValueError(f"Unassigned edge types: {sorted(unassigned)}")

    node_ids = tuple(nodes["node_id"].astype(str))
    node_to_idx = {node_id: idx for idx, node_id in enumerate(node_ids)}
    views: Dict[str, GraphView] = {}
    for view_name, edge_types in view_spec.items():
        edge_types_tuple = tuple(str(item) for item in edge_types)
        selected = edges[edges["edge_type"].astype(str).isin(edge_types_tuple)]
        if require_nonempty and selected.empty:
            raise ValueError(
                f"Required graph view {view_name!r} is empty; expected edge types "
                f"{list(edge_types_tuple)}"
            )
        adjacency = np.zeros((len(node_ids), len(node_ids)), dtype=np.float32)
        included_pairs = set()
        for row in selected.itertuples(index=False):
            src, dst = str(row.src), str(row.dst)
            if src not in node_to_idx or dst not in node_to_idx or src == dst:
                continue
            i, j = node_to_idx[src], node_to_idx[dst]
            weight = float(row.weight)
            adjacency[i, j] = max(adjacency[i, j], weight)
            adjacency[j, i] = max(adjacency[j, i], weight)
            included_pairs.add(tuple(sorted((i, j))))
        if require_nonempty and not included_pairs:
            raise ValueError(f"Required graph view {view_name!r} has no valid node-pair edges")
        views[str(view_name)] = GraphView(
            name=str(view_name),
            edge_types=edge_types_tuple,
            adjacency=adjacency,
            node_ids=node_ids,
            n_edges=len(included_pairs),
        )
    return views


def compute_spectral_embedding(view: GraphView, embedding_dim: int = 64) -> np.ndarray:
    """Return a deterministic normalized-adjacency SVD embedding."""

    if embedding_dim <= 0:
        raise ValueError("embedding_dim must be positive")
    adjacency = np.asarray(view.adjacency, dtype=np.float64)
    degrees = adjacency.sum(axis=1)
    inv_sqrt = np.zeros_like(degrees)
    nonzero = degrees > 0
    inv_sqrt[nonzero] = 1.0 / np.sqrt(degrees[nonzero])
    normalized = inv_sqrt[:, None] * adjacency * inv_sqrt[None, :]
    u, singular_values, _ = np.linalg.svd(normalized, full_matrices=False)
    width = min(embedding_dim, u.shape[1])
    embedding = u[:, :width] * np.sqrt(np.maximum(singular_values[:width], 0.0))
    if width < embedding_dim:
        embedding = np.pad(embedding, ((0, 0), (0, embedding_dim - width)))
    return embedding.astype(np.float32)


def label_defining_genes() -> Set[str]:
    """Return every gene used to construct any NK immune-state label score."""

    return set(NK_MARKERS) | set(NK_CYTOTOXICITY_GENES) | set(NK_DYSFUNCTION_GENES) | set(CAF_ECM_TGFB_GENES)


def filter_expression_features(
    expression: pd.DataFrame,
    *,
    mode: str,
) -> Tuple[pd.DataFrame, Set[str]]:
    """Apply the pre-registered full or label-masked feature policy."""

    if mode not in {"masked", "full"}:
        raise ValueError("mode must be 'masked' or 'full'")
    if mode == "full":
        return expression.copy(), set()
    removed = set(expression.columns) & label_defining_genes()
    return expression.drop(columns=sorted(removed)).copy(), removed


def fit_standardizer(train: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fit feature scaling parameters on training data only."""

    values = np.asarray(train, dtype=np.float32)
    mean = values.mean(axis=0)
    scale = values.std(axis=0)
    scale = np.where(scale > 1e-8, scale, 1.0)
    return mean.astype(np.float32), scale.astype(np.float32)


def transform_standardized(values: np.ndarray, mean: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """Apply already-fitted training statistics to validation or external data."""

    return ((np.asarray(values, dtype=np.float32) - mean) / scale).astype(np.float32)


def raw_graph_projection(
    expression: pd.DataFrame,
    embedding: np.ndarray,
    node_to_idx: Mapping[str, int],
) -> Tuple[np.ndarray, Tuple[str, ...]]:
    """Project only the expression columns supplied by the caller."""

    used_genes = tuple(str(gene) for gene in expression.columns if str(gene) in node_to_idx)
    if not used_genes:
        raise ValueError("No expression genes overlap the graph nodes")
    values = expression.loc[:, list(used_genes)].to_numpy(dtype=np.float32)
    indices = [node_to_idx[gene] for gene in used_genes]
    projection = values @ np.asarray(embedding, dtype=np.float32)[indices]
    return projection.astype(np.float32), used_genes


def permute_view_node_labels(view: GraphView, rng: np.random.RandomState) -> GraphView:
    """Randomize gene assignment while preserving graph degree/weight distributions.

    Simultaneously permuting adjacency rows and columns retains the complete
    weighted topology and therefore the degree sequence, but breaks its
    assignment to named genes.  This is the appropriate null for asking
    whether pre-registered gene modules align with that topology.
    """

    permutation = rng.permutation(len(view.node_ids))
    adjacency = view.adjacency[np.ix_(permutation, permutation)].copy()
    return GraphView(
        name=view.name,
        edge_types=view.edge_types,
        adjacency=adjacency,
        node_ids=view.node_ids,
        n_edges=view.n_edges,
    )


def _mean_cross_module_cosine(
    embedding: np.ndarray,
    node_ids: Sequence[str],
    source_genes: Set[str],
    target_genes: Set[str],
) -> float:
    node_to_idx = {node_id: idx for idx, node_id in enumerate(node_ids)}
    source_idx = [node_to_idx[gene] for gene in sorted(source_genes) if gene in node_to_idx]
    target_idx = [node_to_idx[gene] for gene in sorted(target_genes) if gene in node_to_idx]
    if not source_idx or not target_idx:
        raise ValueError("Both calibration modules must overlap the graph")
    left = embedding[source_idx]
    right = embedding[target_idx]
    left_norm = left / np.maximum(np.linalg.norm(left, axis=1, keepdims=True), 1e-8)
    right_norm = right / np.maximum(np.linalg.norm(right, axis=1, keepdims=True), 1e-8)
    return float((left_norm @ right_norm.T).mean())


def calibrate_module_coupling(
    view: GraphView,
    *,
    source_genes: Set[str],
    target_genes: Set[str],
    embedding_dim: int = 64,
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Tuple[Dict[str, float | int | str], np.ndarray]:
    """Compare observed module embedding coupling with a node-label null."""

    if n_randomizations <= 0:
        raise ValueError("n_randomizations must be positive")
    observed = _mean_cross_module_cosine(
        compute_spectral_embedding(view, embedding_dim),
        view.node_ids,
        source_genes,
        target_genes,
    )
    rng = np.random.RandomState(seed)
    null_values = np.empty(n_randomizations, dtype=np.float64)
    for index in range(n_randomizations):
        randomized = permute_view_node_labels(view, rng)
        null_values[index] = _mean_cross_module_cosine(
            compute_spectral_embedding(randomized, embedding_dim),
            randomized.node_ids,
            source_genes,
            target_genes,
        )
    empirical_p = (1.0 + float(np.sum(null_values >= observed))) / (n_randomizations + 1.0)
    result: Dict[str, float | int | str] = {
        "view": view.name,
        "observed_coupling": observed,
        "null_mean": float(null_values.mean()),
        "null_std": float(null_values.std(ddof=0)),
        "null_ci_low": float(np.percentile(null_values, 2.5)),
        "null_ci_high": float(np.percentile(null_values, 97.5)),
        "empirical_p": empirical_p,
        "n_randomizations": int(n_randomizations),
        "null_type": "node_label_permutation_degree_sequence_preserved",
    }
    return result, null_values


def _combine_graph_views(views: Mapping[str, GraphView], name: str) -> GraphView:
    first = next(iter(views.values()))
    adjacency = np.zeros_like(first.adjacency)
    edge_types = []
    for view in views.values():
        if view.node_ids != first.node_ids:
            raise ValueError("All graph views must share an identical node ordering")
        adjacency += view.adjacency
        edge_types.extend(view.edge_types)
    return GraphView(
        name=name,
        edge_types=tuple(edge_types),
        adjacency=adjacency,
        node_ids=first.node_ids,
        n_edges=int(np.count_nonzero(np.triu(adjacency, 1))),
    )


def calibrate_view_in_context(
    views: Mapping[str, GraphView],
    *,
    randomized_view_name: str,
    source_genes: Set[str],
    target_genes: Set[str],
    embedding_dim: int = 64,
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Tuple[Dict[str, float | int | str], np.ndarray]:
    """Randomize one view while keeping every other view fixed."""

    if randomized_view_name not in views:
        raise ValueError(f"Unknown view: {randomized_view_name}")
    if n_randomizations <= 0:
        raise ValueError("n_randomizations must be positive")
    observed_view = _combine_graph_views(views, "all_views")
    observed = _mean_cross_module_cosine(
        compute_spectral_embedding(observed_view, embedding_dim),
        observed_view.node_ids,
        source_genes,
        target_genes,
    )
    rng = np.random.RandomState(seed)
    null_values = np.empty(n_randomizations, dtype=np.float64)
    for index in range(n_randomizations):
        randomized_views = dict(views)
        randomized_views[randomized_view_name] = permute_view_node_labels(
            views[randomized_view_name], rng
        )
        combined = _combine_graph_views(randomized_views, "randomized_all_views")
        null_values[index] = _mean_cross_module_cosine(
            compute_spectral_embedding(combined, embedding_dim),
            combined.node_ids,
            source_genes,
            target_genes,
        )
    empirical_p = (1.0 + float(np.sum(null_values >= observed))) / (n_randomizations + 1.0)
    result: Dict[str, float | int | str] = {
        "view": randomized_view_name,
        "observed_coupling": observed,
        "null_mean": float(null_values.mean()),
        "null_std": float(null_values.std(ddof=0)),
        "null_ci_low": float(np.percentile(null_values, 2.5)),
        "null_ci_high": float(np.percentile(null_values, 97.5)),
        "empirical_p": empirical_p,
        "n_randomizations": int(n_randomizations),
        "null_type": "target_view_node_labels_permuted_other_views_fixed",
    }
    return result, null_values


class MultiViewFusionClassifier(nn.Module):
    """Sample classifier with fixed-uniform or learned-global view fusion."""

    def __init__(
        self,
        expr_dim: int,
        emb_dim: int,
        n_views: int,
        *,
        mode: str,
        hidden_dims: Tuple[int, int] = (128, 64),
        dropout: float = 0.5,
        n_classes: int = 2,
    ) -> None:
        super().__init__()
        if mode not in {"uniform", "learned"}:
            raise ValueError("mode must be 'uniform' or 'learned'")
        if n_views <= 0:
            raise ValueError("n_views must be positive")
        self.mode = mode
        self.view_logits = nn.Parameter(torch.zeros(n_views), requires_grad=mode == "learned")
        dims = [expr_dim + emb_dim, *hidden_dims, n_classes]
        layers = []
        for idx in range(len(dims) - 1):
            layers.append(nn.Linear(dims[idx], dims[idx + 1]))
            if idx < len(dims) - 2:
                layers.extend([nn.ReLU(), nn.Dropout(dropout)])
        self.head = nn.Sequential(*layers)

    def view_weights(self) -> torch.Tensor:
        return torch.softmax(self.view_logits, dim=0)

    def fused_projection(self, view_projections: torch.Tensor) -> torch.Tensor:
        if view_projections.ndim != 3:
            raise ValueError("view_projections must have shape (views, samples, embedding_dim)")
        return torch.einsum("v,vnd->nd", self.view_weights(), view_projections)

    def forward(self, expression: torch.Tensor, view_projections: torch.Tensor) -> torch.Tensor:
        return self.head(torch.cat([expression, self.fused_projection(view_projections)], dim=1))
