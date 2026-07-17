"""External-cohort benchmark for lightweight multi-view graph fusion.

The primary analysis removes every gene used to construct the NK-state label.
TCGA-STAD is split into training/validation subsets for early stopping; TCGA-
LIHC is loaded only as an untouched external evaluation cohort.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split

from src.models.multiview_fusion import (
    DEFAULT_VIEW_SPEC,
    GraphView,
    MultiViewFusionClassifier,
    build_graph_views,
    compute_spectral_embedding,
    filter_expression_features,
    fit_standardizer,
    raw_graph_projection,
    transform_standardized,
)


DEFAULT_SEEDS = (1234, 2345, 3456, 4567, 5678, 6789, 7890, 8901, 9012, 1357)
FEATURE_MODES = ("masked", "full")


@dataclass(frozen=True)
class BenchmarkConfig:
    graph_dir: Path
    labels_path: Path
    expression_dir: Path
    output_dir: Path
    seeds: Tuple[int, ...] = DEFAULT_SEEDS
    feature_modes: Tuple[str, ...] = FEATURE_MODES
    embedding_dim: int = 64
    max_epochs: int = 300
    patience: int = 20
    batch_size: int = 32
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    n_bootstrap: int = 2000
    validation_fraction: float = 0.2
    coexpression_edges_path: Path | None = None


@dataclass
class BenchmarkResult:
    per_seed: pd.DataFrame
    summary: pd.DataFrame
    comparisons: pd.DataFrame
    weights: pd.DataFrame
    leave_one_out: pd.DataFrame
    predictions: pd.DataFrame
    verdict: str


@dataclass(frozen=True)
class ModelVariantSpec:
    """One pre-registered benchmark variant."""

    views: Tuple[str, ...]
    fusion_mode: str | None
    representation: str


def model_variant_specs(view_names: Sequence[str]) -> Dict[str, ModelVariantSpec]:
    """Return the complete no-graph, single-view, fusion, and ablation matrix."""

    ordered = tuple(str(name) for name in view_names)
    specs: Dict[str, ModelVariantSpec] = {
        "no_graph": ModelVariantSpec(views=(), fusion_mode=None, representation="expression"),
        "merged_svd": ModelVariantSpec(views=ordered, fusion_mode="merged", representation="graph"),
    }
    for name in ordered:
        specs[f"single__{name}"] = ModelVariantSpec(
            views=(name,), fusion_mode="uniform", representation="graph"
        )
    specs["uniform_multiview"] = ModelVariantSpec(
        views=ordered, fusion_mode="uniform", representation="graph"
    )
    specs["learned_multiview"] = ModelVariantSpec(
        views=ordered, fusion_mode="learned", representation="graph"
    )
    for omitted in ordered:
        specs[f"leave_out__{omitted}"] = ModelVariantSpec(
            views=tuple(name for name in ordered if name != omitted),
            fusion_mode="learned",
            representation="graph",
        )
    return specs


def validate_result_coverage(
    results: pd.DataFrame,
    *,
    feature_modes: Sequence[str],
    variants: Sequence[str],
    seeds: Sequence[int],
) -> None:
    """Fail if any pre-registered benchmark cell is absent or duplicated."""

    required = {
        (str(mode), str(variant), int(seed))
        for mode in feature_modes
        for variant in variants
        for seed in seeds
    }
    observed_counts = results.groupby(["feature_mode", "variant", "seed"]).size().to_dict()
    observed = {(str(mode), str(variant), int(seed)) for mode, variant, seed in observed_counts}
    missing = required - observed
    duplicated = {key for key, count in observed_counts.items() if count != 1}
    if missing or duplicated:
        raise ValueError(
            f"Missing benchmark rows: {sorted(missing)}; duplicated benchmark rows: {sorted(duplicated)}"
        )


def summarize_weight_stability(weights: pd.DataFrame) -> pd.DataFrame:
    """Summarize learned-view weights and enforce the 8/10 rank gate."""

    learned = weights[
        (weights["feature_mode"] == "masked") & (weights["variant"] == "learned_multiview")
    ].copy()
    if learned.empty:
        raise ValueError("Masked learned_multiview weights are required")
    top_rows = learned.loc[learned.groupby("seed")["weight"].idxmax(), ["seed", "view"]]
    top_counts = top_rows["view"].value_counts()
    summary = (
        learned.groupby("view", sort=False)["weight"]
        .agg(weight_mean="mean", weight_std=lambda values: values.std(ddof=0), n_seeds="size")
        .reset_index()
    )
    summary["top_rank_count"] = summary["view"].map(top_counts).fillna(0).astype(int)
    summary["stable_preference"] = summary["top_rank_count"] >= 8
    return summary


def _leave_one_out_bootstrap(
    ensemble_predictions: pd.DataFrame,
    weight_stability: pd.DataFrame,
    *,
    n_bootstrap: int,
) -> pd.DataFrame:
    primary = ensemble_predictions[ensemble_predictions["feature_mode"] == "masked"]
    candidate = primary[primary["variant"] == "learned_multiview"]
    stability = weight_stability.set_index("view")["stable_preference"].to_dict()
    rows = []
    leave_variants = sorted(
        variant for variant in primary["variant"].unique() if variant.startswith("leave_out__")
    )
    for variant in leave_variants:
        omitted = variant.split("__", 1)[1]
        comparator = primary[primary["variant"] == variant]
        aligned = candidate.merge(
            comparator,
            on=["sample_id", "y_true"],
            suffixes=("_candidate", "_leave_out"),
            validate="one_to_one",
        )
        metric_rows = []
        for metric in ("AUROC", "AUPRC"):
            result = paired_stratified_bootstrap(
                aligned["y_true"].to_numpy(),
                aligned["probability_candidate"].to_numpy(),
                aligned["probability_leave_out"].to_numpy(),
                metric=metric.lower(),
                n_bootstrap=n_bootstrap,
                seed=42,
            )
            row = {
                "omitted_view": omitted,
                "metric": metric,
                **result,
                "stable_weight_preference": bool(stability.get(omitted, False)),
            }
            metric_rows.append(row)
            rows.append(row)
        supported = all(row["observed_delta"] > 0 and row["ci_low"] > 0 for row in metric_rows) and bool(
            stability.get(omitted, False)
        )
        for row in rows[-2:]:
            row["view_contribution_verdict"] = (
                "stable_predictive_contribution" if supported else "not_stably_supported"
            )
    return pd.DataFrame(rows)


def _metric_value(y_true: np.ndarray, probabilities: np.ndarray, metric: str) -> float:
    if metric.lower() == "auroc":
        return float(roc_auc_score(y_true, probabilities))
    if metric.lower() == "auprc":
        return float(average_precision_score(y_true, probabilities))
    raise ValueError("metric must be 'auroc' or 'auprc'")


def paired_stratified_bootstrap(
    y_true: np.ndarray,
    candidate_probabilities: np.ndarray,
    comparator_probabilities: np.ndarray,
    *,
    metric: str,
    n_bootstrap: int = 2000,
    seed: int = 42,
) -> Dict[str, float]:
    """Bootstrap a paired metric difference while retaining both classes."""

    y = np.asarray(y_true, dtype=np.int64)
    candidate = np.asarray(candidate_probabilities, dtype=np.float64)
    comparator = np.asarray(comparator_probabilities, dtype=np.float64)
    if not (len(y) == len(candidate) == len(comparator)):
        raise ValueError("y_true and both probability vectors must have equal length")
    class_indices = [np.flatnonzero(y == label) for label in (0, 1)]
    if any(len(indices) == 0 for indices in class_indices):
        raise ValueError("paired stratified bootstrap requires both outcome classes")

    observed = _metric_value(y, candidate, metric) - _metric_value(y, comparator, metric)
    rng = np.random.RandomState(seed)
    deltas = np.empty(n_bootstrap, dtype=np.float64)
    for iteration in range(n_bootstrap):
        sampled = np.concatenate(
            [rng.choice(indices, size=len(indices), replace=True) for indices in class_indices]
        )
        deltas[iteration] = _metric_value(y[sampled], candidate[sampled], metric) - _metric_value(
            y[sampled], comparator[sampled], metric
        )
    ci_low, ci_high = np.percentile(deltas, [2.5, 97.5])
    return {
        "observed_delta": float(observed),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "n_bootstrap": int(n_bootstrap),
    }


def determine_external_gain_verdict(comparisons: pd.DataFrame) -> str:
    """Apply the pre-registered two-metric, three-comparator gain gate."""

    required_comparators = {"no_graph", "merged_svd", "uniform_multiview"}
    required_metrics = {"AUROC", "AUPRC"}
    present_pairs = set(zip(comparisons["comparator"], comparisons["metric"]))
    required_pairs = {(comparator, metric) for comparator in required_comparators for metric in required_metrics}
    if not required_pairs <= present_pairs:
        return "no_stable_external_gain"
    required_rows = comparisons[
        comparisons["comparator"].isin(required_comparators)
        & comparisons["metric"].isin(required_metrics)
    ]
    if bool(((required_rows["observed_delta"] > 0) & (required_rows["ci_low"] > 0)).all()):
        return "stable_external_gain"
    return "no_stable_external_gain"


class _ExpressionClassifier(nn.Module):
    def __init__(self, input_dim: int, dropout: float = 0.5) -> None:
        super().__init__()
        self.head = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 2),
        )

    def forward(self, expression: torch.Tensor) -> torch.Tensor:
        return self.head(expression)


def _set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(1)


def _load_labels(path: Path) -> pd.DataFrame:
    labels = pd.read_csv(path, sep="\t", index_col=0, comment="#")
    labels = labels[labels.index.notna()].copy()
    required = {"dataset", "nk_immune_state"}
    if not required <= set(labels.columns):
        raise ValueError(f"Label table is missing columns: {sorted(required - set(labels.columns))}")
    return labels


def _load_cohort(
    labels: pd.DataFrame,
    expression_dir: Path,
    dataset: str,
) -> Tuple[pd.DataFrame, np.ndarray]:
    filenames = {
        "TCGA-STAD": "tcga_stad_expression.tsv",
        "TCGA-LIHC": "tcga_lihc_expression.tsv",
    }
    path = expression_dir / filenames[dataset]
    if not path.exists():
        raise FileNotFoundError(f"Missing expression matrix: {path}")
    expression = pd.read_csv(path, sep="\t", index_col=0, comment="#")
    cohort_labels = labels[labels["dataset"] == dataset]
    common = expression.index.intersection(cohort_labels.index)
    if len(common) == 0:
        raise ValueError(f"No samples overlap expression and labels for {dataset}")
    expression = expression.loc[common].apply(pd.to_numeric, errors="coerce")
    if expression.isna().any().any():
        raise ValueError(f"Non-numeric or missing expression values in {dataset}")
    y = (cohort_labels.loc[common, "nk_immune_state"] == "NK-hot-cytotoxic").astype(np.int64).to_numpy()
    if len(np.unique(y)) != 2:
        raise ValueError(f"Both outcome classes are required in {dataset}")
    return expression, y


def _load_graph(config: BenchmarkConfig) -> Tuple[pd.DataFrame, pd.DataFrame]:
    nodes_path = config.graph_dir / "nodes.tsv"
    edges_path = config.graph_dir / "edges.tsv"
    if not nodes_path.exists() or not edges_path.exists():
        raise FileNotFoundError(f"Graph directory must contain nodes.tsv and edges.tsv: {config.graph_dir}")
    nodes = pd.read_csv(nodes_path, sep="\t")
    edges = pd.read_csv(edges_path, sep="\t")
    coexpression_path = config.coexpression_edges_path
    if coexpression_path is None:
        candidate = config.graph_dir / "coexpression_edges.tsv"
        coexpression_path = candidate if candidate.exists() else None
    if coexpression_path is not None and coexpression_path.exists() and "coexpression" not in set(edges["edge_type"]):
        coexpression = pd.read_csv(coexpression_path, sep="\t")
        edges = pd.concat([edges, coexpression], ignore_index=True)
    return nodes, edges


def _merged_view(views: Dict[str, GraphView]) -> GraphView:
    first = next(iter(views.values()))
    adjacency = np.zeros_like(first.adjacency)
    for view in views.values():
        adjacency += view.adjacency
    adjacency = np.maximum(adjacency, adjacency.T)
    n_edges = int(np.count_nonzero(np.triu(adjacency, 1)))
    return GraphView(
        name="merged_svd",
        edge_types=tuple(edge_type for view in views.values() for edge_type in view.edge_types),
        adjacency=adjacency,
        node_ids=first.node_ids,
        n_edges=n_edges,
    )


def _prepare_scaled_features(
    train_expression: pd.DataFrame,
    external_expression: pd.DataFrame,
    embeddings: Dict[str, np.ndarray],
    node_to_idx: Dict[str, int],
    train_indices: np.ndarray,
    validation_indices: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    raw_train_expression = train_expression.to_numpy(dtype=np.float32)
    raw_external_expression = external_expression.to_numpy(dtype=np.float32)
    expr_mean, expr_scale = fit_standardizer(raw_train_expression[train_indices])
    expression_split = (
        transform_standardized(raw_train_expression[train_indices], expr_mean, expr_scale),
        transform_standardized(raw_train_expression[validation_indices], expr_mean, expr_scale),
        transform_standardized(raw_external_expression, expr_mean, expr_scale),
    )

    train_projections: Dict[str, np.ndarray] = {}
    validation_projections: Dict[str, np.ndarray] = {}
    external_projections: Dict[str, np.ndarray] = {}
    for name, embedding in embeddings.items():
        raw_train_projection, _ = raw_graph_projection(train_expression, embedding, node_to_idx)
        raw_external_projection, _ = raw_graph_projection(external_expression, embedding, node_to_idx)
        mean, scale = fit_standardizer(raw_train_projection[train_indices])
        train_projections[name] = transform_standardized(raw_train_projection[train_indices], mean, scale)
        validation_projections[name] = transform_standardized(raw_train_projection[validation_indices], mean, scale)
        external_projections[name] = transform_standardized(raw_external_projection, mean, scale)
    return (*expression_split, train_projections, validation_projections, external_projections)


def _stack_projections(projections: Dict[str, np.ndarray], view_names: Sequence[str]) -> np.ndarray:
    return np.stack([projections[name] for name in view_names], axis=0).astype(np.float32)


def _validation_auprc(model: nn.Module, expression: torch.Tensor, projections: torch.Tensor | None, y: np.ndarray) -> float:
    model.eval()
    with torch.no_grad():
        logits = model(expression) if projections is None else model(expression, projections)
        probabilities = torch.softmax(logits, dim=1)[:, 1].cpu().numpy()
    return float(average_precision_score(y, probabilities))


def _train_one(
    model: nn.Module,
    train_expression: np.ndarray,
    validation_expression: np.ndarray,
    y_train: np.ndarray,
    y_validation: np.ndarray,
    *,
    train_projections: np.ndarray | None,
    validation_projections: np.ndarray | None,
    config: BenchmarkConfig,
) -> Tuple[nn.Module, int, float]:
    x_train = torch.from_numpy(train_expression).float()
    x_validation = torch.from_numpy(validation_expression).float()
    p_train = None if train_projections is None else torch.from_numpy(train_projections).float()
    p_validation = None if validation_projections is None else torch.from_numpy(validation_projections).float()
    y_tensor = torch.from_numpy(y_train.astype(np.int64))
    n_pos = max(int((y_train == 1).sum()), 1)
    n_neg = max(int((y_train == 0).sum()), 1)
    total = n_pos + n_neg
    criterion = nn.CrossEntropyLoss(
        weight=torch.tensor([total / (2 * n_neg), total / (2 * n_pos)], dtype=torch.float32)
    )
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay
    )
    best_score = -np.inf
    best_state = copy.deepcopy(model.state_dict())
    best_epoch = 0
    stale_epochs = 0
    for epoch in range(config.max_epochs):
        model.train()
        order = torch.randperm(len(y_train))
        for start in range(0, len(y_train), config.batch_size):
            indices = order[start : start + config.batch_size]
            optimizer.zero_grad()
            logits = (
                model(x_train[indices])
                if p_train is None
                else model(x_train[indices], p_train[:, indices, :])
            )
            loss = criterion(logits, y_tensor[indices])
            loss.backward()
            optimizer.step()
        score = _validation_auprc(model, x_validation, p_validation, y_validation)
        if score > best_score + 1e-8:
            best_score = score
            best_state = copy.deepcopy(model.state_dict())
            best_epoch = epoch + 1
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= config.patience:
                break
    model.load_state_dict(best_state)
    return model, best_epoch, float(best_score)


def _predict(
    model: nn.Module,
    expression: np.ndarray,
    projections: np.ndarray | None,
) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        x = torch.from_numpy(expression).float()
        logits = model(x) if projections is None else model(x, torch.from_numpy(projections).float())
        return torch.softmax(logits, dim=1)[:, 1].cpu().numpy()


def _classification_metrics(y_true: np.ndarray, probabilities: np.ndarray) -> Dict[str, float]:
    predictions = (probabilities >= 0.5).astype(np.int64)
    return {
        "AUROC": float(roc_auc_score(y_true, probabilities)),
        "AUPRC": float(average_precision_score(y_true, probabilities)),
        "MCC": float(matthews_corrcoef(y_true, predictions)),
        "BalancedAccuracy": float(balanced_accuracy_score(y_true, predictions)),
    }


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_external_benchmark(config: BenchmarkConfig) -> BenchmarkResult:
    """Run the complete pre-registered dual-track external benchmark."""

    config.output_dir.mkdir(parents=True, exist_ok=True)
    labels = _load_labels(config.labels_path)
    stad_expression, stad_y = _load_cohort(labels, config.expression_dir, "TCGA-STAD")
    lihc_expression, lihc_y = _load_cohort(labels, config.expression_dir, "TCGA-LIHC")
    nodes, edges = _load_graph(config)
    views = build_graph_views(edges, nodes, DEFAULT_VIEW_SPEC, strict=True, require_nonempty=True)
    view_names = tuple(views)
    variants = model_variant_specs(view_names)
    embeddings = {
        name: compute_spectral_embedding(view, config.embedding_dim) for name, view in views.items()
    }
    embeddings["merged_svd"] = compute_spectral_embedding(_merged_view(views), config.embedding_dim)
    node_ids = tuple(nodes["node_id"].astype(str))
    node_to_idx = {node_id: index for index, node_id in enumerate(node_ids)}

    common_genes = [
        gene for gene in node_ids if gene in stad_expression.columns and gene in lihc_expression.columns
    ]
    if not common_genes:
        raise ValueError("No common graph genes across STAD and LIHC expression matrices")
    stad_expression = stad_expression.loc[:, common_genes]
    lihc_expression = lihc_expression.loc[:, common_genes]

    per_seed_rows = []
    prediction_rows = []
    weight_rows = []
    for feature_mode in config.feature_modes:
        stad_features, removed = filter_expression_features(stad_expression, mode=feature_mode)
        lihc_features, _ = filter_expression_features(lihc_expression, mode=feature_mode)
        if list(stad_features.columns) != list(lihc_features.columns):
            raise ValueError("STAD and LIHC feature columns diverged after masking")
        if stad_features.shape[1] == 0:
            raise ValueError(f"Feature mode {feature_mode!r} removed every graph-expression feature")

        for seed in config.seeds:
            _set_seed(seed)
            train_indices, validation_indices = train_test_split(
                np.arange(len(stad_y)),
                test_size=config.validation_fraction,
                stratify=stad_y,
                random_state=seed,
            )
            (
                train_expression,
                validation_expression,
                external_expression,
                train_projections,
                validation_projections,
                external_projections,
            ) = _prepare_scaled_features(
                stad_features,
                lihc_features,
                embeddings,
                node_to_idx,
                train_indices,
                validation_indices,
            )
            y_train, y_validation = stad_y[train_indices], stad_y[validation_indices]

            for variant_name, variant in variants.items():
                _set_seed(seed)
                selected_views = variant.views
                if variant_name == "no_graph":
                    model: nn.Module = _ExpressionClassifier(train_expression.shape[1])
                    train_stack = validation_stack = external_stack = None
                else:
                    if variant.fusion_mode == "merged":
                        selected_views = ("merged_svd",)
                        fusion_mode = "uniform"
                    else:
                        fusion_mode = str(variant.fusion_mode)
                    train_stack = _stack_projections(train_projections, selected_views)
                    validation_stack = _stack_projections(validation_projections, selected_views)
                    external_stack = _stack_projections(external_projections, selected_views)
                    model = MultiViewFusionClassifier(
                        expr_dim=train_expression.shape[1],
                        emb_dim=config.embedding_dim,
                        n_views=len(selected_views),
                        mode=fusion_mode,
                    )
                model, best_epoch, validation_auprc = _train_one(
                    model,
                    train_expression,
                    validation_expression,
                    y_train,
                    y_validation,
                    train_projections=train_stack,
                    validation_projections=validation_stack,
                    config=config,
                )
                probabilities = _predict(model, external_expression, external_stack)
                metrics = _classification_metrics(lihc_y, probabilities)
                per_seed_rows.append(
                    {
                        "feature_mode": feature_mode,
                        "variant": variant_name,
                        "seed": seed,
                        "n_train": len(train_indices),
                        "n_validation": len(validation_indices),
                        "n_external": len(lihc_y),
                        "n_features": stad_features.shape[1],
                        "n_masked_genes": len(removed),
                        "best_epoch": best_epoch,
                        "validation_AUPRC": validation_auprc,
                        **metrics,
                    }
                )
                for sample_id, y_value, probability in zip(lihc_features.index, lihc_y, probabilities):
                    prediction_rows.append(
                        {
                            "feature_mode": feature_mode,
                            "variant": variant_name,
                            "seed": seed,
                            "sample_id": sample_id,
                            "y_true": int(y_value),
                            "probability": float(probability),
                        }
                    )
                if isinstance(model, MultiViewFusionClassifier):
                    weights = model.view_weights().detach().cpu().numpy()
                    for view_name, weight in zip(selected_views, weights):
                        weight_rows.append(
                            {
                                "feature_mode": feature_mode,
                                "variant": variant_name,
                                "seed": seed,
                                "view": view_name,
                                "weight": float(weight),
                            }
                        )

    per_seed = pd.DataFrame(per_seed_rows)
    validate_result_coverage(
        per_seed,
        feature_modes=config.feature_modes,
        variants=tuple(variants),
        seeds=config.seeds,
    )
    metric_columns = ["AUROC", "AUPRC", "MCC", "BalancedAccuracy"]
    if not np.isfinite(per_seed[metric_columns].to_numpy()).all():
        raise ValueError("Benchmark produced non-finite metrics")
    summary = (
        per_seed.groupby(["feature_mode", "variant"], sort=False)[metric_columns]
        .agg(["mean", lambda values: values.std(ddof=0)])
        .reset_index()
    )
    summary.columns = [
        "_".join(str(part) for part in column if part).replace("<lambda_0>", "std")
        if isinstance(column, tuple)
        else str(column)
        for column in summary.columns
    ]

    predictions = pd.DataFrame(prediction_rows)
    comparisons_rows = []
    primary = predictions[predictions["feature_mode"] == "masked"]
    candidate = (
        primary[primary["variant"] == "learned_multiview"]
        .groupby(["sample_id", "y_true"], sort=False)["probability"]
        .mean()
        .reset_index()
    )
    for comparator_name in ("no_graph", "merged_svd", "uniform_multiview"):
        comparator = (
            primary[primary["variant"] == comparator_name]
            .groupby(["sample_id", "y_true"], sort=False)["probability"]
            .mean()
            .reset_index()
        )
        aligned = candidate.merge(
            comparator,
            on=["sample_id", "y_true"],
            suffixes=("_candidate", "_comparator"),
            validate="one_to_one",
        )
        for metric in ("AUROC", "AUPRC"):
            bootstrap = paired_stratified_bootstrap(
                aligned["y_true"].to_numpy(),
                aligned["probability_candidate"].to_numpy(),
                aligned["probability_comparator"].to_numpy(),
                metric=metric.lower(),
                n_bootstrap=config.n_bootstrap,
                seed=42,
            )
            comparisons_rows.append(
                {"comparator": comparator_name, "metric": metric, **bootstrap}
            )
    comparisons = pd.DataFrame(comparisons_rows)
    verdict = determine_external_gain_verdict(comparisons)
    comparisons["verdict"] = verdict
    weights = pd.DataFrame(weight_rows)
    leave_one_out = per_seed[per_seed["variant"].str.startswith("leave_out__")].copy()
    ensemble_predictions = (
        predictions.groupby(["feature_mode", "variant", "sample_id", "y_true"], sort=False)[
            "probability"
        ]
        .mean()
        .reset_index()
    )

    weight_stability = summarize_weight_stability(weights)
    leave_one_out_bootstrap = _leave_one_out_bootstrap(
        ensemble_predictions, weight_stability, n_bootstrap=config.n_bootstrap
    )
    output_tables = {
        "multiview_benchmark_per_seed.tsv": per_seed,
        "multiview_benchmark_summary.tsv": summary,
        "multiview_bootstrap_comparisons.tsv": comparisons,
        "multiview_fusion_weights.tsv": weights,
        "multiview_weight_stability.tsv": weight_stability,
        "multiview_leave_one_out.tsv": leave_one_out,
        "multiview_leave_one_out_bootstrap.tsv": leave_one_out_bootstrap,
        "multiview_external_predictions.tsv": predictions,
        "multiview_external_predictions_ensemble.tsv": ensemble_predictions,
    }
    for filename, table in output_tables.items():
        table.to_csv(config.output_dir / filename, sep="\t", index=False)

    provenance = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            **asdict(config),
            "graph_dir": str(config.graph_dir),
            "labels_path": str(config.labels_path),
            "expression_dir": str(config.expression_dir),
            "output_dir": str(config.output_dir),
            "coexpression_edges_path": (
                str(config.coexpression_edges_path) if config.coexpression_edges_path else None
            ),
        },
        "input_sha256": {
            "nodes": _file_sha256(config.graph_dir / "nodes.tsv"),
            "edges": _file_sha256(config.graph_dir / "edges.tsv"),
            "labels": _file_sha256(config.labels_path),
            "stad_expression": _file_sha256(config.expression_dir / "tcga_stad_expression.tsv"),
            "lihc_expression": _file_sha256(config.expression_dir / "tcga_lihc_expression.tsv"),
        },
        "view_spec": {name: list(edge_types) for name, edge_types in DEFAULT_VIEW_SPEC.items()},
        "verdict": verdict,
    }
    (config.output_dir / "multiview_benchmark_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return BenchmarkResult(per_seed, summary, comparisons, weights, leave_one_out, predictions, verdict)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--expression-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--coexpression-edges", type=Path)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    parser.add_argument("--feature-modes", nargs="+", choices=FEATURE_MODES, default=list(FEATURE_MODES))
    parser.add_argument("--max-epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=20)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    result = run_external_benchmark(
        BenchmarkConfig(
            graph_dir=args.graph_dir,
            labels_path=args.labels,
            expression_dir=args.expression_dir,
            output_dir=args.output_dir,
            coexpression_edges_path=args.coexpression_edges,
            seeds=tuple(args.seeds),
            feature_modes=tuple(args.feature_modes),
            max_epochs=args.max_epochs,
            patience=args.patience,
            n_bootstrap=args.n_bootstrap,
        )
    )
    print(f"Multiview benchmark complete: verdict={result.verdict}")


if __name__ == "__main__":
    main()
