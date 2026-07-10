"""
GC-NKGraph-Atlas: Bayesian Hyperparameter Optimization (Optuna).

Tunes the GNN model (GeneGraphEncoder + NKStateClassifier) using
Tree-structured Parzen Estimator (TPE) sampling.

Search space:
  - embedding_dim: 32-256
  - hidden_dims: architectural choices
  - dropout: 0.1-0.7
  - learning_rate: 1e-5 - 1e-2 (log-uniform)
  - weight_decay: 1e-6 - 1e-1 (log-uniform)
  - batch_size: 8-64
  - edge_weights: tuned per edge type

Usage:
    python src/models/hyperparam_search.py --n-trials 100 --n-ensemble 5
    python src/models/hyperparam_search.py --study-name my_study --n-trials 50
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ.setdefault("OPTUNA_WARN_EXPERIMENTAL", "0")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.log_utils import Logger
from src.common.io_utils import ensure_dir, load_table, load_config
from src.common.seed import set_seed
from src.common.sst_config import load_sst_modules, get_sst_genes

# Import model components
from src.models.gc_nkgraph_atlas import (
    GeneGraphEncoder,
    NKStateClassifier,
    select_informative_genes,
    evaluate,
    _load_training_data,
    _build_nx_graph,
)

import optuna
from optuna.trial import TrialState

logger = Logger()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# =========================================================================
# Trial-space definitions
# =========================================================================

# Hidden layer architecture candidates
HIDDEN_ARCHITECTURES = [
    [128, 64],          # current default (compact)
    [256, 128],         # wider
    [64, 32],           # narrower
    [128, 64, 32],      # deeper-narrow
    [256, 128, 64],     # deeper-wide
    [128, 128],         # equal-width
    [256, 256, 128],    # very wide
    [64, 64],           # narrow-equal
]

# Edge type weight configurations (relative importance per edge type)
EDGE_WEIGHT_CONFIGS = [
    # (ppi, metabolic_crosstalk, sm_topology_axis, ligand_receptor, dysfunction_correlation, tf_target)
    [1.0, 2.0, 2.0, 1.5, 1.0, 0.5],   # default (current)
    [1.0, 3.0, 3.0, 2.0, 1.0, 0.5],   # sst-heavy (更重视 SST 轴)
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],   # uniform (所有边权重相等)
    [1.5, 2.0, 2.0, 1.0, 0.5, 0.5],   # ppi-up
    [1.0, 2.0, 2.0, 2.0, 1.5, 1.0],   # interaction-heavy (更重视配体-受体)
    [2.0, 2.0, 2.0, 1.0, 0.5, 0.3],   # core-only (PPI + SST 为核心)
]

# Edge type order matching _build_nx_graph output
EDGE_TYPES_ORDER = ["ppi", "metabolic_crosstalk", "sm_topology_axis",
                    "ligand_receptor", "dysfunction_correlation", "tf_target"]


# =========================================================================
# Bayesian search objective
# =========================================================================

@dataclass
class SearchContext:
    """Shared context across trials to avoid redundant loading."""
    X_expr: np.ndarray
    y_full: np.ndarray
    graph_dir: str
    selected_genes: List[str]
    n_folds: int = 3        # fewer folds for speed during search
    n_ensemble: int = 3     # fewer ensemble models for speed


def _build_weighted_adj(graph_dir: str,
                        edge_weights: Dict[str, float]) -> Tuple[np.ndarray, int]:
    """Build weighted adjacency matrix for gene graph with given edge weights."""
    graph = _build_nx_graph(graph_dir)
    n_nodes = graph["num_nodes"]
    adj_combined = np.zeros((n_nodes, n_nodes), dtype=np.float32)

    for etype, mat in graph["adj"].items():
        w = edge_weights.get(etype, 1.0)
        adj_combined += w * mat

    # Symmetric normalized Laplacian
    deg = adj_combined.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)
    deg_inv_sqrt = np.diag(1.0 / np.sqrt(deg))
    adj_norm = deg_inv_sqrt @ adj_combined @ deg_inv_sqrt

    adj_norm = adj_norm.astype(np.float64)
    return adj_norm, n_nodes


def svd_embed(adj_norm: np.ndarray, embedding_dim: int,
              n_nodes: int) -> np.ndarray:
    """SVD spectral embedding for given adjacency matrix."""
    from scipy.sparse.linalg import svds

    k = min(embedding_dim, n_nodes - 2)
    if k < 2:
        # fallback: random embeddings
        rng = np.random.RandomState(42)
        return rng.randn(n_nodes, embedding_dim).astype(np.float32)

    u, s, vt = svds(adj_norm, k=k)
    idx = np.argsort(s)[::-1]
    emb = u[:, idx] @ np.diag(np.sqrt(np.maximum(s[idx], 0)))

    if emb.shape[1] < embedding_dim:
        pad = np.zeros((emb.shape[0], embedding_dim - emb.shape[1]),
                       dtype=np.float64)
        emb = np.hstack([emb, pad])
    else:
        emb = emb[:, :embedding_dim]

    return emb.astype(np.float32)


def trial_objective(trial: optuna.Trial, ctx: SearchContext) -> float:
    """Optuna objective: train with suggested params, return mean MCC."""

    # --- Sample hyperparameters ---
    embedding_dim = trial.suggest_categorical("embedding_dim", [32, 48, 64, 96, 128, 192, 256])

    # Architecture
    arch_idx = trial.suggest_categorical("hidden_architecture",
                                         list(range(len(HIDDEN_ARCHITECTURES))))
    hidden_dims = HIDDEN_ARCHITECTURES[arch_idx]

    dropout = trial.suggest_float("dropout", 0.1, 0.7, step=0.1)

    learning_rate = trial.suggest_float("learning_rate", 1e-5, 1e-2, log=True)
    weight_decay = trial.suggest_float("weight_decay", 1e-6, 1e-1, log=True)
    batch_size = trial.suggest_categorical("batch_size", [8, 16, 32, 64])

    # Edge weight config
    edge_wt_idx = trial.suggest_categorical("edge_weight_config",
                                            list(range(len(EDGE_WEIGHT_CONFIGS))))
    edge_wts = EDGE_WEIGHT_CONFIGS[edge_wt_idx]
    edge_weight_dict = dict(zip(EDGE_TYPES_ORDER, edge_wts))

    # --- Build gene embeddings with tuned edge weights ---
    try:
        adj_norm, n_nodes = _build_weighted_adj(ctx.graph_dir, edge_weight_dict)
        gene_embeddings_full = svd_embed(adj_norm, embedding_dim, n_nodes)

        # Map to selected genes
        graph = _build_nx_graph(ctx.graph_dir)
        gene_to_idx = graph["node_to_idx"]
        gene_embeddings = np.zeros((len(ctx.selected_genes), embedding_dim),
                                   dtype=np.float32)
        for i, g in enumerate(ctx.selected_genes):
            if g in gene_to_idx:
                gene_embeddings[i] = gene_embeddings_full[gene_to_idx[g]]
    except Exception as e:
        log(f"  Trial {trial.number}: embedding failed ({e}), pruned")
        raise optuna.exceptions.TrialPruned()

    # --- Cross-validation ---
    from sklearn.model_selection import StratifiedKFold

    skf = StratifiedKFold(n_splits=ctx.n_folds, shuffle=True, random_state=42)
    ensemble_seeds = [42, 123, 456][:ctx.n_ensemble]

    all_mccs = []

    for ensemble_seed in ensemble_seeds:
        for fold, (train_idx, test_idx) in enumerate(
            skf.split(ctx.X_expr, ctx.y_full)
        ):
            n_train = len(train_idx)
            val_size = max(1, int(n_train * 0.2))
            rng = np.random.RandomState(ensemble_seed + fold)
            val_idx = rng.choice(train_idx, size=val_size, replace=False)
            train_clean = np.setdiff1d(train_idx, val_idx)

            classifier = NKStateClassifier(
                input_dim=ctx.X_expr.shape[1],
                embedding_dim=embedding_dim,
                hidden_dims=hidden_dims,
                num_classes=2,
                dropout=dropout,
                learning_rate=learning_rate,
                weight_decay=weight_decay,
            )

            try:
                classifier.fit(
                    ctx.X_expr[train_clean], ctx.y_full[train_clean],
                    gene_embeddings,
                    X_val_expr=ctx.X_expr[val_idx], y_val=ctx.y_full[val_idx],
                    epochs=200,  # fewer epochs during search
                    batch_size=batch_size,
                    verbose=False,
                )

                y_pred = classifier.predict(ctx.X_expr[test_idx], gene_embeddings)
                y_prob = classifier.predict_proba(ctx.X_expr[test_idx], gene_embeddings)
                metrics = evaluate(ctx.y_full[test_idx], y_pred, y_prob)
                all_mccs.append(metrics["MCC"])
            except Exception as e:
                log(f"  Trial {trial.number}: training failed ({e}), pruned")
                raise optuna.exceptions.TrialPruned()

    mean_mcc = float(np.mean(all_mccs))

    # --- Report intermediate values for pruning ---
    trial.report(mean_mcc, step=0)

    # Prune if significantly below median
    # (MedianPruner will handle this automatically)

    return mean_mcc


# =========================================================================
# Full evaluation with best parameters
# =========================================================================

def evaluate_best_params(
    ctx: SearchContext,
    best_params: Dict,
    output_dir: str,
    n_folds: int = 5,
    n_ensemble: int = 5,
) -> pd.DataFrame:
    """Full 5-fold CV evaluation using the best-found hyperparameters."""

    embedding_dim = best_params["embedding_dim"]
    hidden_dims = HIDDEN_ARCHITECTURES[best_params["hidden_architecture"]]
    dropout = best_params["dropout"]
    learning_rate = best_params["learning_rate"]
    weight_decay = best_params["weight_decay"]
    batch_size = best_params["batch_size"]
    edge_wts = EDGE_WEIGHT_CONFIGS[best_params["edge_weight_config"]]
    edge_weight_dict = dict(zip(EDGE_TYPES_ORDER, edge_wts))

    log(f"\n{'='*60}")
    log("FINAL EVALUATION — Best Parameters")
    log(f"{'='*60}")
    log(f"  embedding_dim:    {embedding_dim}")
    log(f"  hidden_dims:      {hidden_dims}")
    log(f"  dropout:          {dropout}")
    log(f"  learning_rate:    {learning_rate:.2e}")
    log(f"  weight_decay:     {weight_decay:.2e}")
    log(f"  batch_size:       {batch_size}")
    log(f"  edge_weights:     {edge_weight_dict}")
    log(f"  n_folds:          {n_folds}")
    log(f"  n_ensemble:       {n_ensemble}")

    # Build embeddings
    adj_norm, n_nodes = _build_weighted_adj(ctx.graph_dir, edge_weight_dict)
    gene_embeddings_full = svd_embed(adj_norm, embedding_dim, n_nodes)
    graph = _build_nx_graph(ctx.graph_dir)
    gene_to_idx = graph["node_to_idx"]
    gene_embeddings = np.zeros((len(ctx.selected_genes), embedding_dim),
                               dtype=np.float32)
    for i, g in enumerate(ctx.selected_genes):
        if g in gene_to_idx:
            gene_embeddings[i] = gene_embeddings_full[gene_to_idx[g]]

    # Cross-validate
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    ensemble_seeds = [42, 123, 456, 789, 1024][:n_ensemble]

    all_results = []

    for seed_idx, ensemble_seed in enumerate(ensemble_seeds):
        log(f"\n--- Ensemble Model {seed_idx+1}/{n_ensemble} (seed={ensemble_seed}) ---")

        for fold, (train_idx, test_idx) in enumerate(skf.split(ctx.X_expr, ctx.y_full)):
            n_train = len(train_idx)
            val_size = max(1, int(n_train * 0.2))
            rng = np.random.RandomState(ensemble_seed + fold)
            val_idx = rng.choice(train_idx, size=val_size, replace=False)
            train_clean = np.setdiff1d(train_idx, val_idx)

            classifier = NKStateClassifier(
                input_dim=ctx.X_expr.shape[1],
                embedding_dim=embedding_dim,
                hidden_dims=hidden_dims,
                num_classes=2,
                dropout=dropout,
                learning_rate=learning_rate,
                weight_decay=weight_decay,
            )

            classifier.fit(
                ctx.X_expr[train_clean], ctx.y_full[train_clean],
                gene_embeddings,
                X_val_expr=ctx.X_expr[val_idx], y_val=ctx.y_full[val_idx],
                epochs=300,
                batch_size=batch_size,
                verbose=(fold == 0 and seed_idx == 0),
            )

            y_pred = classifier.predict(ctx.X_expr[test_idx], gene_embeddings)
            y_prob = classifier.predict_proba(ctx.X_expr[test_idx], gene_embeddings)
            metrics = evaluate(ctx.y_full[test_idx], y_pred, y_prob)
            metrics["fold"] = fold
            metrics["ensemble_seed"] = ensemble_seed
            all_results.append(metrics)

            log(f"  Fold {fold}: ACC={metrics['Accuracy']:.3f} "
                f"F1={metrics['MacroF1']:.3f} MCC={metrics['MCC']:.3f} "
                f"AUROC={metrics.get('AUROC', 0):.3f}")

    # Summarize
    res_df = pd.DataFrame(all_results)
    print(f"\n{'='*60}")
    print(" BAYESIAN-OPTIMIZED GNN RESULTS")
    print(f"{'='*60}")

    mean_metrics = res_df.drop(columns=["fold", "ensemble_seed"]).mean()
    std_metrics = res_df.drop(columns=["fold", "ensemble_seed"]).std()
    for k in mean_metrics.index:
        print(f"  {k:<20} {mean_metrics[k]:.4f} ± {std_metrics[k]:.4f}")

    # Save
    res_path = os.path.join(output_dir, "gc_nkgraph_gnn_bayesian_results.tsv")
    res_df.to_csv(res_path, sep="\t", index=False)
    log(f"\nResults saved: {res_path}")

    # Also save best params
    params_path = os.path.join(output_dir, "gc_nkgraph_best_hyperparams.tsv")
    params_df = pd.DataFrame([{
        "embedding_dim": embedding_dim,
        "hidden_dims": str(hidden_dims),
        "dropout": dropout,
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
        "batch_size": batch_size,
        "edge_weights": str(edge_weight_dict),
        "best_mcc": mean_metrics["MCC"],
        "best_auroc": mean_metrics["AUROC"],
        "best_f1": mean_metrics["MacroF1"],
    }])
    params_df.to_csv(params_path, sep="\t", index=False)
    log(f"Best params saved: {params_path}")

    return res_df


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bayesian Hyperparameter Optimization for GC-NKGraph-Atlas GNN"
    )
    parser.add_argument("--graph-dir", default="data/processed/graph")
    parser.add_argument("--output-dir", default="results/tables")
    parser.add_argument("--n-trials", type=int, default=100,
                       help="Number of Optuna trials (default: 100)")
    parser.add_argument("--n-folds", type=int, default=3,
                       help="CV folds during search (default: 3, for speed)")
    parser.add_argument("--n-ensemble", type=int, default=3,
                       help="Ensemble models during search (default: 3, for speed)")
    parser.add_argument("--n-final-folds", type=int, default=5,
                       help="CV folds for final evaluation (default: 5)")
    parser.add_argument("--n-final-ensemble", type=int, default=5,
                       help="Ensemble models for final evaluation (default: 5)")
    parser.add_argument("--study-name", default=None,
                       help="Optuna study name (for resume/continue)")
    parser.add_argument("--storage", default=None,
                       help="Optuna storage URL (e.g., sqlite:///optuna.db)")
    parser.add_argument("--pruning", action="store_true", default=True,
                       help="Enable MedianPruner (default: on)")
    parser.add_argument("--no-pruning", action="store_false", dest="pruning",
                       help="Disable MedianPruner")
    parser.add_argument("--timeout", type=int, default=None,
                       help="Max search time in seconds")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    log("=" * 60)
    log("GC-NKGraph-Atlas BAYESIAN HYPERPARAMETER SEARCH")
    log(f"  Trials: {args.n_trials}")
    log(f"  Search folds: {args.n_folds}")
    log(f"  Search ensemble: {args.n_ensemble}")
    log(f"  Final folds: {args.n_final_folds}")
    log(f"  Final ensemble: {args.n_final_ensemble}")
    log("=" * 60)

    out_dir = ensure_dir(args.output_dir)

    # --- Validate graph exists ---
    graph_dir = args.graph_dir
    if not os.path.exists(os.path.join(graph_dir, "nodes.tsv")):
        log("Graph not found — cannot optimize GNN hyperparameters")
        log("Run the 'graph' pipeline phase first.")
        sys.exit(1)

    # --- Load data ---
    log("\n=== Loading data ===")
    expr, labels = _load_training_data()
    y_full = (labels["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values
    n_pos = y_full.sum()
    n_neg = len(y_full) - n_pos
    log(f"  Samples: {len(y_full)} (positive={n_pos}, negative={n_neg})")

    # --- Gene selection ---
    graph = _build_nx_graph(graph_dir)
    graph_genes = set(graph["node_to_idx"].keys())
    selected_genes = select_informative_genes(expr, graph_genes)
    X_expr = expr[selected_genes].values.astype(np.float32)
    log(f"  Features: {X_expr.shape[1]} genes")

    # --- Setup search context ---
    ctx = SearchContext(
        X_expr=X_expr,
        y_full=y_full,
        graph_dir=graph_dir,
        selected_genes=selected_genes,
        n_folds=args.n_folds,
        n_ensemble=args.n_ensemble,
    )

    # --- Create Optuna study ---
    study_name = args.study_name or f"gc_nkgraph_bayesian_{int(time.time())}"
    storage = args.storage

    if storage:
        log(f"\nUsing storage: {storage}")
        study = optuna.create_study(
            study_name=study_name,
            storage=storage,
            direction="maximize",
            sampler=optuna.samplers.TPESampler(seed=args.seed),
            pruner=optuna.pruners.MedianPruner(
                n_startup_trials=10,
                n_warmup_steps=1,
                interval_steps=1,
            ) if args.pruning else optuna.pruners.NopPruner(),
            load_if_exists=True,
        )
    else:
        study = optuna.create_study(
            study_name=study_name,
            direction="maximize",
            sampler=optuna.samplers.TPESampler(seed=args.seed),
            pruner=optuna.pruners.MedianPruner(
                n_startup_trials=10,
                n_warmup_steps=1,
                interval_steps=1,
            ) if args.pruning else optuna.pruners.NopPruner(),
        )

    log(f"\nStudy: {study_name}")
    log(f"Direction: maximize MCC")
    log(f"Sampler: TPE (seed={args.seed})")
    log(f"Pruner: {'MedianPruner' if args.pruning else 'None'}")

    # --- Run optimization ---
    log(f"\n=== Running {args.n_trials} trials ===\n")

    def objective(trial):
        return trial_objective(trial, ctx)

    start_time = time.time()
    study.optimize(objective, n_trials=args.n_trials, timeout=args.timeout)
    elapsed = time.time() - start_time

    # --- Print search summary ---
    log(f"\n{'='*60}")
    log("SEARCH COMPLETE")
    log(f"{'='*60}")
    log(f"  Elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    log(f"  Completed trials: {len(study.trials)}")
    log(f"  Pruned trials: {sum(1 for t in study.trials if t.state == TrialState.PRUNED)}")
    log(f"  Failed trials: {sum(1 for t in study.trials if t.state == TrialState.FAIL)}")

    # Best trial
    log(f"\nBest trial (#{study.best_trial.number}):")
    log(f"  MCC: {study.best_value:.4f}")
    for k, v in study.best_params.items():
        if k == "hidden_architecture":
            log(f"  {k}: {HIDDEN_ARCHITECTURES[v]}")
        elif k == "edge_weight_config":
            log(f"  {k}: {dict(zip(EDGE_TYPES_ORDER, EDGE_WEIGHT_CONFIGS[v]))}")
        else:
            log(f"  {k}: {v}")

    # --- Save trial history ---
    trials_df = study.trials_dataframe()
    trials_path = os.path.join(out_dir, "gc_nkgraph_bayesian_trials.tsv")
    trials_df.to_csv(trials_path, sep="\t", index=False)
    log(f"\nTrial history saved: {trials_path}")

    # --- Top-10 trials ---
    log(f"\n=== Top 10 Trials ===")
    completed = [t for t in study.trials if t.state == TrialState.COMPLETE]
    top = sorted(completed, key=lambda t: t.value, reverse=True)[:10]
    for rank, t in enumerate(top, 1):
        arch = HIDDEN_ARCHITECTURES[t.params["hidden_architecture"]]
        ew = EDGE_WEIGHT_CONFIGS[t.params["edge_weight_config"]]
        log(f"  #{rank}: MCC={t.value:.4f}  "
            f"emb={t.params['embedding_dim']}  "
            f"arch={arch}  "
            f"drop={t.params['dropout']:.1f}  "
            f"lr={t.params['learning_rate']:.2e}  "
            f"wd={t.params['weight_decay']:.2e}  "
            f"bs={t.params['batch_size']}")

    # --- Final evaluation with best params ---
    log(f"\n=== Final Evaluation (5-fold, 5-ensemble) ===")
    final_results = evaluate_best_params(
        ctx=ctx,
        best_params=study.best_params,
        output_dir=out_dir,
        n_folds=args.n_final_folds,
        n_ensemble=args.n_final_ensemble,
    )

    log(f"\n{'='*60}")
    log("BAYESIAN HYPERPARAMETER OPTIMIZATION COMPLETE!")
    log(f"{'='*60}")


if __name__ == "__main__":
    main()
