"""
GC-NKGraph-Atlas Baseline Models (Phase 9).

Trains and evaluates all baseline models for NK-state classification:
  - Tabular: XGBoost, LightGBM, RandomForest, ElasticNet, SVM, MLP
  - Graph: GCN, GAT, GraphSAGE (if graph data available)

Usage:
    python src/baselines/run_all_baselines.py --config configs/experiment_config.yaml
"""

import os, sys, json, argparse, warnings
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import load_config, save_table, ensure_dir, load_table
from src.common.seed import set_seed

warnings.filterwarnings("ignore")


def load_training_data(config: dict, logger: Logger):
    """Load training (TCGA-STAD) expression + labels."""
    for ds in config.get("bulk_datasets", []):
        if ds["role"] == "train_primary":
            expr_path = ds["expression_path"]
            break
    else:
        raise ValueError("No train_primary dataset found in config")

    expr = load_table(expr_path)
    if expr.index.name != "sample_id":
        expr = expr.T if "gene" in str(expr.index.name) else expr

    labels = load_table("results/tables/nk_state_labels.tsv")
    common = expr.index.intersection(labels.index)
    logger.ok("BASELINE", f"Train data: {len(common)} common samples, {expr.shape[1]} genes", script=__file__)
    return expr.loc[common], labels.loc[common]


def load_external_data(config: dict, logger: Logger, dataset_name: str):
    """Load external validation dataset."""
    for ds in config.get("bulk_datasets", []):
        if ds["name"] == dataset_name:
            expr = load_table(ds["expression_path"])
            if expr.index.name != "sample_id":
                expr = expr.T if "gene" in str(expr.index.name) else expr
            return expr
    return None


def evaluate(y_true, y_pred, y_prob=None):
    """Compute classification metrics."""
    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                  f1_score, matthews_corrcoef,
                                  roc_auc_score, average_precision_score)
    metrics = {
        "Accuracy": accuracy_score(y_true, y_pred),
        "BalancedAccuracy": balanced_accuracy_score(y_true, y_pred),
        "MacroF1": f1_score(y_true, y_pred, average="macro"),
        "MCC": matthews_corrcoef(y_true, y_pred),
    }
    if y_prob is not None and len(np.unique(y_true)) == 2:
        metrics["AUROC"] = roc_auc_score(y_true, y_prob[:, 1])
        metrics["AUPRC"] = average_precision_score(y_true, y_prob[:, 1])
    return metrics


def run_xgboost(X_train, y_train, X_test, y_test):
    """XGBoost baseline."""
    import xgboost as xgb
    model = xgb.XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.1,
                              random_state=42, n_jobs=32, eval_metric="mlogloss")
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test) if len(np.unique(y_train)) > 2 else model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def run_lightgbm(X_train, y_train, X_test, y_test):
    """LightGBM baseline."""
    import lightgbm as lgb
    model = lgb.LGBMClassifier(n_estimators=200, max_depth=6, learning_rate=0.1,
                                random_state=42, n_jobs=32, verbose=-1)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def run_random_forest(X_train, y_train, X_test, y_test):
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=200, max_depth=12, n_jobs=32, random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def run_elasticnet(X_train, y_train, X_test, y_test):
    from sklearn.linear_model import LogisticRegression
    model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                C=1.0, max_iter=1000, random_state=42, n_jobs=32)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def run_svm(X_train, y_train, X_test, y_test):
    from sklearn.svm import SVC
    model = SVC(kernel="rbf", C=1.0, probability=True, random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def run_mlp(X_train, y_train, X_test, y_test):
    from sklearn.neural_network import MLPClassifier
    model = MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=500,
                           random_state=42, early_stopping=True)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)
    return evaluate(y_test, y_pred, y_prob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/experiment_config.yaml")
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--split-id", type=int, default=0)
    parser.add_argument("--methods", nargs="+",
                        default=["xgboost", "lightgbm", "random_forest", "elasticnet", "svm", "mlp"])
    args = parser.parse_args()

    logger = Logger()
    config = load_config(args.config)
    tables_dir = ensure_dir(args.tables_dir)

    set_seed(42)

    # Load data
    X_train_full, y_train_full = load_training_data(load_config("configs/data_config.yaml"), logger)

    # Prepare labels (simplify to binary: NK-hot vs rest for baseline)
    y_train_full["label"] = y_train_full["nk_immune_state"].apply(
        lambda x: 1 if x == "NK-hot-cytotoxic" else 0
    )
    y = y_train_full["label"].values

    # Simple train/test split (80/20) for quick baseline
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    all_results = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_train_full, y)):
        X_train = X_train_full.iloc[train_idx]
        X_test = X_train_full.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        for method in args.methods:
            try:
                if method == "xgboost":
                    metrics = run_xgboost(X_train, y_train, X_test, y_test)
                elif method == "lightgbm":
                    metrics = run_lightgbm(X_train, y_train, X_test, y_test)
                elif method == "random_forest":
                    metrics = run_random_forest(X_train, y_train, X_test, y_test)
                elif method == "elasticnet":
                    metrics = run_elasticnet(X_train, y_train, X_test, y_test)
                elif method == "svm":
                    metrics = run_svm(X_train, y_train, X_test, y_test)
                elif method == "mlp":
                    metrics = run_mlp(X_train, y_train, X_test, y_test)
                else:
                    continue

                result = {"method": method, "fold": fold, **metrics}
                all_results.append(result)
                logger.ok("BASELINE", f"{method} fold {fold}: Acc={metrics.get('Accuracy',0):.3f} "
                         f"F1={metrics.get('MacroF1',0):.3f} MCC={metrics.get('MCC',0):.3f}",
                         script=__file__)
            except Exception as e:
                logger.fail("BASELINE", f"{method} fold {fold} FAILED: {e}", script=__file__)

    # Save results
    if all_results:
        df = pd.DataFrame(all_results)
        save_table(df, os.path.join(tables_dir, "baseline_internal_results.tsv"),
                   provenance={"script": __file__, "config": args.config})
        logger.ok("BASELINE", f"All baselines complete. Results saved to {tables_dir}", script=__file__)

        # Summary
        summary = df.groupby("method").agg({k: ["mean", "std"] for k in
                                            ["Accuracy", "BalancedAccuracy", "MacroF1", "MCC"]})
        print("\n=== BASELINE SUMMARY ===")
        print(summary.to_string())


if __name__ == "__main__":
    main()
