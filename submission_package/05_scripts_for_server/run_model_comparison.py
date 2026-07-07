"""
Model comparison harness: GC-NKGraph-Atlas GNN vs. tabular baselines (Blocker B).

PROBLEM
-------
Methods 2.7 + Table 2 promise a GNN-vs-baseline comparison, but
`baseline_internal_results.tsv` was never produced. The GNN's solid internal
numbers (`gc_nkgraph_internal_results.tsv`, acc~0.87/AUROC~0.94) stand alone and
therefore mean nothing without a baseline reference.

VALIDITY NOTE (verified against source)
---------------------------------------
`gc_nkgraph_atlas.py` and `run_all_baselines.py` both define the label as
`nk_immune_state == "NK-hot-cytotoxic"` and both use
`StratifiedKFold(n_splits=5, shuffle=True, random_state=42)` over the same
`_load_training_data()` ordering. Folds are therefore identical by construction,
so per-fold paired tests between the GNN and each baseline are valid — provided
you did not change the data ordering or the label between the two runs.

WHAT THIS DOES
--------------
1. Runs the 6 tabular baselines on the seed-42 folds (reusing the existing,
   tested `run_all_baselines` functions).
2. Loads the existing GNN per-fold results and merges them in as one method.
3. Writes `results/tables/model_comparison.tsv` (per-fold) and
   `results/tables/model_comparison_summary.tsv` (mean+-std per method),
   plus paired Wilcoxon + paired t-tests of the GNN vs each baseline on MCC and
   AUROC across the 5 folds.

USAGE (run on server)
---------------------
    python src/baselines/run_model_comparison.py \
        --gnn-results results/tables/gc_nkgraph_internal_results.tsv

If the GNN results file is absent, run `src/models/gc_nkgraph_atlas.py` first.

NOTE: untested locally (no raw expression matrix + no xgboost/lightgbm here).
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import load_config, ensure_dir
from src.common.seed import set_seed
from src.baselines.run_all_baselines import (
    load_training_data, run_xgboost, run_lightgbm, run_random_forest,
    run_elasticnet, run_svm, run_mlp,
)

BASELINES = {
    "XGBoost": run_xgboost,
    "LightGBM": run_lightgbm,
    "RandomForest": run_random_forest,
    "ElasticNet": run_elasticnet,
    "SVM": run_svm,
    "MLP": run_mlp,
}
GNN_NAME = "GC-NKGraph-Atlas"


def run_baselines(logger):
    from sklearn.model_selection import StratifiedKFold
    set_seed(42)
    X_full, y_df = load_training_data(load_config("configs/data_config.yaml"), logger)
    y = (y_df["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    rows = []
    for fold, (tr, te) in enumerate(skf.split(X_full, y)):
        X_tr, X_te = X_full.iloc[tr], X_full.iloc[te]
        y_tr, y_te = y[tr], y[te]
        for name, fn in BASELINES.items():
            try:
                m = fn(X_tr, y_tr, X_te, y_te)
                rows.append({"method": name, "fold": fold, **m})
                logger.ok("COMPARISON",
                          f"{name} fold {fold}: MCC={m.get('MCC',float('nan')):.3f} "
                          f"AUROC={m.get('AUROC',float('nan')):.3f}", script=__file__)
            except Exception as e:
                logger.fail("COMPARISON", f"{name} fold {fold} FAILED: {e}", script=__file__)
    return pd.DataFrame(rows)


def load_gnn(path, logger):
    if not os.path.exists(path):
        logger.fail("COMPARISON",
                    f"GNN results not found at {path}; run gc_nkgraph_atlas.py first",
                    script=__file__)
        return pd.DataFrame()
    g = pd.read_csv(path, sep="\t")
    # normalize column names to match baseline metric keys
    rename = {"fold": "fold", "Accuracy": "Accuracy", "BalancedAccuracy": "BalancedAccuracy",
              "MacroF1": "MacroF1", "MCC": "MCC", "AUROC": "AUROC", "AUPRC": "AUPRC"}
    g = g.rename(columns=rename)
    g["method"] = GNN_NAME
    logger.ok("COMPARISON", f"loaded GNN results: {len(g)} folds from {path}", script=__file__)
    return g


def paired_tests(df, metric, logger):
    """Paired GNN-vs-baseline tests across folds for one metric."""
    from scipy import stats
    out = []
    gnn = df[df["method"] == GNN_NAME].sort_values("fold")[metric].values
    if len(gnn) == 0:
        return pd.DataFrame()
    for name in BASELINES:
        base = df[df["method"] == name].sort_values("fold")[metric].values
        if len(base) != len(gnn) or len(base) < 2:
            continue
        try:
            w_stat, w_p = stats.wilcoxon(gnn, base)
        except Exception:
            w_stat, w_p = np.nan, np.nan
        t_stat, t_p = stats.ttest_rel(gnn, base)
        out.append({
            "metric": metric, "baseline": name,
            "gnn_mean": float(np.mean(gnn)), "baseline_mean": float(np.mean(base)),
            "delta": float(np.mean(gnn) - np.mean(base)),
            "wilcoxon_p": w_p, "ttest_p": t_p,
        })
    res = pd.DataFrame(out)
    if not res.empty:
        logger.ok("COMPARISON",
                  f"{metric}: GNN vs baselines paired tests computed "
                  f"({len(res)} comparisons)", script=__file__)
    return res


def main():
    ap = argparse.ArgumentParser(description="GNN vs baselines comparison (Blocker B)")
    ap.add_argument("--gnn-results", default="results/tables/gc_nkgraph_internal_results.tsv")
    ap.add_argument("--tables-dir", default="results/tables")
    ap.add_argument("--skip-baselines", action="store_true",
                    help="reuse existing baseline_internal_results.tsv instead of rerunning")
    args = ap.parse_args()

    logger = Logger()
    tables_dir = ensure_dir(args.tables_dir)

    if args.skip_baselines and os.path.exists(os.path.join(tables_dir, "baseline_internal_results.tsv")):
        base_df = pd.read_csv(os.path.join(tables_dir, "baseline_internal_results.tsv"), sep="\t")
    else:
        base_df = run_baselines(logger)
        base_df.to_csv(os.path.join(tables_dir, "baseline_internal_results.tsv"),
                       sep="\t", index=False)

    gnn_df = load_gnn(args.gnn_results, logger)
    metric_cols = ["Accuracy", "BalancedAccuracy", "MacroF1", "MCC", "AUROC", "AUPRC"]
    keep = ["method", "fold"] + [c for c in metric_cols if c in base_df.columns]
    combined = pd.concat([base_df[keep], gnn_df[[c for c in keep if c in gnn_df.columns]]],
                         ignore_index=True)
    combined.to_csv(os.path.join(tables_dir, "model_comparison.tsv"), sep="\t", index=False)

    # per-method summary (mean +/- std)
    summary = (combined.groupby("method")[[c for c in metric_cols if c in combined.columns]]
               .agg(["mean", "std"]).round(4))
    summary.to_csv(os.path.join(tables_dir, "model_comparison_summary.tsv"), sep="\t")

    # paired significance on the headline metrics
    stats_frames = [paired_tests(combined, m, logger) for m in ("MCC", "AUROC")]
    stats_frames = [s for s in stats_frames if not s.empty]
    if stats_frames:
        pd.concat(stats_frames, ignore_index=True).to_csv(
            os.path.join(tables_dir, "model_comparison_stats.tsv"), sep="\t", index=False)

    print("\n=== MODEL COMPARISON (mean) ===")
    print(summary.xs("mean", axis=1, level=1).to_string())
    logger.ok("COMPARISON", f"comparison complete -> {tables_dir}/model_comparison*.tsv",
              script=__file__)


if __name__ == "__main__":
    main()
