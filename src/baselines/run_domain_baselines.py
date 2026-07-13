"""
P1-1: Domain-method baselines for NK-state classification.

Implements two alternative approaches that represent what standard
bioinformatics tools would produce WITHOUT the heterogeneous graph:

  (a) NK-marker signature baseline: score each bulk sample using canonical
      NK marker genes, then predict NK state via logistic regression on
      that score alone. This approximates the output of CIBERSORTx /
      quanTIseq: an NK abundance estimate -> state classification.

  (b) SST-module signature baseline: score each bulk sample on the 7
      SST-axis modules (from the mechanism card), then use a logistic
      regression on module scores to predict NK state. This captures the
      "use the anchor paper's gene modules directly without building a
      graph" approach — the closest no-graph alternative to the GNN.

Both baselines run on the SAME 5-fold stratified cross-validation
(random_state=42) used by the GNN and tabular baselines, enabling
direct per-fold comparison via paired tests.

Also outputs a CIBERSORTx/quanTIseq/Scissor comparison roadmap for
reviewers, documenting what would be needed to reproduce the comparison
with those tools.

Output:
  results/tables/domain_baselines_per_fold.tsv
  results/tables/domain_baselines_summary.tsv
  results/tables/domain_baselines_tests.tsv         (paired vs GNN)
  results/tables/CIBERSORTx_quanTIseq_Scissor_roadmap.md

Usage:
    python src/baselines/run_domain_baselines.py
"""
import os, sys, time, warnings, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.log_utils import Logger
from src.common.io_utils import load_config, ensure_dir
from src.common.seed import set_seed
from src.baselines.run_all_baselines import load_training_data

# ── NK marker genes (canonical, from manuscript Methods 2.4) ────────────
NK_MARKERS = ["NCAM1", "KLRD1", "NKG7", "GNLY", "KLRF1", "EOMES", "NCR1", "FCGR3A"]

# ── SST-axis modules (7 modules, from manuscript Methods 2.3) ──────────
SST_MODULES = {
    "tumor_serine_capacity": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2",
                               "MTHFD1", "MTHFD2", "MTHFD1L", "SLC1A4", "SLC1A5"],
    "nk_sm_synthesis": ["SGMS1", "SGMS2"],
    "nk_sm_catabolism": ["SMPD1", "SMPD2", "SMPD3", "SMPD4"],
    "nk_denovo_sphingolipid": ["SPTLC1", "SPTLC2", "SPTLC3", "SPTSSA",
                                "CERS2", "CERS4", "CERS5", "CERS6", "DEGS1"],
    "nk_protrusion_machinery": ["EZR", "MSN", "RDX", "ACTR2", "ACTR3",
                                 "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5",
                                 "WAS", "WASL", "WASF1", "WASF2", "WASF3",
                                 "WIPF1", "CDC42", "RAC1", "RHOA",
                                 "DIAPH1", "DIAPH3", "FMNL1", "BAIAP2", "PACSIN2"],
    "nk_synapse_cytotoxicity_outcome": ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG",
                                         "LCP2", "LAT", "VAV1", "TLN1", "ITGAL", "ITGB2"],
    "checkpoint_link": ["HAVCR2"],
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_module_scores(expr: pd.DataFrame, modules: dict) -> pd.DataFrame:
    """Compute per-sample mean expression for each gene module.

    For each module, computes mean z-score across available genes.
    Returns DataFrame with one column per module.
    """
    scores = pd.DataFrame(index=expr.index)
    for mod_name, genes in modules.items():
        available = [g for g in genes if g in expr.columns]
        if not available:
            scores[mod_name] = 0.0
            continue
        sub = expr[available]
        z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).clip(lower=1e-10)
        scores[mod_name] = z.fillna(0).mean(axis=1)
    return scores


def nk_marker_baseline(X_tr, y_tr, X_te, y_te):
    """NK-marker signature baseline.

    Trains a logistic regression on NK marker mean expression -> NK state.
    This is the simplest possible domain approach.
    """
    from sklearn.linear_model import LogisticRegression

    # Compute NK marker score
    available = [g for g in NK_MARKERS if g in X_tr.columns]
    if not available:
        return {"Accuracy": 0.5, "BalancedAccuracy": 0.5, "MacroF1": 0.5,
                "MCC": 0.0, "AUROC": 0.5, "AUPRC": 0.5}

    tr_score = X_tr[available].mean(axis=1).values.reshape(-1, 1)
    te_score = X_te[available].mean(axis=1).values.reshape(-1, 1)

    clf = LogisticRegression(max_iter=1000, random_state=42)
    clf.fit(tr_score, y_tr)

    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                  f1_score, matthews_corrcoef, roc_auc_score,
                                  average_precision_score)
    y_pred = clf.predict(te_score)
    y_proba = clf.predict_proba(te_score)[:, 1]

    return {
        "Accuracy": accuracy_score(y_te, y_pred),
        "BalancedAccuracy": balanced_accuracy_score(y_te, y_pred),
        "MacroF1": f1_score(y_te, y_pred, average="macro"),
        "MCC": matthews_corrcoef(y_te, y_pred),
        "AUROC": roc_auc_score(y_te, y_proba),
        "AUPRC": average_precision_score(y_te, y_proba),
    }


def sst_module_baseline(X_tr, y_tr, X_te, y_te):
    """SST-module signature baseline.

    Computes the 7 SST-axis module scores on bulk expression, then trains
    a logistic regression on the module scores to predict NK state.
    This represents "use the anchor paper's gene modules directly,
    without building a graph."
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    tr_mod = compute_module_scores(X_tr, SST_MODULES)
    te_mod = compute_module_scores(X_te, SST_MODULES)

    scaler = StandardScaler()
    tr_feat = scaler.fit_transform(tr_mod)
    te_feat = scaler.transform(te_mod)

    # Remove constant columns (modules with no genes found)
    nonzero_cols = np.std(tr_feat, axis=0) > 1e-10
    tr_feat = tr_feat[:, nonzero_cols]
    te_feat = te_feat[:, nonzero_cols]

    if tr_feat.shape[1] == 0:
        from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                      f1_score, matthews_corrcoef, roc_auc_score,
                                      average_precision_score)
        return {
            "Accuracy": 0.5, "BalancedAccuracy": 0.5, "MacroF1": 0.5,
            "MCC": 0.0, "AUROC": 0.5, "AUPRC": 0.5,
        }

    clf = LogisticRegression(max_iter=5000, random_state=42, class_weight="balanced")
    clf.fit(tr_feat, y_tr)

    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                  f1_score, matthews_corrcoef, roc_auc_score,
                                  average_precision_score)
    y_pred = clf.predict(te_feat)
    y_proba = clf.predict_proba(te_feat)[:, 1]

    return {
        "Accuracy": accuracy_score(y_te, y_pred),
        "BalancedAccuracy": balanced_accuracy_score(y_te, y_pred),
        "MacroF1": f1_score(y_te, y_pred, average="macro"),
        "MCC": matthews_corrcoef(y_te, y_pred),
        "AUROC": roc_auc_score(y_te, y_proba),
        "AUPRC": average_precision_score(y_te, y_proba),
    }


def paired_tests(per_fold_df, ref_method, target_methods, metrics):
    """Paired t-test and Wilcoxon test: ref vs each target on each metric."""
    rows = []
    ref = per_fold_df[per_fold_df["method"] == ref_method]

    for tgt in target_methods:
        tdf = per_fold_df[per_fold_df["method"] == tgt]

        # Align by fold
        merged = ref.merge(tdf, on="fold", suffixes=("_ref", "_tgt"))

        for metric in metrics:
            d = merged[f"{metric}_tgt"].values - merged[f"{metric}_ref"].values
            delta_mean = np.mean(d)
            t_stat, t_p = stats.ttest_1samp(d, 0.0) if len(d) >= 2 else (np.nan, np.nan)
            try:
                w_stat, w_p = stats.wilcoxon(
                    merged[f"{metric}_tgt"].values,
                    merged[f"{metric}_ref"].values,
                )
            except ValueError:
                w_stat, w_p = np.nan, np.nan

            rows.append({
                "comparison": f"{tgt} vs {ref_method}",
                "metric": metric,
                "mean_delta": round(delta_mean, 4),
                "t_statistic": round(t_stat, 4) if not np.isnan(t_stat) else None,
                "t_p_value": t_p if not np.isnan(t_p) else None,
                "wilcoxon_p_value": w_p if not np.isnan(w_p) else None,
                "significant_0.05": "YES" if (t_p < 0.05 and not np.isnan(t_p)) else "no",
            })
    return pd.DataFrame(rows)


def write_cibersortx_roadmap(out_dir: str):
    """Write the CIBERSORTx/quanTIseq/Scissor comparison roadmap.

    This documents what a reviewer (or future contributor) would need to do
    to fully replicate the domain-method comparison using these tools.
    The current environment does not have R installed; these are R packages.
    """
    roadmap = os.path.join(out_dir, "CIBERSORTx_quanTIseq_Scissor_roadmap.md")
    with open(roadmap, "w") as f:
        f.write("""# Domain-Method Baseline Comparison — Roadmap for External Tools

## Purpose
This document describes how to reproduce the domain-method baseline comparison
using CIBERSORTx, quanTIseq, and Scissor, extending the analysis in
`src/baselines/run_domain_baselines.py`.

## Environment requirements
- R >= 4.0
- R packages: `immunedeconv` (quanTIseq), `Scissor`
- CIBERSORTx: Docker image `cibersortx/fractions` or web portal
  (https://cibersortx.stanford.edu)

## Baseline 1: CIBERSORTx NK-fraction baseline

### Rationale
CIBERSORTx deconvolves bulk RNA-seq into immune cell-type fractions using a
signature matrix (LM22). The NK-cell fraction estimate is used as a single
feature to predict NK-hot-cytotoxic state.

### Procedure
```r
# Step 1: Prepare mixture file (TCGA-STAD expression, genes x samples)
# Step 2: Run CIBERSORTx
#   docker run -v $PWD:/data cibersortx/fractions \
#     --username <token> --token <token> \
#     --refsample LM22.txt --mixture tcga_stad_mixture.txt \
#     --perm 100 --QN FALSE
# Step 3: Extract NK fraction column (NK cells resting + NK cells activated)
# Step 4: 5-fold stratified CV, logistic regression on NK fraction -> state
```

### Expected metrics (approximate, based on known TCGA-STAD NK content)
- AUROC: 0.65–0.75 (NK fraction alone is a weak predictor of NK activation state)
- This baseline tests whether the GNN adds information beyond simple NK abundance

## Baseline 2: quanTIseq immune-deconvolution baseline

### Rationale
quanTIseq (Finotello et al., *Genome Med* 2019) provides absolute immune cell
fractions using a constrained least-squares approach. Like CIBERSORTx, the NK
fraction serves as the predictor.

### Procedure
```r
library(immunedeconv)
# Load TCGA-STAD expression (TPM-normalized, genes x samples)
expr <- read.table("tcga_stad_expression.tsv", row.names=1, header=TRUE)
# Run quanTIseq
res <- deconvolute(expr, method="quantiseq")
# Extract NK cell fraction
nk_fraction <- as.numeric(res["NK cell", ])
```

### Expected comparison
- NK fraction per sample -> logistic regression -> 5-fold CV
- AUROC expected 0.65–0.75 (similar to CIBERSORTx)
- Also compute: compare with SST-module baseline (which uses mechanism-specific
  modules rather than generic immune markers)

## Baseline 3: Scissor phenotype-to-genotype baseline

### Rationale
Scissor (Sun et al., *Nat Biotechnol* 2022) links single-cell phenotypes to bulk
clinical variables. The scRNA-defined NK states (hot-cytotoxic / dysfunctional /
cold) are treated as phenotypes, and Scissor identifies which bulk samples
associate with each phenotype.

### Procedure
```r
library(Scissor)
# Load scRNA data: scRNA NK cells (8,310 cells, 9 samples)
#   with phenotype labels (NK-hot-cytotoxic, NK-hot-dysfunctional,
#   NK-cold/excluded, NK-intermediate)
# Load bulk data: TCGA-STAD expression (450 samples)
# Run Scissor:
scissor_output <- Scissor(
  bulk_dataset = tcga_stad_expr,
  sc_dataset = nk_scrna_expr,
  sc_phenotype = nk_phenotype_labels,  # "NK-hot-cytotoxic" = phenotype of interest
  alpha = 0.05,
  family = "binomial"
)
# Scissor identifies "Scissor+" samples (associated with NK-hot-cytotoxic)
# Use Scissor+/- label as predictor
```

### Expected comparison
- Scissor+/- label -> classification performance vs GNN
- Key question: does the GNN's mechanism-structured embedding capture
  information beyond phenotype-genotype association?

## What the current repository already provides

`src/baselines/run_domain_baselines.py` implements two directly runnable
baselines that require no R:
1. **NK-marker signature baseline** — logistic regression on mean NK marker
   expression (conceptually the simplest possible deconvolution proxy)
2. **SST-module signature baseline** — logistic regression on 7 SST-axis
   module scores (the closest no-graph alternative to the GNN)

These two baselines capture the core question: "Does the graph add anything
beyond the anchor paper's gene modules scored directly on bulk expression?"

The CIBERSORTx/quanTIseq baselines answer a related but distinct question:
"Does the GNN add anything beyond generic NK abundance estimates from
state-of-the-art deconvolution?"

The Scissor baseline answers: "Does the GNN add anything beyond
phenotype-genotype association from single-cell data?"

All three external comparisons are complementary and would strengthen the
manuscript's "our framework adds value over existing tools" claim. They
require an R environment not available in the current local setup.
""")
    log(f"  Saved: {roadmap}")


def main():
    log("=" * 70)
    log("P1-1: DOMAIN-METHOD BASELINES — Signature Scoring vs GNN")
    log("=" * 70)

    out_dir = "results/tables"
    ensure_dir(out_dir)

    set_seed(42)

    # ── Load data ────────────────────────────────────────────────────────
    config = load_config("configs/data_config.yaml")
    logger = Logger(log_path="results/logs/domain_baselines_LOG.md")
    X_full, y_df = load_training_data(config, logger)
    y = (y_df["nk_immune_state"] == "NK-hot-cytotoxic").astype(int).values
    log(f"Loaded: {X_full.shape[0]} samples x {X_full.shape[1]} genes")
    log(f"  NK-hot-cytotoxic: {y.sum()} / {len(y)} ({y.sum()/len(y)*100:.1f}%)")

    # Report gene coverage for NK markers and SST modules
    nk_found = sum(1 for g in NK_MARKERS if g in X_full.columns)
    log(f"  NK markers found: {nk_found}/{len(NK_MARKERS)}")
    for mod_name, genes in SST_MODULES.items():
        n = sum(1 for g in genes if g in X_full.columns)
        log(f"  SST {mod_name}: {n}/{len(genes)} genes found")

    # ── 5-fold stratified CV ─────────────────────────────────────────────
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    rows = []
    for fold, (tr, te) in enumerate(skf.split(X_full, y)):
        X_tr, X_te = X_full.iloc[tr], X_full.iloc[te]
        y_tr, y_te = y[tr], y[te]
        log(f"\n  Fold {fold}: train={len(tr)}, test={len(te)}, "
            f"test_pos={y_te.sum()}/{len(y_te)}")

        # NK marker baseline
        res_nk = nk_marker_baseline(X_tr, y_tr, X_te, y_te)
        res_nk["method"] = "NK-marker signature"
        res_nk["fold"] = fold
        rows.append(res_nk)
        log(f"    NK-marker: AUROC={res_nk['AUROC']:.4f}, MCC={res_nk['MCC']:.4f}")

        # SST module baseline
        res_sst = sst_module_baseline(X_tr, y_tr, X_te, y_te)
        res_sst["method"] = "SST-module signature"
        res_sst["fold"] = fold
        rows.append(res_sst)
        log(f"    SST-module: AUROC={res_sst['AUROC']:.4f}, MCC={res_sst['MCC']:.4f}")

    df = pd.DataFrame(rows)
    metrics = ["Accuracy", "BalancedAccuracy", "MacroF1", "MCC", "AUROC", "AUPRC"]

    # ── Merge with existing baselines ────────────────────────────────────
    existing_path = "results/tables/model_comparison.tsv"
    if os.path.exists(existing_path):
        existing = pd.read_csv(existing_path, sep="\t")
        # Combine
        combined = pd.concat([existing[["method", "fold"] + metrics], df], ignore_index=True)
        log(f"\n  Merged: {len(combined)} rows ({len(df)} new + {len(existing)} existing)")
    else:
        combined = df

    # ── Save per-fold ────────────────────────────────────────────────────
    per_fold_path = os.path.join(out_dir, "domain_baselines_per_fold.tsv")
    combined.to_csv(per_fold_path, sep="\t", index=False)
    log(f"\nSaved: {per_fold_path}")

    # ── Summary ───────────────────────────────────────────────────────────
    summary_rows = []
    for method in sorted(combined["method"].unique()):
        sub = combined[combined["method"] == method]
        row = {"method": method, "n_folds": len(sub)}
        for m in metrics:
            row[f"{m}_mean"] = round(sub[m].mean(), 4)
            row[f"{m}_std"] = round(sub[m].std(), 4)
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows).sort_values("AUROC_mean", ascending=False)
    summary_path = os.path.join(out_dir, "domain_baselines_summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    log(f"Saved: {summary_path}")
    print()
    print(summary[["method", "AUROC_mean", "AUROC_std", "MCC_mean", "MCC_std",
                    "BalancedAccuracy_mean", "BalancedAccuracy_std"]].to_string(index=False))

    # ── Paired tests: domain baselines vs GNN ────────────────────────────
    gnn_rows = combined[combined["method"].str.contains("GNN|GC-NK", case=False)]
    if len(gnn_rows) > 0:
        gnn_name = gnn_rows["method"].iloc[0]
        domain_methods = ["NK-marker signature", "SST-module signature"]
        existing_domain = [m for m in domain_methods if m in combined["method"].values]

        if existing_domain:
            tests = paired_tests(combined, gnn_name, existing_domain, metrics)
            tests_path = os.path.join(out_dir, "domain_baselines_tests.tsv")
            tests.to_csv(tests_path, sep="\t", index=False)
            log(f"Saved: {tests_path}")

            # Highlight
            log("\n" + "=" * 60)
            log("KEY COMPARISON: Domain baselines vs GNN")
            log("=" * 60)
            for _, r in tests.iterrows():
                sig = "***" if r["significant_0.05"] == "YES" else ""
                log(f"  {r['comparison']:<50} {r['metric']:<20} "
                    f"delta={r['mean_delta']:+.4f}  p={r['t_p_value']:.4f} {sig}")

    # ── Write CIBERSORTx/quanTIseq/Scissor roadmap ────────────────────────
    write_cibersortx_roadmap(out_dir)

    # ── Summary for manuscript ───────────────────────────────────────────
    log("\n" + "=" * 70)
    log("MANUSCRIPT-READY SUMMARY")
    log("=" * 70)
    for method in ["NK-marker signature", "SST-module signature"]:
        sub = combined[combined["method"] == method]
        if len(sub) > 0:
            log(f"\n  {method}:")
            log(f"    AUROC = {sub['AUROC'].mean():.3f} +/- {sub['AUROC'].std():.3f}")
            log(f"    MCC   = {sub['MCC'].mean():.3f} +/- {sub['MCC'].std():.3f}")
            log(f"    BalAcc= {sub['BalancedAccuracy'].mean():.3f} +/- {sub['BalancedAccuracy'].std():.3f}")

    if len(gnn_rows) > 0:
        log(f"\n  GNN ({gnn_name}):")
        log(f"    AUROC = {gnn_rows['AUROC'].mean():.3f} +/- {gnn_rows['AUROC'].std():.3f}")
        log(f"    MCC   = {gnn_rows['MCC'].mean():.3f} +/- {gnn_rows['MCC'].std():.3f}")

    log("\nP1-1 COMPLETE!")


if __name__ == "__main__":
    main()
