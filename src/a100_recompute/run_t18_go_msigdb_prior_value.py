"""
T18 -- External-sample value of GO-prior / MSigDB-prior edges (exploratory).

Context: Methods 2.5 notes that the ChEA TF-target edge type contributes 0
edges within the ~100-gene axis panel (too sparse to fire). This is an
exploratory follow-up, NOT part of the reported manuscript graph, testing
whether two alternative generic-prior sources -- inspired by
github.com/nblvguohao/CANOPY-Router / CANOPYNet-KBS (same author, unpublished,
different biology) -- add measurable cross-cohort predictive value the same
way T17 tested metabolic_crosstalk:
  - go_prior: GO Biological Process 2023 term co-membership (Jaccard >= 0.05),
    fixed weight 0.2 (mirrors TreeNet's GO-prior augmentation).
  - msigdb_prior: MSigDB C2 curated-pathway co-membership (Jaccard >= 0.05),
    fixed weight 0.2, matrix reused from the CANOPY-Router bundle.

Test design (same as T17, generalized to 4 graph variants):
  1. Graph variants already built by build_heterograph.py into
     data/processed/graph_ablation/{baseline,go,msigdb,both}/.
  2. For each variant: SVD gene embeddings -> train NK-state classifier on
     TCGA-STAD (gastric) -> evaluate cross-cohort on TCGA-LIHC (liver,
     positive control), same direction as T17.
  3. Bootstrap the MCC difference (variant - baseline) over 1000 resamples of
     the LIHC test set; report empirical p-value and 95% CI.

Criterion: variant > baseline in cross-cohort MCC with bootstrap p < 0.05 ->
the added prior edges provide externally measurable predictive value on top
of the manuscript's real-prior graph.

Output: results/tables/t18_go_msigdb_prior_value.tsv

Run:  python src/graph_construction/build_heterograph.py --out-dir data/processed/graph_ablation/baseline
      python src/graph_construction/build_heterograph.py --out-dir data/processed/graph_ablation/go --enable-go-prior
      python src/graph_construction/build_heterograph.py --out-dir data/processed/graph_ablation/msigdb --enable-msigdb-prior
      python src/graph_construction/build_heterograph.py --out-dir data/processed/graph_ablation/both --enable-go-prior --enable-msigdb-prior
      python src/a100_recompute/run_t18_go_msigdb_prior_value.py
"""
import sys, time
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir, load_table
from src.common.seed import set_seed
from src.a100_recompute.run_t17_edge_external_value import (
    build_adj_matrix, compute_svd_embeddings, load_cohort_expression,
    prepare_features, prepare_features_raw, standardize_with_params,
    train_and_evaluate_cross_cohort, bootstrap_mcc_difference,
)

GRAPH_ROOT = "data/processed/graph_ablation"
VARIANTS = ["baseline", "go", "msigdb", "both"]
SEEDS = [1234, 2345, 3456]
RESULTS = "results/tables"
OUT_PATH = f"{RESULTS}/t18_go_msigdb_prior_value.tsv"
OUT_PATH_SEEDS = f"{RESULTS}/t18_go_msigdb_prior_value_by_seed.tsv"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    log("=" * 60)
    log("T18 -- GO-prior / MSigDB-prior EDGE EXTERNAL-SAMPLE VALUE TEST (exploratory)")
    log("=" * 60)

    log("\n[1] Loading NK state labels...")
    labels = pd.read_csv("results/tables/nk_state_labels.tsv", sep="\t", index_col=0, comment="#")
    labels = labels[labels.index.notnull()]

    log("\n[2] Preparing train (STAD -> gastric) and test (LIHC -> liver) data across 3 seeds...")
    log("    (same direction as T17: gastric-trained, liver-evaluated)")
    log("    Train MCC=1.0 is expected (small STAD training set, MLP overfits fully);")
    log("    what matters is TEST (held-out LIHC) MCC variance across seeds.")

    expr_train, y_train_data = load_cohort_expression(labels, "TCGA-STAD")
    expr_test, y_test_data = load_cohort_expression(labels, "TCGA-LIHC")
    if expr_train is None or expr_test is None:
        log("ERROR: missing cohort expression data")
        return

    edge_counts = {}
    seed_results = {v: [] for v in VARIANTS}  # variant -> list of result dicts, one per seed
    per_seed_rows = []

    for variant in VARIANTS:
        vdir = f"{GRAPH_ROOT}/{variant}"
        nodes = pd.read_csv(f"{vdir}/nodes.tsv", sep="\t")
        edges = pd.read_csv(f"{vdir}/edges.tsv", sep="\t")
        adj, node_to_idx, idx_to_node = build_adj_matrix(edges, nodes)
        n_edges = int((adj > 0).sum() / 2)
        edge_counts[variant] = n_edges
        log(f"\n  --- {variant}: {len(nodes)} nodes, {n_edges} edges ---")

        # SVD embedding is deterministic given the (seed-independent) adjacency.
        emb = compute_svd_embeddings(adj, embedding_dim=64)

        for seed in SEEDS:
            set_seed(seed)
            X_train_full, mean_tr, std_tr = prepare_features(expr_train, emb, node_to_idx)
            X_test_full = prepare_features_raw(expr_test, emb, node_to_idx)
            X_test_full_s = standardize_with_params(X_test_full, mean_tr, std_tr)

            result = train_and_evaluate_cross_cohort(
                X_train_full, y_train_data, X_test_full_s, y_test_data,
                hidden_dims=(128, 64), epochs=200, batch_size=32, lr=1e-3,
            )
            log(f"    seed={seed}  Train MCC={result['train_mcc']:.4f}  "
                f"Test MCC={result['test_mcc']:.4f}  AUROC={result['test_auroc']:.4f}")
            seed_results[variant].append(result)
            per_seed_rows.append({
                "variant": variant, "seed": seed, "n_edges": n_edges,
                "train_mcc": round(result["train_mcc"], 4),
                "test_mcc": round(result["test_mcc"], 4),
                "test_auroc": round(result["test_auroc"], 4),
            })

    ensure_dir(RESULTS)
    pd.DataFrame(per_seed_rows).to_csv(OUT_PATH_SEEDS, sep="\t", index=False)
    log(f"\n[3] Per-seed results written to {OUT_PATH_SEEDS}")

    log(f"\n[4] Aggregate (mean +/- std over {len(SEEDS)} seeds) and paired test vs baseline...")

    out_rows = []
    baseline_mccs = np.array([r["test_mcc"] for r in seed_results["baseline"]])
    for v in VARIANTS:
        mccs = np.array([r["test_mcc"] for r in seed_results[v]])
        aurocs = np.array([r["test_auroc"] for r in seed_results[v]])
        row = {
            "variant": v, "n_edges": edge_counts[v], "n_seeds": len(SEEDS),
            "test_mcc_mean": round(mccs.mean(), 4), "test_mcc_std": round(mccs.std(), 4),
            "test_auroc_mean": round(aurocs.mean(), 4), "test_auroc_std": round(aurocs.std(), 4),
            "mean_diff_vs_baseline": "", "wilcoxon_p": "", "verdict": "reference" if v == "baseline" else "",
        }
        if v != "baseline":
            from scipy.stats import wilcoxon
            diff = mccs - baseline_mccs
            row["mean_diff_vs_baseline"] = round(diff.mean(), 4)
            try:
                _, p = wilcoxon(mccs, baseline_mccs)
            except ValueError:
                p = float("nan")  # all-zero differences or n<1 after ties
            row["wilcoxon_p"] = round(p, 4) if p == p else "NA"
            # n=3 seeds: treat as indicative, not definitive (same caveat TreeNet
            # itself makes for its own n=3 paired tests).
            row["verdict"] = (
                "indicative_positive" if diff.mean() > 0 and (p == p and p < 0.2)
                else ("indicative_negative" if diff.mean() < 0 and (p == p and p < 0.2)
                      else "inconclusive")
            )
        out_rows.append(row)
        log(f"  {v}: test_mcc = {row['test_mcc_mean']:.4f} +/- {row['test_mcc_std']:.4f}"
            + (f"   diff_vs_baseline={row['mean_diff_vs_baseline']:.4f}  "
               f"wilcoxon_p={row['wilcoxon_p']}  [{row['verdict']}]" if v != "baseline" else ""))

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(OUT_PATH, sep="\t", index=False)
    log(f"\n[5] Aggregate results written to {OUT_PATH}")

    log("\n" + "=" * 60)
    log("T18 COMPLETE (exploratory -- not part of the reported manuscript graph)")
    log("With only n=3 seeds, treat p-values as indicative per TreeNet's own")
    log("convention for small-n paired tests, not as confirmatory evidence.")
    log("=" * 60)


if __name__ == "__main__":
    main()
