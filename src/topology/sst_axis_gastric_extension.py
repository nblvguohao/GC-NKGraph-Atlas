"""
GC-NKGraph-Atlas SST Axis Gastric Extension (Phase 14R.6).
Arm B — Test the SST axis in gastric cancer as the novel result.

Strategy:
  1. Compare NK SST scores: Gastric Cancer vs Healthy Liver
  2. Within GC NK cells: test H2-H4 axis coherence
  3. Identify tumor-intrinsic SST-axis candidate targets
  4. Report axis recovery status in GC

Usage:
    python src/topology/sst_axis_gastric_extension.py
"""

import os, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

MODULES = {
    "tumor_serine_capacity": ["PHGDH","PSAT1","PSPH","SHMT1","SHMT2",
                               "MTHFD1","MTHFD2","MTHFD1L","SLC1A4","SLC1A5"],
    "nk_sm_synthesis": ["SGMS1","SGMS2"],
    "nk_sm_catabolism": ["SMPD1","SMPD2","SMPD3","SMPD4"],
    "nk_protrusion_machinery": ["EZR","MSN","RDX",
                                 "ACTR2","ACTR3","ARPC1B","ARPC2","ARPC3","ARPC4","ARPC5",
                                 "WAS","WASL","WASF1","WASF2","WASF3","WIPF1",
                                 "CDC42","RAC1","RHOA","DIAPH1","DIAPH3","FMNL1",
                                 "BAIAP2","PACSIN2"],
    "nk_synapse_cytotoxicity_outcome": ["NKG7","GNLY","GZMB","PRF1","IFNG",
                                         "LCP2","LAT","VAV1","TLN1","ITGAL","ITGB2"],
    "checkpoint_link": ["HAVCR2"],
}

SST_CORE_GENES = (MODULES["nk_sm_synthesis"] + MODULES["nk_sm_catabolism"] +
                   MODULES["nk_protrusion_machinery"] + MODULES["nk_synapse_cytotoxicity_outcome"])

def mean_zscore(expr, genes):
    available = [g for g in genes if g in expr.columns]
    if not available:
        return pd.Series(0.0, index=expr.index)
    z = (expr[available] - expr[available].mean(0)) / expr[available].std(0, ddof=0)
    return z.fillna(0).mean(axis=1)

def main():
    log("=" * 60)
    log("SST AXIS — Arm B: GASTRIC CANCER EXTENSION")
    log("=" * 60)

    out_dir = "results/tables"
    fig_dir = "results/figures"
    integrated_path = "data/processed/scrna/gc_integrated.h5ad"
    os.makedirs(out_dir, exist_ok=True)

    import scanpy as sc
    import scipy.sparse
    adata = sc.read(integrated_path)
    log(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    if scipy.sparse.issparse(adata.X):
        expr = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    else:
        expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    is_nk = adata.obs["cell_type"] == "NK"
    is_gc = adata.obs["tissue"] == "gastric_cancer"
    is_hl = adata.obs["tissue"] == "healthy_liver"

    # ---- 1. NK SST scores: GC vs HL ----
    log("\n--- 1. SST scores: Gastric Cancer vs Healthy Liver NK ---")
    nk_expr = expr.loc[is_nk]
    nk_scores = pd.DataFrame(index=nk_expr.index)
    for name in ["nk_sm_synthesis", "nk_sm_catabolism", "nk_protrusion_machinery",
                  "nk_synapse_cytotoxicity_outcome", "checkpoint_link"]:
        nk_scores[f"{name}_score"] = mean_zscore(nk_expr, MODULES[name])
    nk_scores["nk_sm_balance_score"] = nk_scores["nk_sm_synthesis_score"] - nk_scores["nk_sm_catabolism_score"]
    nk_scores["nk_topology_permissive_score"] = (nk_scores["nk_sm_balance_score"] + nk_scores["nk_protrusion_machinery_score"]) / 2
    nk_scores["tissue"] = adata.obs.loc[is_nk, "tissue"].values
    nk_scores["sample_id"] = adata.obs.loc[is_nk, "sample_id"].values

    gc_nk = nk_scores[nk_scores["tissue"] == "gastric_cancer"]
    hl_nk = nk_scores[nk_scores["tissue"] == "healthy_liver"]
    log(f"  GC NK: {len(gc_nk)} cells, HL NK: {len(hl_nk)} cells")

    metrics = ["nk_sm_balance_score", "nk_protrusion_machinery_score",
               "nk_topology_permissive_score", "nk_synapse_cytotoxicity_outcome_score"]
    comparison = []
    for m in metrics:
        gc_mean, hl_mean = gc_nk[m].mean(), hl_nk[m].mean()
        t_stat, p_val = stats.ttest_ind(gc_nk[m], hl_nk[m], equal_var=False)
        direction = "UP" if gc_mean > hl_mean else "DOWN"
        comparison.append({
            "metric": m, "GC_mean": round(gc_mean, 4), "HL_mean": round(hl_mean, 4),
            "delta": round(gc_mean - hl_mean, 4), "direction": direction,
            "p_value": f"{p_val:.4e}", "significant": p_val < 0.05
        })
        log(f"  {m}: GC={gc_mean:.4f} vs HL={hl_mean:.4f} ({direction}), p={p_val:.4e}")

    comp_df = pd.DataFrame(comparison)
    comp_df.to_csv(os.path.join(out_dir, "sst_axis_gc_vs_hl.tsv"), sep="\t", index=False)

    # ---- 2. Within-GC axis coherence (H2-H4) ----
    log("\n--- 2. SST axis coherence WITHIN gastric cancer NK ---")
    for label, nk_subset in [("GC NK", gc_nk), ("HL NK", hl_nk), ("ALL NK", nk_scores)]:
        log(f"\n  {label} (n={len(nk_subset)}):")
        r2, p2 = stats.pearsonr(nk_subset["nk_sm_balance_score"], nk_subset["nk_protrusion_machinery_score"])
        r3, p3 = stats.pearsonr(nk_subset["nk_protrusion_machinery_score"], nk_subset["nk_synapse_cytotoxicity_outcome_score"])
        r4, p4 = stats.pearsonr(nk_subset["nk_topology_permissive_score"], nk_subset["checkpoint_link_score"])
        log(f"    H2 SM->Protrusion: r={r2:.4f}, p={p2:.4e} {'PASS' if r2>0 and p2<0.05 else 'FAIL'}")
        log(f"    H3 Protrusion->Cytotox: r={r3:.4f}, p={p3:.4e} {'PASS' if r3>0 and p3<0.05 else 'FAIL'}")
        log(f"    H4 Topology->Checkpoint: r={r4:.4f}, p={p4:.4e} {'PASS' if r4<0 and p4<0.05 else 'FAIL'}")

    # ---- 3. Tumor-intrinsic SST candidate targets ----
    log("\n--- 3. Tumor-intrinsic SST-axis candidates ---")
    is_other_tumor = (adata.obs["cell_type"] == "Other") & (adata.obs["condition"] == "tumor")
    is_epithelial = pd.Series(False, index=expr.index)
    epi_genes = [g for g in ["EPCAM","KRT19","KRT18","KRT8","CDH1"] if g in expr.columns]
    if epi_genes:
        is_epithelial = expr[epi_genes].mean(axis=1) > expr[epi_genes].mean(axis=1).quantile(0.75)
    is_malignant = is_other_tumor | is_epithelial
    mal_expr = expr.loc[is_malignant]
    log(f"  Malignant cells: {len(mal_expr)}")

    # Score each SST gene for tumor specificity and NK-dysfunction correlation
    candidates = []
    for gene in SST_CORE_GENES:
        if gene not in expr.columns:
            continue
        # Tumor specificity: mean expression in malignant vs non-malignant
        in_mal = np.log1p(mal_expr[gene]).mean()
        in_rest = np.log1p(expr.loc[~is_malignant, gene]).mean()
        tumor_ratio = in_mal / max(in_rest, 0.01)

        # Correlation with NK dysfunction score (across all cells)
        dysf_score = mean_zscore(expr, MODULES["nk_sm_catabolism"] + ["HAVCR2"])
        corr_with_dysf = expr[gene].corr(dysf_score)

        candidates.append({
            "gene": gene, "module": next(n for n, gs in MODULES.items() if gene in gs),
            "tumor_mean_expr": round(in_mal, 3),
            "non_tumor_mean_expr": round(in_rest, 3),
            "tumor_specificity": round(tumor_ratio, 3),
            "nk_dysfunction_correlation": round(corr_with_dysf, 4),
        })
    cand_df = pd.DataFrame(candidates)
    cand_df = cand_df.sort_values("tumor_specificity", ascending=False)
    print("\n" + cand_df.to_string(index=False))
    cand_df.to_csv(os.path.join(out_dir, "sst_axis_gastric_candidates.tsv"), sep="\t", index=False)

    # ---- 4. Summary ----
    log("\n--- 4. Arm B: GASTRIC EXTENSION SUMMARY ---")
    n_axis = sum(1 for _, r in comp_df.iterrows() if r["significant"])
    log(f"  SST score differences (GC vs HL): {n_axis}/{len(metrics)} significant")
    h3_gc_pass = stats.pearsonr(gc_nk["nk_protrusion_machinery_score"], gc_nk["nk_synapse_cytotoxicity_outcome_score"])[1] < 0.05
    log(f"  H3 axis coherence in GC NK: {'PASS' if h3_gc_pass else 'FAIL'}")
    log(f"  Top tumor-specific SST gene: {cand_df.iloc[0]['gene'] if len(cand_df) > 0 else 'N/A'}")
    log(f"  SST-axis candidates: {len(cand_df)} genes")

    log("\n" + "=" * 60)
    log("ARM B: GASTRIC EXTENSION COMPLETE!")
    log("=" * 60)


if __name__ == "__main__":
    main()
