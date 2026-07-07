"""
GC-NKGraph-Atlas SST Axis scRNA analysis (Phase 14R).
Cell-type-specific SST scoring on full integrated data.
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
EPITHELIAL_MARKERS = ["EPCAM", "KRT19", "KRT18", "KRT8", "CDH1"]

def mean_zscore(expr, genes):
    available = [g for g in genes if g in expr.columns]
    if not available:
        return pd.Series(0.0, index=expr.index)
    z = (expr[available] - expr[available].mean(0)) / expr[available].std(0, ddof=0)
    return z.fillna(0).mean(axis=1)

def main():
    log("=" * 60)
    log("SST AXIS scRNA -- Cell-Type-Specific Analysis")
    log("=" * 60)
    out_dir = "results/tables"
    fig_dir = "results/figures"
    integrated_path = "data/processed/scrna/gc_integrated.h5ad"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    import scanpy as sc
    import scipy.sparse
    adata = sc.read(integrated_path)
    log(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    if scipy.sparse.issparse(adata.X):
        expr = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    else:
        expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    log(f"  Cell types: {dict(adata.obs['cell_type'].value_counts())}")

    is_nk = adata.obs["cell_type"] == "NK"
    is_tumor = adata.obs["condition"] == "tumor"
    is_other = adata.obs["cell_type"] == "Other"
    epi_markers = [g for g in EPITHELIAL_MARKERS if g in expr.columns]
    is_epithelial = pd.Series(False, index=expr.index)
    if epi_markers:
        epi_score = expr[epi_markers].mean(axis=1)
        is_epithelial = epi_score > epi_score.quantile(0.75)
    is_malignant = (is_other & is_tumor) | is_epithelial
    log(f"  NK: {is_nk.sum()}, Malignant proxy: {is_malignant.sum()}")

    # NK scores
    nk_expr = expr.loc[is_nk]
    nk_scores = pd.DataFrame(index=nk_expr.index)
    for name in ["nk_sm_synthesis", "nk_sm_catabolism", "nk_protrusion_machinery",
                  "nk_synapse_cytotoxicity_outcome", "checkpoint_link"]:
        nk_scores[f"{name}_score"] = mean_zscore(nk_expr, MODULES[name])
        n_found = sum(1 for g in MODULES[name] if g in nk_expr.columns)
        log(f"  NK {name}: {n_found}/{len(MODULES[name])} genes")

    nk_scores["nk_sm_balance_score"] = nk_scores["nk_sm_synthesis_score"] - nk_scores["nk_sm_catabolism_score"]
    nk_scores["nk_topology_permissive_score"] = (nk_scores["nk_sm_balance_score"] + nk_scores["nk_protrusion_machinery_score"]) / 2
    nk_scores["tissue"] = adata.obs.loc[is_nk, "tissue"].values
    nk_scores["sample_id"] = adata.obs.loc[is_nk, "sample_id"].values

    # Tumor scores
    if is_malignant.sum() > 0:
        tumor_expr = expr.loc[is_malignant]
        tumor_serine = mean_zscore(tumor_expr, MODULES["tumor_serine_capacity"])
        n_found = sum(1 for g in MODULES["tumor_serine_capacity"] if g in tumor_expr.columns)
        log(f"  Tumor serine: {n_found}/{len(MODULES['tumor_serine_capacity'])} genes")

        # H1: Crosstalk
        tumor_mean = tumor_serine.groupby(adata.obs.loc[is_malignant, "sample_id"]).mean()
        nk_bal = nk_scores["nk_sm_balance_score"].groupby(nk_scores["sample_id"]).mean()
        common = tumor_mean.index.intersection(nk_bal.index)
        if len(common) >= 3:
            r, p = stats.pearsonr(tumor_mean.loc[common], nk_bal.loc[common])
            log(f"  H1: r={r:.4f}, p={p:.4e}, sign={'POS' if r>0 else 'NEG'} CALIBRATED")
        else:
            log(f"  H1: insufficient samples ({len(common)})")
    else:
        log("  No malignant cells, skipping H1")

    # H2-H4 within NK
    r, p = stats.pearsonr(nk_scores["nk_sm_balance_score"], nk_scores["nk_protrusion_machinery_score"])
    log(f"  H2: SM->Protrusion r={r:.4f}, p={p:.4e} {'PASS' if r>0 and p<0.05 else 'FAIL'}")

    r, p = stats.pearsonr(nk_scores["nk_protrusion_machinery_score"], nk_scores["nk_synapse_cytotoxicity_outcome_score"])
    log(f"  H3: Protrusion->Cytotox r={r:.4f}, p={p:.4e} {'PASS' if r>0 and p<0.05 else 'FAIL'}")

    r, p = stats.pearsonr(nk_scores["nk_topology_permissive_score"], nk_scores["checkpoint_link_score"])
    log(f"  H4: Topology->Checkpoint r={r:.4f}, p={p:.4e} {'PASS' if r<0 and p<0.05 else 'FAIL'}")

    # Tissue comparison
    log("\nSST scores by tissue:")
    tissue_cols = ["nk_sm_balance_score", "nk_protrusion_machinery_score",
                    "nk_topology_permissive_score", "nk_synapse_cytotoxicity_outcome_score"]
    ts = nk_scores.groupby("tissue")[tissue_cols].mean().round(4)
    print(ts.to_string())
    ts.to_csv(os.path.join(out_dir, "sst_axis_scrna_by_tissue.tsv"), sep="\t")
    nk_scores.to_csv(os.path.join(out_dir, "sst_axis_scrna_nk_scores.tsv"), sep="\t")
    log(f"  Saved: {len(nk_scores)} NK cells")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        for ax, col in zip(axes.flatten(), tissue_cols):
            for tissue in nk_scores["tissue"].unique():
                data = nk_scores[nk_scores["tissue"] == tissue][col]
                ax.hist(data, bins=50, alpha=0.5, label=tissue, density=True)
            ax.set_title(col)
            ax.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "fig9_sst_axis_scrna.pdf"), dpi=150)
        log("  Figure saved")
    except Exception as e:
        log(f"  Plot failed: {e}")

    log("\nSST scRNA ANALYSIS COMPLETE!")

if __name__ == "__main__":
    main()
