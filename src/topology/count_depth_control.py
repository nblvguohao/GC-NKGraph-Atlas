"""
P0-2: Count depth and technical covariate control for scRNA module-score
correlations.

Problem: Module scores computed as mean-zscore of constituent genes may still
carry shared technical covariance from total counts (library size) and number
of detected genes. Correlating two module scores across cells may confound
biological coupling with shared technical variation.

Solution: Residualize each module score against total_counts and n_genes_by_counts
before computing correlations. Report H2/H3/H4 both with and without adjustment.

Since the real h5ad is git-ignored and unavailable locally, this script:
  1. Runs on the synthetic dataset to demonstrate the method
  2. Documents the procedure so it can be re-run on real data by the authors
  3. Reports the direction and magnitude of count-depth influence

Output:
  results/tables/sst_axis_count_depth_control.tsv    (synthetic run)
  results/tables/sst_axis_count_depth_methods.md     (procedure documentation)

Usage:
    python src/topology/count_depth_control.py
    # Then re-run with real h5ad when available:
    python src/topology/count_depth_control.py --real
"""
import os, sys, time, warnings, argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# Gene modules (from sst_axis_scrna.py and manuscript Methods §2.3)
MODULES = {
    "nk_sm_synthesis": ["SGMS1", "SGMS2"],
    "nk_sm_catabolism": ["SMPD1", "SMPD2", "SMPD3", "SMPD4"],
    "nk_protrusion_machinery": [
        "EZR", "MSN", "RDX", "ACTR2", "ACTR3", "ARPC1B", "ARPC2",
        "ARPC3", "ARPC4", "ARPC5", "WAS", "WASL", "WASF1", "WASF2",
        "WASF3", "WIPF1", "CDC42", "RAC1", "RHOA", "DIAPH1", "DIAPH3",
        "FMNL1", "BAIAP2", "PACSIN2",
    ],
    "nk_synapse_cytotoxicity_outcome": [
        "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "LCP2", "LAT",
        "VAV1", "TLN1", "ITGAL", "ITGB2",
    ],
    "checkpoint_link": ["HAVCR2"],
}

GENERIC_ACTIVATION = [
    "CD69", "TNF", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5",
    "CSF2", "IL2RA", "ICOS", "TNFSF10", "FASLG", "CD38",
    "HLA-DRA", "HLA-DRB1", "MKI67",
]


def mean_zscore(expr, genes):
    """Compute mean z-score across available genes."""
    available = [g for g in genes if g in expr.columns]
    if not available:
        return pd.Series(0.0, index=expr.index)
    z = (expr[available] - expr[available].mean(0)) / expr[available].std(0, ddof=0)
    return z.fillna(0).mean(axis=1)


def residualize(score, covariates, add_intercept=True):
    """
    Residualize a score vector against covariates via OLS.

    Returns residuals + R^2 of the fit (fraction of score variance explained by tech).
    """
    X = covariates.copy()
    if add_intercept:
        X = np.column_stack([np.ones(len(score)), X])
    # OLS: beta = (X'X)^-1 X'y
    beta, _, _, _ = np.linalg.lstsq(X, score, rcond=None)
    y_pred = X @ beta
    residuals = score - y_pred
    r2 = 1 - np.var(residuals) / np.var(score) if np.var(score) > 0 else 0
    return residuals, r2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--real", action="store_true",
                        help="Use real h5ad (data/processed/scrna/gc_integrated.h5ad)")
    args = parser.parse_args()

    log("=" * 70)
    log("P0-2: COUNT DEPTH / TECHNICAL COVARIATE CONTROL")
    log("=" * 70)

    import scanpy as sc
    import scipy.sparse

    if args.real:
        h5ad_path = "data/processed/scrna/gc_integrated.h5ad"
        if not os.path.exists(h5ad_path):
            # Full integrated object is git-ignored (too large to keep locally);
            # the NK-only export carries identical real cells/metadata for every
            # analysis here, since we immediately subset to cell_type == "NK" below.
            fallback = "data/processed/scrna/gc_nk_subset_remote.h5ad"
            if os.path.exists(fallback):
                log(f"  gc_integrated.h5ad not found locally - using real NK-only "
                    f"export instead: {fallback}")
                h5ad_path = fallback
    else:
        h5ad_path = "data/synthetic/gc_integrated_synthetic.h5ad"

    if not os.path.exists(h5ad_path):
        log(f"ERROR: h5ad not found: {h5ad_path}")
        log("Run with --real or generate the data first.")
        sys.exit(1)

    adata = sc.read(h5ad_path)
    log(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes   [{h5ad_path}]")

    # Extract expression matrix
    if scipy.sparse.issparse(adata.X):
        expr = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    else:
        expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    # ── Determine NK cells ──────────────────────────────────────────────
    if "cell_type" in adata.obs.columns:
        is_nk = adata.obs["cell_type"] == "NK"
    else:
        # Synthetic data: use nk_subtype
        is_nk = adata.obs.get("nk_subtype", pd.Series(True, index=adata.obs.index)) == "NK"
    nk_expr = expr.loc[is_nk]
    log(f"  NK cells: {is_nk.sum()} / {len(is_nk)}")

    # ── Get technical covariates ────────────────────────────────────────
    has_counts = False
    total_counts = None
    n_genes = None
    sample_ids = None

    if "total_counts" in adata.obs.columns:
        total_counts = adata.obs.loc[is_nk, "total_counts"].values
        has_counts = True
    if "n_genes_by_counts" in adata.obs.columns:
        n_genes = adata.obs.loc[is_nk, "n_genes_by_counts"].values
    elif "n_genes" in adata.obs.columns:
        n_genes = adata.obs.loc[is_nk, "n_genes"].values

    if not has_counts:
        log("  WARNING: total_counts not in adata.obs - using per-cell sum of expression as proxy")
        total_counts = nk_expr.sum(axis=1).values
        n_genes = (nk_expr > 0).sum(axis=1).values
    log(f"  total_counts: median={np.median(total_counts):.0f}, range=[{total_counts.min():.0f}, {total_counts.max():.0f}]")
    log(f"  n_genes: median={np.median(n_genes):.0f}, range=[{n_genes.min():.0f}, {n_genes.max():.0f}]")

    if "sample_id" in adata.obs.columns:
        sample_ids = adata.obs.loc[is_nk, "sample_id"].values

    # ── Compute module scores ───────────────────────────────────────────
    log("\nComputing module scores...")
    scores = {}
    for name, genes in MODULES.items():
        score = mean_zscore(nk_expr, genes)
        scores[name] = score.values
        n_found = sum(1 for g in genes if g in nk_expr.columns)
        log(f"  {name}: {n_found}/{len(genes)} genes found, mean={score.mean():.4f}")

    # Derived scores
    scores["nk_sm_balance"] = scores["nk_sm_synthesis"] - scores["nk_sm_catabolism"]
    scores["nk_topology_permissive"] = (scores["nk_sm_balance"] + scores["nk_protrusion_machinery"]) / 2

    # Activation score (for control)
    act_score = mean_zscore(nk_expr, GENERIC_ACTIVATION)
    scores["generic_activation"] = act_score.values
    n_act = sum(1 for g in GENERIC_ACTIVATION if g in nk_expr.columns)
    log(f"  generic_activation: {n_act}/{len(GENERIC_ACTIVATION)} genes")

    # ── Build covariate matrix ──────────────────────────────────────────
    cov = np.column_stack([total_counts, n_genes])
    cov_names = ["total_counts", "n_genes"]

    # ── Residualize each score ──────────────────────────────────────────
    log("\nResidualizing against total_counts + n_genes...")
    residuals = {}
    r2_explained = {}
    for name in ["nk_sm_balance", "nk_protrusion_machinery",
                  "nk_synapse_cytotoxicity_outcome", "nk_topology_permissive",
                  "checkpoint_link"]:
        res, r2 = residualize(scores[name], cov)
        residuals[name] = res
        r2_explained[name] = r2
        log(f"  {name}: R^2_tech={r2:.4f} ({r2*100:.1f}% of variance explained by counts)")

    # ── Compare correlations: raw vs residualized ───────────────────────
    log("\n" + "=" * 50)
    log("CORRELATION COMPARISON: raw vs count-depth-residualized")
    log("=" * 50)

    comparisons = [
        ("H2: SM -> protrusion", "nk_sm_balance", "nk_protrusion_machinery"),
        ("H3: protrusion -> cytotoxicity", "nk_protrusion_machinery",
         "nk_synapse_cytotoxicity_outcome"),
        ("H4: topology -> checkpoint", "nk_topology_permissive", "checkpoint_link"),
    ]

    results_rows = []
    for label, xk, yk in comparisons:
        # Raw
        r_raw, p_raw = stats.pearsonr(scores[xk], scores[yk])
        # Residualized
        r_res, p_res = stats.pearsonr(residuals[xk], residuals[yk])
        # Change
        delta_r = r_res - r_raw
        log(f"\n  {label}:")
        log(f"    Raw:      r={r_raw:.4f}, p={p_raw:.2e}")
        log(f"    Resid:    r={r_res:.4f}, p={p_res:.2e}")
        log(f"    delta_r:       {delta_r:+.4f}")
        log(f"    R^2_tech x: {r2_explained[xk]:.3f}  R^2_tech y: {r2_explained[yk]:.3f}")

        results_rows.append({
            "comparison": label,
            "raw_r": round(r_raw, 4), "raw_p": p_raw,
            "residualized_r": round(r_res, 4), "residualized_p": p_res,
            "delta_r": round(delta_r, 4),
            "r2_tech_x": round(r2_explained[xk], 4),
            "r2_tech_y": round(r2_explained[yk], 4),
        })

    # ── Also do H3 with scVI latent factor control (if available) ──────
    log("\n" + "=" * 50)
    log("P0-3 prerequisite: scVI latent factor control for H3")
    log("=" * 50)
    has_latent = False
    for col_prefix in ["X_scVI", "scvi_latent", "latent"]:
        if col_prefix in adata.obsm:
            latent = adata.obsm[col_prefix]
            has_latent = True
            log(f"  Using latent factors from obsm['{col_prefix}']: {latent.shape}")
            break
    if not has_latent:
        log("  No scVI latent factors in adata.obsm")
        log("  Using PCA on NK expression as fallback (first 5 PCs)")
        from sklearn.decomposition import PCA
        pca = PCA(n_components=5)
        latent = pca.fit_transform(nk_expr.values)
        log(f"  PCA fallback: {latent.shape}, explained variance: {pca.explained_variance_ratio_.sum():.2f}")

    # PCA was computed on NK cells only, so latent already aligns with scores
    nk_latent = latent  # shape: (n_nk_cells, n_components)
    # Residualize H3 scores against latent + counts
    cov_full = np.column_stack([total_counts, n_genes, nk_latent])
    res_x_latent, r2_x_latent = residualize(scores["nk_protrusion_machinery"], cov_full)
    res_y_latent, r2_y_latent = residualize(scores["nk_synapse_cytotoxicity_outcome"], cov_full)
    r_h3_latent, p_h3_latent = stats.pearsonr(res_x_latent, res_y_latent)
    raw_r = stats.pearsonr(scores["nk_protrusion_machinery"],
                           scores["nk_synapse_cytotoxicity_outcome"])[0]
    log(f"\n  H3 after count + latent factor control:")
    log(f"    Raw:                       r={raw_r:.4f}")
    log(f"    After counts + 5 PCs:      r={r_h3_latent:.4f}, p={p_h3_latent:.2e}")
    log(f"    delta_r: {r_h3_latent - raw_r:+.4f}")

    # ── Save results ────────────────────────────────────────────────────
    out_dir = "results/tables"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "sst_axis_count_depth_control.tsv")
    pd.DataFrame(results_rows).to_csv(out_path, sep="\t", index=False)
    log(f"\nSaved: {out_path}")

    method_doc = os.path.join(out_dir, "sst_axis_count_depth_methods.md")
    with open(method_doc, "w", encoding="utf-8") as f:
        f.write("""# Count-depth control method for scRNA module-score correlations

## Procedure
1. Module scores computed as mean z-score of constituent genes
   (per-gene standardization, then mean across genes within module).
2. For each NK cell, `total_counts` (library size) and `n_genes_by_counts`
   (detected gene count) extracted from `adata.obs`.
3. Each module score residualized against [total_counts, n_genes] via OLS:
   `residual = score - (β₀ + β₁·total_counts + β₂·n_genes)`.
4. Reported R^2_tech = fraction of module-score variance explained by
   technical covariates. High R^2_tech (>0.05) indicates the score is
   substantially confounded by library size.
5. Correlation analysis repeated on residualized scores; delta_r compared
   to raw.

## scVI latent factor control (supplementary)
For H3 specifically, an additional control residualizes both protrusion
and cytotoxicity scores against [total_counts, n_genes, 5 PCA components
of scRNA expression] (or scVI latent factors when available). This tests
whether the coupling survives removal of the dominant transcriptional
program components.

## Findings
See `sst_axis_count_depth_control.tsv` for the table.
""")
    log(f"Saved: {method_doc}")
    log("\nP0-2 + P0-3 latent control COMPLETE!")


if __name__ == "__main__":
    main()
