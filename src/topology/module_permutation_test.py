"""
P0-3: Module membership permutation test - stricter null for effector coupling.

Problem: The H3 correlation (protrusion_machinery ~ cytotoxicity_outcome) may
reflect shared NK-activation state rather than a specific protrusion->cytotoxicity
axis coupling. The existing generic 16-gene activation-signature partial
correlation controls for broad activation but does not exclude the possibility
that any two modules drawn from the NK synapse/effector program would show
similar correlation.

Solution: A permutation-based null distribution:
  1. Define a "universe" of NK-expressed genes (detected in >=50% of NK cells).
  2. For each of N=10,000 permutations, randomly draw two modules of the same
     sizes (25 and 11 genes respectively) from this universe.
  3. Compute mean-zscore correlation for each random pair.
  4. Report where the observed r falls in this empirical null distribution
     (empirical P-value, one-tailed: fraction of permuted r >= observed r).

Additionally:
  5. Repeat controlling for generic activation (partial out 16-gene signature
     before permutation, to test whether specific module composition matters
     beyond activation control).

Output:
  results/tables/h3_module_permutation_test.tsv
  results/tables/h3_empirical_null_summary.md

Usage:
    python src/topology/module_permutation_test.py
"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from numpy.random import default_rng

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# Module sizes under test
MODULE_SIZES = {
    "nk_protrusion_machinery": 25,
    "nk_synapse_cytotoxicity_outcome": 11,
}

GENERIC_ACTIVATION = [
    "CD69", "TNF", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5",
    "CSF2", "IL2RA", "ICOS", "TNFSF10", "FASLG", "CD38",
    "HLA-DRA", "HLA-DRB1", "MKI67",
]


def mean_zscore_genes(expr_matrix, gene_indices):
    """Compute mean z-score across selected gene columns (by index)."""
    sub = expr_matrix[:, gene_indices]
    # Per-gene z-score
    gene_means = sub.mean(axis=0)
    gene_stds = sub.std(axis=0, ddof=0)
    gene_stds[gene_stds == 0] = 1.0
    z = (sub - gene_means) / gene_stds
    return z.mean(axis=1)


def partial_corr_via_residuals(x, y, z_matrix):
    """Partial correlation: residualize x and y against z, then correlate."""
    # OLS residualization
    Xz = np.column_stack([np.ones(len(x)), z_matrix])
    beta_x, _, _, _ = np.linalg.lstsq(Xz, x, rcond=None)
    beta_y, _, _, _ = np.linalg.lstsq(Xz, y, rcond=None)
    res_x = x - Xz @ beta_x
    res_y = y - Xz @ beta_y
    return stats.pearsonr(res_x, res_y)


def main():
    log("=" * 70)
    log("P0-3: MODULE MEMBERSHIP PERMUTATION TEST - H3 effector coupling")
    log("=" * 70)

    # ── Load NK expression ──────────────────────────────────────────────
    import scanpy as sc
    import scipy.sparse

    h5ad_path = "data/processed/scrna/gc_integrated.h5ad"
    if not os.path.exists(h5ad_path):
        # Full integrated object is git-ignored (too large to keep locally);
        # the NK-only export carries identical real cells/metadata for this
        # analysis, since it immediately subsets to cell_type == "NK" below.
        real_nk_fallback = "data/processed/scrna/gc_nk_subset_remote.h5ad"
        if os.path.exists(real_nk_fallback):
            log(f"gc_integrated.h5ad not found locally - using real NK-only "
                f"export instead: {real_nk_fallback}")
            h5ad_path = real_nk_fallback
        else:
            log("Real h5ad not found - using synthetic data for demonstration.")
            h5ad_path = "data/synthetic/gc_integrated_synthetic.h5ad"
            if not os.path.exists(h5ad_path):
                log("ERROR: no h5ad available.")
                sys.exit(1)

    adata = sc.read(h5ad_path)
    log(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes   [{h5ad_path}]")

    if scipy.sparse.issparse(adata.X):
        expr = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    else:
        expr = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    if "cell_type" in adata.obs.columns:
        is_nk = adata.obs["cell_type"] == "NK"
    else:
        is_nk = adata.obs.get("nk_subtype", pd.Series(True, index=adata.obs.index)) == "NK"

    nk_expr = expr.loc[is_nk]
    nk_matrix = nk_expr.values
    gene_names = list(nk_expr.columns)
    log(f"  NK cells: {nk_matrix.shape[0]}, genes: {nk_matrix.shape[1]}")

    # ── Define NK-expressed gene universe ───────────────────────────────
    detection_rate = (nk_matrix > 0).mean(axis=0)
    expressed_mask = detection_rate >= 0.5
    expressed_indices = np.where(expressed_mask)[0]
    n_universe = expressed_mask.sum()
    log(f"  Genes expressed in >=50% NK: {n_universe} / {nk_matrix.shape[1]}")

    if n_universe < 50:
        log("  WARNING: too few expressed genes for meaningful permutation test "
            "(synthetic data). Results are illustrative only - re-run on real data.")

    # ── Identify the actual module genes in this dataset ────────────────
    # Build indices for the real modules (genes present in expressed set)
    protrusion_genes = [
        "EZR", "MSN", "RDX", "ACTR2", "ACTR3", "ARPC1B", "ARPC2",
        "ARPC3", "ARPC4", "ARPC5", "WAS", "WASL", "WASF1", "WASF2",
        "WASF3", "WIPF1", "CDC42", "RAC1", "RHOA", "DIAPH1", "DIAPH3",
        "FMNL1", "BAIAP2", "PACSIN2",
    ]
    cytotoxicity_genes = [
        "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "LCP2", "LAT",
        "VAV1", "TLN1", "ITGAL", "ITGB2",
    ]
    activation_gene_set = GENERIC_ACTIVATION

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    protrusion_idx = [gene_to_idx[g] for g in protrusion_genes if g in gene_to_idx]
    cyto_idx = [gene_to_idx[g] for g in cytotoxicity_genes if g in gene_to_idx]
    activation_idx = [gene_to_idx[g] for g in activation_gene_set if g in gene_to_idx]
    log(f"  Protrusion module: {len(protrusion_idx)}/{len(protrusion_genes)} genes found")
    log(f"  Cytotoxicity module: {len(cyto_idx)}/{len(cytotoxicity_genes)} genes found")
    log(f"  Activation control: {len(activation_idx)}/{len(activation_gene_set)} genes found")

    # ── Observed correlation ────────────────────────────────────────────
    obs_x = mean_zscore_genes(nk_matrix, protrusion_idx)
    obs_y = mean_zscore_genes(nk_matrix, cyto_idx)
    obs_r, obs_p = stats.pearsonr(obs_x, obs_y)
    log(f"\n  Observed H3: r={obs_r:.4f}, p={obs_p:.2e}")

    # Partial out activation
    act_mat = nk_matrix[:, activation_idx]
    act_score = act_mat.mean(axis=1)  # simplified activation score
    obs_r_partial, obs_p_partial = partial_corr_via_residuals(obs_x, obs_y, act_score.reshape(-1, 1))
    log(f"  Observed H3 (activation-controlled): r={obs_r_partial:.4f}, p={obs_p_partial:.2e}")

    # ── Permutation test ────────────────────────────────────────────────
    N_PERM = 10000
    rng = default_rng(42)
    n1 = len(protrusion_idx)
    n2 = len(cyto_idx)

    log(f"\n  Running {N_PERM} permutations...")
    log(f"  Each: draw {n1}+{n2} genes from universe of {n_universe} expressed genes")

    perm_r = np.zeros(N_PERM)
    perm_r_partial = np.zeros(N_PERM)

    for i in range(N_PERM):
        if (i + 1) % 2000 == 0:
            log(f"    {i+1}/{N_PERM}...")
        # Randomly draw two disjoint modules
        draw = rng.choice(expressed_indices, size=n1 + n2, replace=False)
        rand_x_idx = draw[:n1]
        rand_y_idx = draw[n1:]

        rand_x = mean_zscore_genes(nk_matrix, rand_x_idx)
        rand_y = mean_zscore_genes(nk_matrix, rand_y_idx)

        # Raw correlation
        perm_r[i], _ = stats.pearsonr(rand_x, rand_y)

        # Activation-controlled
        r_part, _ = partial_corr_via_residuals(rand_x, rand_y, act_score.reshape(-1, 1))
        perm_r_partial[i] = r_part

    # ── Empirical P-values ──────────────────────────────────────────────
    # One-sided: how often does random draw exceed observed r?
    emp_p_raw = np.mean(perm_r >= obs_r)
    emp_p_partial = np.mean(perm_r_partial >= obs_r_partial)

    log(f"\n  =" * 30)
    log(f"  PERMUTATION RESULTS (N={N_PERM})")
    log(f"  =" * 30)
    log(f"  Observed r (raw):                    {obs_r:.4f}")
    log(f"  Null mean r (raw):                   {perm_r.mean():.4f}")
    log(f"  Null SD r (raw):                     {perm_r.std():.4f}")
    log(f"  Null 95th percentile (raw):          {np.percentile(perm_r, 95):.4f}")
    log(f"  Null 99th percentile (raw):          {np.percentile(perm_r, 99):.4f}")
    log(f"  Empirical P (one-sided, raw):         {emp_p_raw:.4f}")
    log(f"  z-score vs null (raw):               {(obs_r - perm_r.mean()) / perm_r.std():.2f}")
    log(f"")
    log(f"  Observed r (activation-controlled):   {obs_r_partial:.4f}")
    log(f"  Null mean r (activation-controlled):  {perm_r_partial.mean():.4f}")
    log(f"  Null SD r (activation-controlled):    {perm_r_partial.std():.4f}")
    log(f"  Null 95th percentile (act-ctrl):      {np.percentile(perm_r_partial, 95):.4f}")
    log(f"  Null 99th percentile (act-ctrl):      {np.percentile(perm_r_partial, 99):.4f}")
    log(f"  Empirical P (one-sided, act-ctrl):    {emp_p_partial:.4f}")
    log(f"  z-score vs null (act-ctrl):           {(obs_r_partial - perm_r_partial.mean()) / perm_r_partial.std():.2f}")

    # ── Save results ────────────────────────────────────────────────────
    out_dir = "results/tables"
    os.makedirs(out_dir, exist_ok=True)

    # Save null distribution plus observed
    out_df = pd.DataFrame({
        "perm_id": range(N_PERM),
        "perm_r_raw": perm_r,
        "perm_r_partial_activation": perm_r_partial,
    })
    out_path = os.path.join(out_dir, "h3_module_permutation_test.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    log(f"\n  Saved: {out_path}")

    # Summary
    summary_path = os.path.join(out_dir, "h3_empirical_null_summary.md")
    with open(summary_path, "w") as f:
        f.write(f"""# H3 Module Membership Permutation Test - Summary

## Method
- **Universe:** {n_universe} genes detected in >=50% of NK cells
- **Permutations:** {N_PERM} random draws of {n1}+{n2}= {n1+n2} genes from the universe
- **Observed modules:** protrusion_machinery ({n1} genes), cytotoxicity_outcome ({n2} genes)
- **Null:** Two random modules of the same sizes, drawn without replacement
- **One-sided test:** fraction of permuted r >= observed r

## Results
| Metric | Raw | Activation-controlled |
|--------|-----|-----------------------|
| Observed r | {obs_r:.4f} | {obs_r_partial:.4f} |
| Null mean r | {perm_r.mean():.4f} | {perm_r_partial.mean():.4f} |
| Null SD | {perm_r.std():.4f} | {perm_r_partial.std():.4f} |
| Null 95th %ile | {np.percentile(perm_r, 95):.4f} | {np.percentile(perm_r_partial, 95):.4f} |
| Null 99th %ile | {np.percentile(perm_r, 99):.4f} | {np.percentile(perm_r_partial, 99):.4f} |
| z-score | {(obs_r - perm_r.mean()) / perm_r.std():.2f} | {(obs_r_partial - perm_r_partial.mean()) / perm_r_partial.std():.2f} |
| Empirical P | {emp_p_raw:.4f} | {emp_p_partial:.4f} |

## Interpretation
- The null mean of random-module correlations ({perm_r.mean():.4f}) quantifies
  the baseline correlation expected from two random {n1}+{n2}-gene sets drawn
  from the NK transcriptome.
- The empirical P-value answers: "given the observed correlation r={obs_r:.4f},
  how often would we see r >= this from {N_PERM} random module pairs?"
- The activation-controlled version partials out a generic 16-gene NK-activation
  signature before computing both observed and permuted correlations.

## Caveats
- Synthetic data results are illustrative only - gene universe is small ({nk_matrix.shape[1]} genes).
- Real data (76+ genes -> ~20,000+ genes) provides a much richer null distribution.
- Re-run on real data for manuscript numbers.
""")
    log(f"  Saved: {summary_path}")

    log("\nP0-3 COMPLETE!")
    log("NOTE: Real data needed for publishable numbers. Synthetic data shows method.")
    log("Expected on real data: null mean near zero (random gene sets uncorrelated),")
    log("observed r >> null 99th %ile -> effector coupling is specific, not generic.")


if __name__ == "__main__":
    main()
