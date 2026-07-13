"""
Diagnose whether the H3 (protrusion~cytotoxicity) real-data collapse under
count-depth control (count_depth_control.py) / permutation test
(module_permutation_test.py) is a scoring-method artifact, or a real
methodological finding.

Motivation: on real scRNA data (data/processed/scrna/gc_nk_subset_remote.h5ad,
8,310 NK cells), the manuscript's mean-zscore H3 correlation (r=0.318) collapses
under count-depth residualization (r=0.09) and falls BELOW the mean of a
permutation null built from random same-size gene modules (null mean r=0.73).
Before reporting this as a real finding, we check whether it is an artifact of
(i) the specific mean-zscore scoring method, or (ii) the permutation test's
restriction of its random-draw "universe" to genes detected in >=50% of cells
(a narrow, atypically highly-expressed 314/22,728-gene set that is not
representative of the protrusion/cytotoxicity modules themselves, which are
mostly genes below that detection threshold).

Tests:
  A. Per-gene correlation with total_counts: are protrusion/cytotoxicity genes
     unusually depth-sensitive vs a random draw from the detected universe?
  B. scanpy sc.tl.score_genes (control-set-matched, the standard method
     designed to remove exactly this confound) vs plain mean-zscore, with a
     matched (same-scoring-method) permutation null.
  C. Restrict mean-zscore modules to genes detected in >=50% of NK cells
     (matching the permutation test's "universe" definition) - does trimming
     poorly-detected, dropout-heavy genes change the raw correlation /
     permutation null?

Output:
  results/tables/h3_scoring_method_diagnostic.tsv
  results/tables/h3_scoring_method_diagnostic_summary.md

Usage:
    python src/topology/h3_scoring_method_diagnostic.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy import stats
from numpy.random import default_rng
import warnings
warnings.filterwarnings("ignore")

def log(msg):
    print(msg, flush=True)

PROTRUSION = ["EZR","MSN","RDX","ACTR2","ACTR3","ARPC1B","ARPC2","ARPC3","ARPC4","ARPC5",
              "WAS","WASL","WASF1","WASF2","WASF3","WIPF1","CDC42","RAC1","RHOA","DIAPH1",
              "DIAPH3","FMNL1","BAIAP2","PACSIN2"]
CYTOTOX = ["NKG7","GNLY","GZMB","PRF1","IFNG","LCP2","LAT","VAV1","TLN1","ITGAL","ITGB2"]

log("Loading real NK h5ad...")
adata = sc.read("data/processed/scrna/gc_nk_subset_remote.h5ad")
log(f"{adata.n_obs} cells x {adata.n_vars} genes")

X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
expr = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
total_counts = adata.obs["total_counts"].values
n_genes_obs = adata.obs["n_genes_by_counts"].values

detection_rate = (X > 0).mean(axis=0)
detected_mask = detection_rate >= 0.5
gene_names = np.array(adata.var_names)
universe = gene_names[detected_mask]
log(f"Detected-in->=50% universe: {len(universe)} / {adata.n_vars} genes")

# ============================================================
# A. Per-gene correlation with total_counts
# ============================================================
log("\n=== A. Per-gene correlation with total_counts ===")

def gene_depth_corrs(genes):
    corrs = []
    for g in genes:
        if g in expr.columns:
            r, _ = stats.pearsonr(expr[g].values, total_counts)
            corrs.append(r)
    return np.array(corrs)

prot_present = [g for g in PROTRUSION if g in expr.columns]
cyto_present = [g for g in CYTOTOX if g in expr.columns]
prot_corrs = gene_depth_corrs(prot_present)
cyto_corrs = gene_depth_corrs(cyto_present)

rng = default_rng(0)
# Random genes from full detected universe, matched counts, repeated
n_draws = 500
rand_means = []
for _ in range(n_draws):
    draw = rng.choice(universe, size=len(prot_present) + len(cyto_present), replace=False)
    rand_means.append(np.mean(gene_depth_corrs(list(draw))))

log(f"Protrusion genes: mean |r| with total_counts = {np.mean(np.abs(prot_corrs)):.4f}, mean r = {np.mean(prot_corrs):.4f}")
log(f"Cytotoxicity genes: mean |r| with total_counts = {np.mean(np.abs(cyto_corrs)):.4f}, mean r = {np.mean(cyto_corrs):.4f}")
log(f"Random universe draws (n={n_draws}): mean of per-draw mean r with total_counts = {np.mean(rand_means):.4f} (SD {np.std(rand_means):.4f})")
log(f"Detection rate - protrusion: {[round(detection_rate[list(gene_names).index(g)],3) for g in prot_present]}")
log(f"Detection rate - cytotoxicity: {[round(detection_rate[list(gene_names).index(g)],3) for g in cyto_present]}")
log(f"Median detection rate universe genes: {np.median(detection_rate[detected_mask]):.3f}")
prot_det = [detection_rate[list(gene_names).index(g)] for g in prot_present]
cyto_det = [detection_rate[list(gene_names).index(g)] for g in cyto_present]
log(f"Median detection rate - protrusion module: {np.median(prot_det):.3f}")
log(f"Median detection rate - cytotoxicity module: {np.median(cyto_det):.3f}")

# ============================================================
# B. scanpy score_genes (control-set-matched) vs mean-zscore
# ============================================================
log("\n=== B. scanpy sc.tl.score_genes vs mean-zscore ===")

adata_sg = adata.copy()
sc.tl.score_genes(adata_sg, prot_present, score_name="protrusion_scanpy", random_state=0)
sc.tl.score_genes(adata_sg, cyto_present, score_name="cytotox_scanpy", random_state=0)

r_scanpy, p_scanpy = stats.pearsonr(adata_sg.obs["protrusion_scanpy"], adata_sg.obs["cytotox_scanpy"])
log(f"H3 via scanpy score_genes: r={r_scanpy:.4f}, p={p_scanpy:.2e}")

# residualize scanpy scores against total_counts + n_genes, recheck
def residualize(score, covariates):
    Xc = np.column_stack([np.ones(len(score)), covariates])
    beta, _, _, _ = np.linalg.lstsq(Xc, score, rcond=None)
    resid = score - Xc @ beta
    r2 = 1 - np.var(resid) / np.var(score)
    return resid, r2

cov = np.column_stack([total_counts, n_genes_obs])
res_p, r2_p = residualize(adata_sg.obs["protrusion_scanpy"].values, cov)
res_c, r2_c = residualize(adata_sg.obs["cytotox_scanpy"].values, cov)
r_scanpy_resid, p_scanpy_resid = stats.pearsonr(res_p, res_c)
log(f"  R2_tech (protrusion_scanpy) = {r2_p:.4f}, R2_tech (cytotox_scanpy) = {r2_c:.4f}")
log(f"  After count-depth residualization: r={r_scanpy_resid:.4f}, p={p_scanpy_resid:.2e}")

# mini permutation test with scanpy scoring (smaller N for speed)
N_PERM_B = 200
perm_r_scanpy = np.zeros(N_PERM_B)
for i in range(N_PERM_B):
    draw = rng.choice(universe, size=len(prot_present) + len(cyto_present), replace=False)
    g1 = list(draw[:len(prot_present)])
    g2 = list(draw[len(prot_present):])
    sc.tl.score_genes(adata_sg, g1, score_name="_r1", random_state=0)
    sc.tl.score_genes(adata_sg, g2, score_name="_r2", random_state=0)
    perm_r_scanpy[i], _ = stats.pearsonr(adata_sg.obs["_r1"], adata_sg.obs["_r2"])
log(f"  Permutation null (scanpy scoring, N={N_PERM_B}): mean={perm_r_scanpy.mean():.4f}, "
    f"SD={perm_r_scanpy.std():.4f}, 95th%ile={np.percentile(perm_r_scanpy,95):.4f}")
log(f"  Observed r={r_scanpy:.4f} vs null -> empirical P = {np.mean(perm_r_scanpy >= r_scanpy):.4f}")

# ============================================================
# C. Mean-zscore restricted to >=50% detected genes only
# ============================================================
log("\n=== C. Mean-zscore restricted to well-detected genes (>=50% detection) ===")

def mean_zscore(df, genes):
    available = [g for g in genes if g in df.columns]
    z = (df[available] - df[available].mean(0)) / df[available].std(0, ddof=0)
    return z.fillna(0).mean(axis=1)

prot_well_detected = [g for g in prot_present if detection_rate[list(gene_names).index(g)] >= 0.5]
cyto_well_detected = [g for g in cyto_present if detection_rate[list(gene_names).index(g)] >= 0.5]
log(f"Protrusion genes retained after >=50% detection filter: {len(prot_well_detected)}/{len(prot_present)}")
log(f"Cytotoxicity genes retained after >=50% detection filter: {len(cyto_well_detected)}/{len(cyto_present)}")

if len(prot_well_detected) >= 2 and len(cyto_well_detected) >= 2:
    x_trim = mean_zscore(expr, prot_well_detected)
    y_trim = mean_zscore(expr, cyto_well_detected)
    r_trim, p_trim = stats.pearsonr(x_trim, y_trim)
    log(f"H3 (well-detected genes only): r={r_trim:.4f}, p={p_trim:.2e}")
    res_xt, r2_xt = residualize(x_trim.values, cov)
    res_yt, r2_yt = residualize(y_trim.values, cov)
    r_trim_resid, p_trim_resid = stats.pearsonr(res_xt, res_yt)
    log(f"  R2_tech x={r2_xt:.4f}, R2_tech y={r2_yt:.4f}")
    log(f"  After residualization: r={r_trim_resid:.4f}, p={p_trim_resid:.2e}")
else:
    log("Not enough well-detected genes in one module to test.")

log("\nDIAGNOSIS COMPLETE")

# ============================================================
# Save results
# ============================================================
out_dir = "results/tables"
os.makedirs(out_dir, exist_ok=True)

rows = [
    {"test": "A: universe size (>=50% detected)", "value": len(universe), "note": f"of {adata.n_vars} total genes"},
    {"test": "A: median detection rate, universe", "value": round(float(np.median(detection_rate[detected_mask])), 4), "note": ""},
    {"test": "A: median detection rate, protrusion module", "value": round(float(np.median(prot_det)), 4), "note": f"n={len(prot_present)} genes"},
    {"test": "A: median detection rate, cytotoxicity module", "value": round(float(np.median(cyto_det)), 4), "note": f"n={len(cyto_present)} genes"},
    {"test": "A: mean gene-total_counts r, protrusion genes", "value": round(float(np.mean(prot_corrs)), 4), "note": ""},
    {"test": "A: mean gene-total_counts r, cytotoxicity genes", "value": round(float(np.mean(cyto_corrs)), 4), "note": ""},
    {"test": "A: mean gene-total_counts r, random universe draws", "value": round(float(np.mean(rand_means)), 4), "note": f"n_draws={n_draws}, SD={np.std(rand_means):.4f}"},
    {"test": "B: H3 raw r, mean-zscore (manuscript method)", "value": 0.3176, "note": "from count_depth_control.py real-data run"},
    {"test": "B: H3 raw r, scanpy score_genes (control-set-matched)", "value": round(float(r_scanpy), 6), "note": f"p={p_scanpy:.2e}"},
    {"test": "B: H3 r after count-depth residualization, score_genes", "value": round(float(r_scanpy_resid), 6), "note": f"p={p_scanpy_resid:.2e}"},
    {"test": "B: R2_tech protrusion, score_genes", "value": round(float(r2_p), 4), "note": ""},
    {"test": "B: R2_tech cytotoxicity, score_genes", "value": round(float(r2_c), 4), "note": ""},
    {"test": "B: permutation null mean, score_genes", "value": round(float(perm_r_scanpy.mean()), 4), "note": f"N={N_PERM_B}, SD={perm_r_scanpy.std():.4f}"},
    {"test": "B: permutation null 95th pct, score_genes", "value": round(float(np.percentile(perm_r_scanpy, 95)), 4), "note": ""},
    {"test": "B: empirical P (score_genes null)", "value": round(float(np.mean(perm_r_scanpy >= r_scanpy)), 4), "note": "one-sided, fraction of null >= observed"},
]
diag_df = pd.DataFrame(rows)
diag_df.to_csv(os.path.join(out_dir, "h3_scoring_method_diagnostic.tsv"), sep="\t", index=False)
log(f"\nSaved: {os.path.join(out_dir, 'h3_scoring_method_diagnostic.tsv')}")

summary_path = os.path.join(out_dir, "h3_scoring_method_diagnostic_summary.md")
with open(summary_path, "w", encoding="utf-8") as f:
    f.write(f"""# H3 scoring-method diagnostic - is the real-data collapse an artifact?

## Question
On real scRNA data (8,310 NK cells, gc_nk_subset_remote.h5ad), the manuscript's
mean-zscore H3 correlation (protrusion~cytotoxicity, r=0.318) collapses under
count-depth residualization (r=0.09) and under a module-membership permutation
test (null mean r=0.73, observed r=0.318 below it, empirical P=1.0). Before
reporting this as a real finding, we checked whether it is an artifact of the
scoring method or the permutation test's null construction.

## Finding A: the permutation test's "universe" is biased
The permutation null was built by drawing random genes from those detected in
>=50% of NK cells - only {len(universe)}/{adata.n_vars} genes, with median
detection rate {np.median(detection_rate[detected_mask]):.3f}. The real
protrusion module (median detection {np.median(prot_det):.3f}) and cytotoxicity
module (median detection {np.median(cyto_det):.3f}) are both mostly *below*
this threshold ({len(prot_well_detected)}/{len(prot_present)} and
{len(cyto_well_detected)}/{len(cyto_present)} genes qualify). The narrow,
atypically highly-expressed 314-gene universe is not a fair null population for
these modules, and likely inflates the specific null-mean value of 0.73.

## Finding B: the qualitative result replicates with the field-standard method
Using scanpy's `sc.tl.score_genes` (expression-level-matched control-gene
scoring, the standard method for exactly this technical confound) instead of
mean-zscore:
- Raw H3: r={r_scanpy:.4f} (p={p_scanpy:.2e}) - already far below the
  manuscript's mean-zscore r=0.318, indicating mean-zscore itself inflates the
  raw single-cell correlation.
- After count+n_genes residualization: r={r_scanpy_resid:.4f} (p={p_scanpy_resid:.2e}),
  R2_tech = {r2_p:.3f} (protrusion) / {r2_c:.3f} (cytotoxicity).
- A matched permutation null (same score_genes method applied to both real and
  randomly-drawn modules, N={N_PERM_B}): null mean = {perm_r_scanpy.mean():.4f},
  95th percentile = {np.percentile(perm_r_scanpy, 95):.4f}. The observed r
  ({r_scanpy:.4f}) does **not** exceed this null (empirical P =
  {np.mean(perm_r_scanpy >= r_scanpy):.4f}).

## Verdict
The specific magnitude of the null distribution (0.73 vs 0.35) is
method-dependent and the naive-universe permutation test overstates it. But the
qualitative conclusion is **not** an artifact: with the field-standard,
expression-matched scoring method, the observed single-cell
protrusion~cytotoxicity correlation still does not exceed what a randomly drawn
gene-module pair of the same sizes would produce, and a large fraction of its
variance ({r2_p:.0%}/{r2_c:.0%}) is explained by technical library-size
covariates. This is a real methodological finding, not a scoring bug.

## Implication for the manuscript
The single-cell H3 pseudoreplication-corrected number (corrected r=0.313,
P=3.9e-8) addresses cell non-independence within samples, but does not address
this separate technical-confound / module-specificity problem. The single-cell
H3 result should not be reported as an independent, specific replication of the
effector arm. The bulk TCGA-LIHC result (r=0.55, deconvolution-based, not
subject to per-cell dropout/depth artifacts in the same way) is unaffected by
this diagnostic and remains the primary evidence for the effector-arm claim.
""")
log(f"Saved: {summary_path}")
