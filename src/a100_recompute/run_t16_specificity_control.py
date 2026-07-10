"""
T16 — Effector-arm MECHANISM-SPECIFICITY control (peer-review R1/R3).

Rationale (from the Nature-style pre-submission review):
  T14 already showed H3 (protrusion-machinery ~ cytotoxicity-output) survives a
  partial correlation controlling for a 16-gene NK-activation signature
  (partial r = 0.286). Reviewers noted this answers "is the coupling orthogonal
  to a generic activation axis?" but NOT "is the coupling SPECIFIC to the
  serine->SM->topology mechanism, rather than a property of ANY two
  activation-correlated NK modules of the same size?"

Approach (matched-null permutation):
  1. Load the real protrusion & cytotoxicity module gene lists from the canonical
     SST config (src/common/sst_config.py — single source of truth).
  2. Recompute the OBSERVED effect: partial r of protrusion~cytotoxicity |
     activation, exactly as in T14, so the observed value is self-consistent.
  3. Build a NULL POOL of NK-expressed genes EXCLUDING every SST-axis gene and
     every activation gene.
  4. Draw N>=1000 random "pseudo-protrusion" (size = |protrusion present|) and
     "pseudo-cytotoxicity" (size = |cytotoxicity present|) module pairs, disjoint,
     MATCHED on activation loading (each pseudo-module's |corr to activation| is
     required to fall within a tolerance band around the real module's loading).
  5. For each pair compute the same partial r | activation -> null distribution.
  6. Empirical p = (1 + #{null >= observed}) / (1 + N). If p < 0.05 the observed
     coupling is stronger than matched random NK module pairs -> MECHANISM-SPECIFIC.
     Otherwise the "effector arm recovers" wording is downgraded to
     "consistent with, but not specific beyond, generic activation-module covariance".

Input:
  data/processed/scrna/gc_nk_subset.h5ad  (8,310 NK x 22,728 genes, log-normalized;
                                           the same gene-level file T14 used)

Produces:
  results/tables/t16_specificity_control.tsv      (observed vs null summary, verdict)
  results/tables/t16_null_distribution.tsv        (per-permutation null r; feeds §3.2 suppl. fig)

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_t16_specificity_control.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from scipy.stats import linregress

warnings.filterwarnings("ignore")

# Make src importable when run from project root or from src/a100_recompute
_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.common.sst_config import load_sst_modules, get_sst_genes  # noqa: E402

T = "results/tables/"
H5AD = "data/processed/scrna/gc_nk_subset.h5ad"
N_PERM = 1000
RNG = np.random.default_rng(42)          # reproducible, same seed convention as the folds
ACT_TOL = 0.5                            # matched-null tolerance: |corr-to-activation| within +/-50% of real
MIN_DETECT_FRAC = 0.05                  # a null-pool gene must be detected in >=5% of NK cells

# --- Generic NK-activation signature (mirrors T14 run_h3_activation_control.py) ---
# Kept in sync with T14 by hand rather than imported: that module runs its analysis
# at import time, so importing it would execute the whole T14 script.
NK_ACTIVATION_GENES = [
    "CD69", "TNF", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5", "CSF2",
    "IL2RA", "ICOS", "TNFSF10", "FASLG", "CD38", "HLA-DRA", "HLA-DRB1", "MKI67",
]


def zscore_cols(mat: np.ndarray) -> np.ndarray:
    """Column-wise (per-gene) z-score of a cells x genes matrix."""
    mu = mat.mean(axis=0)
    sd = mat.std(axis=0)
    return (mat - mu) / (sd + 1e-10)


def partial_r(x: np.ndarray, y: np.ndarray, a: np.ndarray):
    """Partial Pearson r of x~y controlling for a (residualize both on a)."""
    sx, ix, *_ = linregress(a, x)
    sy, iy, *_ = linregress(a, y)
    xr = x - (sx * a + ix)
    yr = y - (sy * a + iy)
    return stats.pearsonr(xr, yr)


def main() -> int:
    # --- 1. Module gene lists from the canonical config ---
    modules = load_sst_modules()
    if "nk_protrusion_machinery" not in modules or "nk_synapse_cytotoxicity_outcome" not in modules:
        print("ERROR: expected modules missing from sst_axis_config.yaml")
        print("Available:", list(modules))
        return 1
    prot_genes = list(modules["nk_protrusion_machinery"]["genes"])
    cyto_genes = list(modules["nk_synapse_cytotoxicity_outcome"]["genes"])
    sst_all = set(get_sst_genes(modules))
    print(f"Protrusion module: {len(prot_genes)} genes; cytotoxicity module: {len(cyto_genes)} genes")

    # --- 2. Load gene-level NK expression from the h5ad (same file T14 used) ---
    import anndata as ad
    adata = ad.read_h5ad(H5AD)
    if "cell_type" in adata.obs.columns:
        adata = adata[adata.obs["cell_type"] == "NK"].copy()
    var = list(map(str, adata.var_names))
    var_idx = {g: i for i, g in enumerate(var)}
    X = adata.X
    X = X.tocsr() if sp.issparse(X) else np.asarray(X)
    print(f"scNK cells: {adata.n_obs}; genes: {adata.n_vars}")

    prot_present = [g for g in prot_genes if g in var_idx]
    cyto_present = [g for g in cyto_genes if g in var_idx]
    act_present = [g for g in NK_ACTIVATION_GENES if g in var_idx]
    print(f"Present -> protrusion {len(prot_present)}/{len(prot_genes)}, "
          f"cytotoxicity {len(cyto_present)}/{len(cyto_genes)}, activation {len(act_present)}")
    if len(prot_present) < 3 or len(cyto_present) < 3 or len(act_present) < 2:
        print("ERROR: too few genes present to run the specificity control")
        return 1

    def dense(genes):
        idx = [var_idx[g] for g in genes]
        sub = X[:, idx]
        return sub.toarray() if sp.issparse(sub) else np.asarray(sub)

    # --- 3. Activation score + observed partial r (self-consistent with T14) ---
    x_obs = zscore_cols(dense(prot_present)).mean(axis=1)
    y_obs = zscore_cols(dense(cyto_present)).mean(axis=1)
    activation = zscore_cols(dense(act_present)).mean(axis=1)
    valid = ~(np.isnan(x_obs) | np.isnan(y_obs) | np.isnan(activation))
    x_obs, y_obs, a = x_obs[valid], y_obs[valid], activation[valid]
    n_valid = len(a)

    r_raw, _ = stats.pearsonr(x_obs, y_obs)
    r_obs, p_obs = partial_r(x_obs, y_obs, a)
    # Real modules' activation loading (target for the matched null)
    load_prot = abs(stats.pearsonr(x_obs, a)[0])
    load_cyto = abs(stats.pearsonr(y_obs, a)[0])
    print(f"\nOBSERVED  raw r={r_raw:.4f}  partial r|act={r_obs:.4f} (p={p_obs:.2e})  n={n_valid}")
    print(f"Activation loading  protrusion={load_prot:.3f}  cytotoxicity={load_cyto:.3f}")

    # --- 4. Build null pool: NK-expressed genes, excluding SST + activation ---
    exclude = sst_all | set(NK_ACTIVATION_GENES)
    detect = np.asarray((X > 0).mean(axis=0)).ravel()
    pool = [i for i, g in enumerate(var) if g not in exclude and detect[i] >= MIN_DETECT_FRAC]
    print(f"Null pool: {len(pool)} NK-expressed genes (excluded {len(exclude)} SST/activation genes)")
    if len(pool) < (len(prot_present) + len(cyto_present)) * 3:
        print("ERROR: null pool too small for disjoint sampling")
        return 1

    pool_mat = X[:, pool]
    pool_mat = (pool_mat.toarray() if sp.issparse(pool_mat) else np.asarray(pool_mat))[valid]
    pool_z = zscore_cols(pool_mat)

    lo_p, hi_p = (1 - ACT_TOL) * load_prot, (1 + ACT_TOL) * load_prot
    lo_c, hi_c = (1 - ACT_TOL) * load_cyto, (1 + ACT_TOL) * load_cyto
    k_prot, k_cyto = len(prot_present), len(cyto_present)

    null_r = []
    attempts, max_attempts = 0, N_PERM * 200
    while len(null_r) < N_PERM and attempts < max_attempts:
        attempts += 1
        idx = RNG.choice(len(pool), size=k_prot + k_cyto, replace=False)
        pi, ci = idx[:k_prot], idx[k_prot:]
        xs = np.nanmean(pool_z[:, pi], axis=1)
        ys = np.nanmean(pool_z[:, ci], axis=1)
        if np.nanstd(xs) < 1e-8 or np.nanstd(ys) < 1e-8:
            continue
        # matched-null: require similar activation loading to the real modules
        lx = abs(stats.pearsonr(xs, a)[0])
        ly = abs(stats.pearsonr(ys, a)[0])
        if not (lo_p <= lx <= hi_p and lo_c <= ly <= hi_c):
            continue
        rp, _ = partial_r(xs, ys, a)
        null_r.append(rp)

    null_r = np.asarray(null_r)
    if len(null_r) < 100:
        print(f"WARNING: only {len(null_r)} matched null pairs after {attempts} attempts; "
              f"consider relaxing ACT_TOL (currently {ACT_TOL}).")

    # --- 5/6. Empirical p and verdict ---
    n_ge = int(np.sum(null_r >= r_obs))
    p_emp = (1 + n_ge) / (1 + len(null_r))
    print(f"\nNull (matched): n={len(null_r)}  mean={null_r.mean():.4f}  sd={null_r.std():.4f}  "
          f"95pct={np.percentile(null_r, 95):.4f}")
    print(f"Empirical p (observed >= matched-random upper tail): {p_emp:.4g}")

    print("\n" + "=" * 60)
    if p_emp < 0.05 and r_obs > null_r.mean():
        verdict = "MECHANISM-SPECIFIC"
        print("VERDICT: MECHANISM-SPECIFIC")
        print("-> protrusion~cytotoxicity partial coupling exceeds matched random NK")
        print("   module pairs of the same size and activation loading (p<0.05).")
        print("-> the 'effector arm recovers' claim is specific, not generic activation covariance.")
    else:
        verdict = "NOT_SPECIFIC"
        print("VERDICT: NOT SPECIFIC beyond matched activation-module covariance")
        print("-> observed coupling is within the matched-null distribution.")
        print("-> DOWNGRADE §3.2 to: 'consistent with, but not specific beyond, the covariance")
        print("   of any two similarly activation-loaded NK modules.'")

    # --- Save ---
    pd.DataFrame({
        "metric": ["raw_r", "partial_r_given_activation", "null_mean", "null_sd",
                   "null_p95", "empirical_p", "n_null", "n_cells", "verdict"],
        "value": [r_raw, r_obs, null_r.mean(), null_r.std(),
                  np.percentile(null_r, 95), p_emp, len(null_r), n_valid, verdict],
    }).to_csv(T + "t16_specificity_control.tsv", sep="\t", index=False)
    pd.DataFrame({"null_partial_r": null_r}).to_csv(
        T + "t16_null_distribution.tsv", sep="\t", index=False)
    print(f"\nSaved {T}t16_specificity_control.tsv and {T}t16_null_distribution.tsv")
    print("T16 PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
