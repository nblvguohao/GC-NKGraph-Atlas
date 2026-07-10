"""
T16 robustness — stress-test the "effector-arm coupling is NOT mechanism-specific"
verdict before any manuscript wording is changed.

Two worries about the matched-null permutation in run_t16_specificity_control.py:
  (A) The null may be "too strong": random NK module averages co-vary through
      global covariates (sequencing depth, cell size / total transcriptional
      activity) that a single 16-gene activation signature does not remove. If so,
      additionally controlling for depth might make the real coupling stand out.
  (B) The verdict may hinge on the activation-loading match tolerance (ACT_TOL).

This script recomputes the observed statistic and the matched null under:
  - control set = {activation}          vs {activation + log1p(total_counts)}
  - matching band ACT_TOL in {0.25, 0.5, 1.0} and an UNMATCHED null
and reports empirical p for each cell. Verdict is robust only if observed stays
inside the null across variants.

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_t16_robustness.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

warnings.filterwarnings("ignore")

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.common.sst_config import load_sst_modules, get_sst_genes  # noqa: E402

T = "results/tables/"
H5AD = "data/processed/scrna/gc_nk_subset.h5ad"
N_PERM = 1000
RNG = np.random.default_rng(42)
MIN_DETECT_FRAC = 0.05

NK_ACTIVATION_GENES = [
    "CD69", "TNF", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5", "CSF2",
    "IL2RA", "ICOS", "TNFSF10", "FASLG", "CD38", "HLA-DRA", "HLA-DRB1", "MKI67",
]


def zscore_cols(mat):
    return (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-10)


def residualize(v, C):
    """Residual of v after OLS on design C (intercept added)."""
    A = np.column_stack([np.ones(len(v)), C])
    beta, *_ = np.linalg.lstsq(A, v, rcond=None)
    return v - A @ beta


def partial_r_multi(x, y, C):
    xr, yr = residualize(x, C), residualize(y, C)
    return stats.pearsonr(xr, yr)[0]


def main() -> int:
    import anndata as ad
    modules = load_sst_modules()
    prot = [g for g in modules["nk_protrusion_machinery"]["genes"]]
    cyto = [g for g in modules["nk_synapse_cytotoxicity_outcome"]["genes"]]
    sst_all = set(get_sst_genes(modules))

    adata = ad.read_h5ad(H5AD)
    if "cell_type" in adata.obs.columns:
        adata = adata[adata.obs["cell_type"] == "NK"].copy()
    var = list(map(str, adata.var_names))
    vi = {g: i for i, g in enumerate(var)}
    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)

    prot_p = [g for g in prot if g in vi]
    cyto_p = [g for g in cyto if g in vi]
    act_p = [g for g in NK_ACTIVATION_GENES if g in vi]

    def dense(genes):
        sub = X[:, [vi[g] for g in genes]]
        return sub.toarray() if sp.issparse(sub) else np.asarray(sub)

    x_obs = zscore_cols(dense(prot_p)).mean(axis=1)
    y_obs = zscore_cols(dense(cyto_p)).mean(axis=1)
    activation = zscore_cols(dense(act_p)).mean(axis=1)

    # depth covariate
    if "total_counts" in adata.obs.columns:
        depth = np.log1p(adata.obs["total_counts"].to_numpy(dtype=float))
    else:
        depth = np.log1p(np.asarray(X.sum(axis=1)).ravel())
    depth = (depth - depth.mean()) / (depth.std() + 1e-10)

    valid = ~(np.isnan(x_obs) | np.isnan(y_obs) | np.isnan(activation) | np.isnan(depth))
    x_obs, y_obs, act, dep = x_obs[valid], y_obs[valid], activation[valid], depth[valid]
    load_prot = abs(stats.pearsonr(x_obs, act)[0])
    load_cyto = abs(stats.pearsonr(y_obs, act)[0])
    print(f"cells={valid.sum()}  activation loading prot={load_prot:.3f} cyto={load_cyto:.3f}")

    # null pool
    exclude = sst_all | set(NK_ACTIVATION_GENES)
    detect = np.asarray((X > 0).mean(axis=0)).ravel()
    pool = [i for i, g in enumerate(var) if g not in exclude and detect[i] >= MIN_DETECT_FRAC]
    pmat = X[:, pool]
    pmat = (pmat.toarray() if sp.issparse(pmat) else np.asarray(pmat))[valid]
    pool_z = zscore_cols(pmat)
    kP, kC = len(prot_p), len(cyto_p)

    controls = {"activation": np.column_stack([act]),
                "activation+depth": np.column_stack([act, dep])}
    tols = {"match0.25": 0.25, "match0.5": 0.5, "match1.0": 1.0, "unmatched": None}

    rows = []
    for cname, C in controls.items():
        r_obs = partial_r_multi(x_obs, y_obs, C)
        for tname, tol in tols.items():
            if tol is None:
                lo_p = hi_p = lo_c = hi_c = None
            else:
                lo_p, hi_p = (1 - tol) * load_prot, (1 + tol) * load_prot
                lo_c, hi_c = (1 - tol) * load_cyto, (1 + tol) * load_cyto
            null, attempts, cap = [], 0, N_PERM * 300
            while len(null) < N_PERM and attempts < cap:
                attempts += 1
                idx = RNG.choice(len(pool), size=kP + kC, replace=False)
                xs = pool_z[:, idx[:kP]].mean(axis=1)
                ys = pool_z[:, idx[kP:]].mean(axis=1)
                if xs.std() < 1e-8 or ys.std() < 1e-8:
                    continue
                if tol is not None:
                    lx, ly = abs(stats.pearsonr(xs, act)[0]), abs(stats.pearsonr(ys, act)[0])
                    if not (lo_p <= lx <= hi_p and lo_c <= ly <= hi_c):
                        continue
                null.append(partial_r_multi(xs, ys, C))
            null = np.asarray(null)
            p_emp = (1 + int((null >= r_obs).sum())) / (1 + len(null))
            verdict = "SPECIFIC" if (p_emp < 0.05 and r_obs > null.mean()) else "not_spec"
            rows.append(dict(control=cname, null=tname, r_obs=round(r_obs, 4),
                             null_mean=round(float(null.mean()), 4),
                             null_p95=round(float(np.percentile(null, 95)), 4),
                             n_null=len(null), emp_p=round(p_emp, 4), verdict=verdict))
            print(f"  [{cname:17s} | {tname:9s}] r_obs={r_obs:.3f} "
                  f"null_mean={null.mean():.3f} p95={np.percentile(null,95):.3f} "
                  f"n={len(null):4d} emp_p={p_emp:.3f} -> {verdict}")

    out = pd.DataFrame(rows)
    out.to_csv(T + "t16_robustness.tsv", sep="\t", index=False)
    print(f"\nSaved {T}t16_robustness.tsv")
    n_spec = (out["verdict"] == "SPECIFIC").sum()
    print(f"\nSPECIFIC in {n_spec}/{len(out)} variants.")
    if n_spec == 0:
        print("ROBUST: NOT_SPECIFIC holds across all control sets and matching bands.")
    elif n_spec == len(out):
        print("REVERSED: SPECIFIC across all variants -> original null was too strong.")
    else:
        print("MIXED: verdict depends on control/matching -> report sensitivity explicitly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
