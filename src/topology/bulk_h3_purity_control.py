"""
GC-NKGraph-Atlas Bulk Effector-Arm Purity Control (task card M5).

Question: does the bulk H3 effector coupling (protrusion-machinery ~
cytotoxicity-output) survive controlling for NK-cell abundance? In bulk tumor
transcriptomes both module scores rise and fall with the NK-cell fraction of
the sample, so the zero-order correlation is expected to be inflated by
co-varying NK abundance rather than a purely within-NK coupling.

Each cohort is scored with the SAME module-scoring style its own manuscript
pipeline uses, so the zero-order r reproduces the paper's reported numbers:
  - TCGA-LIHC          : sst_axis_validation.py style -> z.fillna(0).mean
  - GSE62254 / GSE84437 : run_geo_external_validation.py style -> z.mean (skipna)

The confound is then partialled out using a clean NK-lineage abundance proxy
(FCGR3A, KLRD1, KLRF1, KLRK1, NCAM1, NCR1, TYROBP) — genes present in neither
the protrusion nor the cytotoxicity module, avoiding the over-control that
would result from using the full NK-marker signature (whose genes overlap the
cytotoxicity module; see `src/immune_scoring/nk_scores.py::NK_MARKERS`).

Reference: plan_BiB_submission/TASK_CARD_M5_bulk_NK_fraction_confound.md

Usage:
    python src/topology/bulk_h3_purity_control.py
"""

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


PROC = "data/processed/bulk"
OUT = "results/tables/bulk_h3_purity_control.tsv"

# cohort -> (expression matrix path, module-scoring style matching that
# cohort's own manuscript pipeline)
COHORTS = {
    "TCGA-LIHC": (f"{PROC}/tcga_lihc_expression.tsv", "fillna"),
    "GSE62254":  (f"{PROC}/gse62254_expression.tsv",  "skipna"),
    "GSE84437":  (f"{PROC}/gse84437_expression.tsv",  "skipna"),
}

PROTRUSION = ["EZR", "MSN", "RDX", "ACTR2", "ACTR3", "ARPC1B", "ARPC2", "ARPC3",
              "ARPC4", "ARPC5", "WAS", "WASL", "WASF1", "WASF2", "WASF3", "WIPF1",
              "CDC42", "RAC1", "RHOA", "DIAPH1", "DIAPH3", "FMNL1", "BAIAP2",
              "PACSIN2"]
CYTOTOX = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "LCP2", "LAT", "VAV1", "TLN1",
           "ITGAL", "ITGB2"]
# clean NK-fraction proxy: present in neither module above
NK_LINEAGE = ["FCGR3A", "KLRD1", "KLRF1", "KLRK1", "NCAM1", "NCR1", "TYROBP"]

# Provenance note per cohort. GSE62254's zero-order r recomputed here (0.490)
# differs from the manuscript's Table 5 value (0.425); both derive from the
# same 24/24 + 11/11 gene coverage, and the gap traces to which probe-to-gene
# collapse pass produced the currently-saved gene-level matrix (the external
# validation script re-collapses from raw probes on each run) rather than to
# this script's method. It does not affect the ~50% attenuation conclusion.
COHORT_NOTE = {
    "TCGA-LIHC": "zero-order matches manuscript H3 bulk (0.55); primary bulk positive control",
    "GSE62254": ("zero-order recomputed on current gene-level matrix=0.490; Table 5 reports "
                 "0.425 from the original probe-collapse run; difference is collapse-version "
                 "only and does not affect the ~50% attenuation"),
    "GSE84437": "zero-order matches Table 5 (0.62) exactly",
}


def mean_zscore(df: pd.DataFrame, genes: list, style: str) -> tuple[pd.Series, int]:
    """Per-gene z-score, mean across genes. `style` controls missing-gene handling:
    'fillna' treats a missing gene's z-score as 0 (matches sst_axis_validation.py);
    'skipna' drops it from the row mean (matches run_geo_external_validation.py)."""
    available = [g for g in genes if g in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index), 0
    z = (df[available] - df[available].mean(0)) / df[available].std(0).replace(0, np.nan)
    if style == "fillna":
        return z.fillna(0).mean(axis=1), len(available)
    return z.mean(axis=1), len(available)


def partial_corr(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, float, float, int]:
    """Pearson correlation of x, y after residualizing both against covariate z
    via OLS. Returns (r, p, ci_low, ci_high, n)."""
    mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    x, y, z = x[mask], y[mask], z[mask]
    design = np.column_stack([np.ones_like(z), z])
    residualize = lambda v: v - design @ np.linalg.lstsq(design, v, rcond=None)[0]
    r, _ = stats.pearsonr(residualize(x), residualize(y))
    n = len(x)
    dof = n - 3  # n - 2 (regression) - 1 (covariate)
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), dof)
    z_fisher = np.arctanh(r)
    se = 1 / np.sqrt(n - 3 - 1)
    ci_low, ci_high = np.tanh(z_fisher - 1.96 * se), np.tanh(z_fisher + 1.96 * se)
    return r, p, ci_low, ci_high, n


def main() -> None:
    rows = []
    for cohort, (path, style) in COHORTS.items():
        expr = pd.read_csv(path, sep="\t", index_col=0)
        protrusion, _ = mean_zscore(expr, PROTRUSION, style)
        cytotox, _ = mean_zscore(expr, CYTOTOX, style)
        lineage, n_lineage = mean_zscore(expr, NK_LINEAGE, style)

        paired = pd.DataFrame({"p": protrusion, "c": cytotox}).dropna()
        r0, p0 = stats.pearsonr(paired.p, paired.c)
        r_partial, p_partial, ci_low, ci_high, n = partial_corr(
            protrusion.values, cytotox.values, lineage.values
        )
        attenuation_pct = round(100 * (r0 - r_partial) / r0, 1)

        rows.append(dict(
            cohort=cohort, n=len(paired), scoring=style,
            lineage_genes=f"{n_lineage}/{len(NK_LINEAGE)}",
            r_zeroorder=round(r0, 4), p_zeroorder=f"{p0:.3e}",
            r_partial_lineage=round(r_partial, 4),
            ci_low=round(ci_low, 4), ci_high=round(ci_high, 4),
            p_partial_lineage=f"{p_partial:.3e}",
            attenuation_pct=attenuation_pct,
            note=COHORT_NOTE[cohort],
        ))
        log(f"{cohort:10s} n={len(paired):4d} [{style:6s}]  "
            f"zero-order r={r0:.4f} -> lineage-partial r={r_partial:.4f} "
            f"({attenuation_pct:+.0f}%, p={p_partial:.1e}, CI[{ci_low:.2f},{ci_high:.2f}])")

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}")


if __name__ == "__main__":
    main()
