"""
GC-NKGraph-Atlas Second Mechanism Card, End-to-End Recovery: TGFβ→SMAD→NK exclusion.

Runs the `tgfb_nk_exclusion` mechanism card's pre-registered hypotheses (H2-H5)
end-to-end on the card's positive-control cohort (TCGA-STAD) plus the two
external gastric cohorts, using the same module-scoring convention as the
Zheng card (`src/topology/sst_axis_validation.py`) and the same NK-fraction
purity control as the M5 task card
(`src/topology/bulk_h3_purity_control.py::partial_corr`).

Purpose: test whether the Zheng card's transcriptional-reach boundary
generalizes to a second, independently authored mechanism, not to show that a
transcriptional (SMAD) mechanism recovers more cleanly than a metabolite-level
one. It does not: the effector-arm coupling (H3) recovers, but the
mechanism-specific causal predictions (H2: TGFβ signaling suppresses activating
receptors; H4: CAF-ECM excludes cytotoxic NK) both fail in bulk across all
three cohorts, mirroring the Zheng card's upstream metabolic coupling failure.
Note that H3's magnitude is partly tautological -- the card's
"activating receptors" module (KLRK1/NCR1-3/CD226/...) is itself part of the
NK activation program and therefore co-regulated with the cytotoxicity
outcome -- so the informative cross-card signal is the *consistent failure*
of H2/H4, not the size of H3 (see manuscript §4.2, Limitations item 12).

Usage:
    python src/interpretation/run_tgfb_card_recovery.py
"""

import sys
import time
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.topology.bulk_h3_purity_control import mean_zscore, partial_corr  # noqa: E402


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


PROC = "data/processed/bulk"
OUT = "results/tables/mechanism_card_tgfb_recovery.tsv"

# cohort -> (expression matrix path, module-scoring style). TCGA-STAD is the
# card's own positive-control cohort (fillna, matching sst_axis_validation.py);
# the two external cohorts use skipna (matching run_geo_external_validation.py).
# For this card's gene sets, coverage is complete in both GSE cohorts so the
# style choice does not change any value -- it is kept consistent with
# bulk_h3_purity_control.py purely for convention.
COHORTS = {
    "TCGA-STAD": (f"{PROC}/tcga_stad_expression.tsv", "fillna"),
    "GSE62254":  (f"{PROC}/gse62254_expression.tsv",  "skipna"),
    "GSE84437":  (f"{PROC}/gse84437_expression.tsv",  "skipna"),
}

MODULES = {
    "nk_tgfb_signaling": ["TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4",
                           "SMAD7", "SKIL"],
    "nk_activating_receptors": ["KLRK1", "NCR1", "NCR2", "NCR3", "CD226", "SLAMF6",
                                 "SLAMF7", "CRTAM", "KLRF1"],
    "caf_ecm_program": ["COL1A1", "COL1A2", "COL3A1", "FN1", "FAP", "ACTA2",
                         "TGFBR1", "TGFBR2", "LOX", "LOXL2", "POSTN", "THBS2",
                         "SPARC", "VCAN"],
    "nk_cytolytic_machinery": ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "LCP2",
                                "LAT", "VAV1", "TLN1", "ITGAL", "ITGB2"],
}
ILC1_UP = ["IL7R", "KIT", "RORA", "CXCR6"]
ILC1_DOWN = ["EOMES", "TBX21"]
# clean NK-fraction proxy: absent from every module tested here
NK_LINEAGE_CLEAN = ["NCAM1", "FCGR3A", "TYROBP"]

# (hypothesis id, module A, module B, expected sign) per the card's
# `validation.positive_control.prereg_hypotheses`
HYPOTHESES = [
    ("H2", "nk_tgfb_signaling", "nk_activating_receptors", "-"),
    ("H3", "nk_activating_receptors", "nk_cytolytic_machinery", "+"),
    ("H4", "caf_ecm_program", "nk_cytolytic_machinery", "-"),
]


def main() -> None:
    rows = []
    for cohort, (path, style) in COHORTS.items():
        expr = pd.read_csv(path, sep="\t", index_col=0)
        scores = {name: mean_zscore(expr, genes, style)[0] for name, genes in MODULES.items()}
        ilc1_up, _ = mean_zscore(expr, ILC1_UP, style)
        ilc1_down, _ = mean_zscore(expr, ILC1_DOWN, style)
        scores["nk_ilc1_score"] = ilc1_up - ilc1_down
        lineage, n_lineage = mean_zscore(expr, NK_LINEAGE_CLEAN, style)

        n_passed = n_tested = 0
        for hid, a, b, expected in HYPOTHESES:
            r, p = stats.pearsonr(scores[a].values, scores[b].values)
            passed = (r > 0) == (expected == "+") and p < 0.05
            n_passed += passed
            n_tested += 1
            rows.append(dict(cohort=cohort, hypothesis=f"{hid} {a}~{b}",
                              expected=expected, r=round(r, 3), p=f"{p:.2e}",
                              passed=bool(passed)))

        # H5: nk_ilc1_score ~ cytolytic (expected -)
        r, p = stats.pearsonr(scores["nk_ilc1_score"].values,
                               scores["nk_cytolytic_machinery"].values)
        passed = r < 0 and p < 0.05
        n_passed += passed
        n_tested += 1
        rows.append(dict(cohort=cohort, hypothesis="H5 ilc1~cytolytic",
                          expected="-", r=round(r, 3), p=f"{p:.2e}", passed=bool(passed)))

        # Purity-controlled H3 (the card's effector-arm headline), same
        # NK-fraction confound check as the Zheng card (M5).
        r0, _ = stats.pearsonr(scores["nk_activating_receptors"].values,
                                scores["nk_cytolytic_machinery"].values)
        r_partial, p_partial, ci_low, ci_high, _ = partial_corr(
            scores["nk_activating_receptors"].values,
            scores["nk_cytolytic_machinery"].values,
            lineage.values,
        )
        rows.append(dict(cohort=cohort, hypothesis="H3 purity-controlled (lineage)",
                          expected="+", r=round(r_partial, 3), p=f"{p_partial:.2e}",
                          passed=bool(r_partial > 0 and p_partial < 0.05),
                          note=f"zero-order {r0:.3f} -> partial {r_partial:.3f}; "
                               f"lineage {n_lineage}/{len(NK_LINEAGE_CLEAN)}"))

        log(f"{cohort:10s}  H2-H5 passed {n_passed}/{n_tested}  | "
            f"H3 zero-order {r0:.3f} -> purity-controlled {r_partial:.3f} (p={p_partial:.1e})")

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}")


if __name__ == "__main__":
    main()
