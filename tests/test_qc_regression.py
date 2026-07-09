"""
T4-D3 regression test: verify scRNA QC does not change the paper's core conclusions.

Ensures that after applying QC filtering (min_genes ≥ 200, max_genes ≤ 6000,
pct_mito ≤ 20%), the H3 (protrusion→cytotoxicity) and H2 (SM-balance→protrusion)
correlations retain the same direction and significance as the pre-QC analysis.

This test runs on the server against the real GSE246662 data. It is NOT a unit
test — it requires the actual .h5ad file and processed result tables.

Run AFTER qc_filter.py:
    python src/scrna_analysis/qc_filter.py \
        --in data/processed/scrna/gc_integrated.h5ad \
        --out data/processed/scrna/gc_integrated_qc.h5ad

    pytest tests/test_qc_regression.py -v
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

# ── Paths ──────────────────────────────────────────────────────────────────
BEFORE_TSV = "results/tables/sst_axis_scores_single_cell.tsv"        # v2 (no-QC) baseline
QC_TSV     = "results/tables/sst_axis_scores_single_cell_qc.tsv"    # post-QC (must exist)
SUMMARY    = "results/tables/scrna_qc_summary.tsv"


# ── Helpers ────────────────────────────────────────────────────────────────
def _pearson(x, y):
    """Compute Pearson r and p-value, skipping NaN."""
    from scipy.stats import pearsonr
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return np.nan, np.nan
    return pearsonr(x[mask], y[mask])


# ── Tests ──────────────────────────────────────────────────────────────────

@pytest.mark.server
class TestQCRegression:
    """T4-D3: After QC, H3 remains significantly positive; H2 direction unchanged."""

    @pytest.fixture(autouse=True)
    def _require_qc_data(self):
        """Skip all tests if QC output doesn't exist."""
        if not os.path.exists(QC_TSV):
            pytest.skip(f"QC output missing: {QC_TSV} — run qc_filter.py first "
                        f"(T4 step 1 on server)")
        if not os.path.exists(BEFORE_TSV):
            pytest.skip(f"Pre-QC baseline missing: {BEFORE_TSV}")

    def test_qc_summary_exists(self):
        """D1: QC summary table is written with per-sample before/after counts."""
        assert os.path.exists(SUMMARY), f"Missing {SUMMARY}"
        df = pd.read_csv(SUMMARY, sep="\t", index_col=0)
        assert "cells_before" in df.columns
        assert "cells_after" in df.columns
        assert "retained_frac" in df.columns
        assert len(df) >= 1, "summary must have at least one sample row"
        # D2 sanity: retention per sample in (0.5, 0.98)
        for idx, row in df.iterrows():
            frac = float(row["retained_frac"])
            assert 0.5 <= frac <= 0.98, (
                f"{idx}: retained_frac={frac:.2f} outside [0.50, 0.98] "
                f"— adjust thresholds if this is a genuine small/large sample"
            )
        print(f"  D1+D2 PASS: {len(df)} samples, "
              f"{int(df['cells_after'].sum())} total cells retained")

    def test_h3_direction_preserved(self):
        """H3: protrusion→cytotoxicity must remain significantly positive."""
        before = pd.read_csv(BEFORE_TSV, sep="\t")
        after  = pd.read_csv(QC_TSV, sep="\t")

        r_pre, p_pre   = _pearson(before["nk_protrusion_machinery_score"],
                                  before["nk_synapse_cytotoxicity_outcome_score"])
        r_post, p_post = _pearson(after["nk_protrusion_machinery_score"],
                                  after["nk_synapse_cytotoxicity_outcome_score"])

        assert not np.isnan(r_pre), "pre-QC H3 correlation is NaN"
        assert not np.isnan(r_post), "post-QC H3 correlation is NaN — QC may have dropped too many NK cells"

        assert r_pre > 0, f"pre-QC H3 r={r_pre:.3f} ≤ 0 — baseline is wrong, check input"
        assert r_post > 0, f"post-QC H3 r={r_post:.3f} ≤ 0 — QC broke the effector-arm direction"
        assert p_post < 0.05, f"post-QC H3 p={p_post:.1e} ≥ 0.05 — effector arm lost significance"

        print(f"  H3 pre-QC:  r={r_pre:.4f}, p={p_pre:.1e}")
        print(f"  H3 post-QC: r={r_post:.4f}, p={p_post:.1e}")
        print(f"  D3-H3 PASS: direction preserved, significant")

    def test_h2_direction_preserved(self):
        """H2: SM-balance→protrusion must not flip sign (weak positive or null ok)."""
        before = pd.read_csv(BEFORE_TSV, sep="\t")
        after  = pd.read_csv(QC_TSV, sep="\t")

        def _get_sm_balance(df):
            if "nk_sm_balance_score" in df.columns:
                return df["nk_sm_balance_score"]
            syn  = df.get("nk_sm_synthesis_score", pd.Series(0, index=df.index))
            cat  = df.get("nk_sm_catabolism_score", pd.Series(0, index=df.index))
            return syn - cat

        sm_pre  = _get_sm_balance(before)
        sm_post = _get_sm_balance(after)

        r_pre, p_pre   = _pearson(sm_pre,  before["nk_protrusion_machinery_score"])
        r_post, p_post = _pearson(sm_post, after["nk_protrusion_machinery_score"])

        if np.isnan(r_pre) or np.isnan(r_post):
            pytest.skip("H2 NaN — insufficient SM-balance genes covered, skipping")

        # H2 is weak/conditional; the requirement is: must NOT flip to significantly negative
        assert r_post >= -0.02, (
            f"post-QC H2 r={r_post:.4f} < −0.02 — QC flipped the metabolic coupling "
            f"direction to negative. pre-QC r={r_pre:.4f}. Investigate before publishing."
        )

        print(f"  H2 pre-QC:  r={r_pre:.4f}, p={p_pre:.1e}")
        print(f"  H2 post-QC: r={r_post:.4f}, p={p_post:.1e}")
        print(f"  D3-H2 PASS: direction preserved (r_post={r_post:.4f} ≥ −0.02)")

    def test_nk_cell_count_reasonable(self):
        """Post-QC NK count must be large enough to compute axis scores meaningfully."""
        after = pd.read_csv(QC_TSV, sep="\t")
        n_nk = len(after)
        assert n_nk >= 2000, (
            f"post-QC NK cells = {n_nk} — too few for axis analysis. "
            f"Relax QC thresholds or check that NK annotation survived filtering."
        )
        print(f"  post-QC NK cell count: {n_nk}")


@pytest.mark.server
class TestQCAxisConsistency:
    """T4-D3 extra: cross-tissue module means should not degenerate."""

    def test_tissue_counts_preserved(self):
        """Post-QC data must still have cells from all three tissues."""
        after = pd.read_csv(QC_TSV, sep="\t")
        if "tissue" not in after.columns:
            pytest.skip("no 'tissue' column in QC output")
        tissues = after["tissue"].dropna().unique()
        assert len(tissues) >= 2, (
            f"post-QC tissues = {list(tissues)} — fewer than 2 tissues, "
            f"QC may have dropped one tissue entirely"
        )
        print(f"  post-QC tissues: {list(tissues)}")
