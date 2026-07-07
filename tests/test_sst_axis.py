"""
Tests for src/topology/sst_axis.py — SST-axis score computation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.topology.sst_axis import (
    mean_zscore,
    compute_sst_scores,
    compute_sample_summary,
)
from src.common.sst_config import load_sst_modules


# ---------------------------------------------------------------------------
# mean_zscore
# ---------------------------------------------------------------------------

class TestMeanZscore:
    """Test the core scoring helper."""

    def test_returns_correct_shape(self, sample_expression_df):
        genes = ["PHGDH", "PSAT1", "PSPH"]
        result = mean_zscore(sample_expression_df, genes)
        assert len(result) == len(sample_expression_df)

    def test_missing_genes_returns_zeros(self, sample_expression_df):
        result = mean_zscore(sample_expression_df, ["NOT_A_GENE"])
        assert (result == 0.0).all()

    def test_all_missing_returns_zeros(self, sample_expression_df):
        result = mean_zscore(sample_expression_df, [])
        assert (result == 0.0).all()

    def test_single_gene_matches_its_zscore(self):
        """With one gene, result should equal that gene's z-score."""
        df = pd.DataFrame({
            "A": [1.0, 2.0, 3.0, 4.0, 5.0],
        })
        result = mean_zscore(df, ["A"])
        expected = (df["A"] - df["A"].mean()) / df["A"].std(ddof=0)
        np.testing.assert_array_almost_equal(result.values, expected.values)

    def test_mean_of_two_genes(self):
        """Mean z-score of two genes = average of individual z-scores."""
        rng = np.random.RandomState(123)
        df = pd.DataFrame(rng.randn(50, 2), columns=["A", "B"])
        result = mean_zscore(df, ["A", "B"])
        z_a = (df["A"] - df["A"].mean()) / df["A"].std(ddof=0)
        z_b = (df["B"] - df["B"].mean()) / df["B"].std(ddof=0)
        expected = (z_a + z_b) / 2
        np.testing.assert_array_almost_equal(result.values, expected.values)


# ---------------------------------------------------------------------------
# compute_sst_scores
# ---------------------------------------------------------------------------

class TestComputeSstScores:
    """Integration test for the full SST scoring pipeline."""

    @pytest.fixture
    def adata(self, sample_expression_df):
        """Create a minimal AnnData from the sample expression DataFrame."""
        try:
            import anndata
            adata = anndata.AnnData(
                X=sample_expression_df.values,
                obs=pd.DataFrame(index=sample_expression_df.index),
                var=pd.DataFrame(index=sample_expression_df.columns),
            )
            adata.obs["sample_id"] = [f"S{i % 4}" for i in range(len(adata))]
            adata.obs["tissue"] = "tumor"
            adata.obs["condition"] = [
                "tumor" if i < 70 else "normal" for i in range(len(adata))
            ]
            return adata
        except ImportError:
            pytest.skip("anndata not installed")

    def test_output_is_dataframe(self, adata, sample_sst_modules):
        scores = compute_sst_scores(adata, modules=sample_sst_modules)
        assert isinstance(scores, pd.DataFrame)

    def test_output_has_expected_columns(self, adata, sample_sst_modules):
        scores = compute_sst_scores(adata, modules=sample_sst_modules)
        assert "nk_sm_balance_score" in scores.columns
        assert "nk_topology_permissive_score" in scores.columns
        assert "sst_axis_score" in scores.columns
        assert "sample_id" in scores.columns
        assert "tissue" in scores.columns
        assert "condition" in scores.columns

    def test_output_same_index_as_input(self, adata, sample_sst_modules):
        scores = compute_sst_scores(adata, modules=sample_sst_modules)
        assert scores.index.equals(adata.obs_names)

    def test_balance_is_synthesis_minus_catabolism(self, adata, sample_sst_modules):
        scores = compute_sst_scores(adata, modules=sample_sst_modules)
        diff = (
            scores["nk_sm_synthesis_score"]
            - scores["nk_sm_catabolism_score"]
        )
        np.testing.assert_array_almost_equal(
            scores["nk_sm_balance_score"].values, diff.values
        )

    def test_no_nan_in_scores(self, adata, sample_sst_modules):
        scores = compute_sst_scores(adata, modules=sample_sst_modules)
        score_cols = [c for c in scores.columns if c.endswith("_score")]
        for c in score_cols:
            assert not scores[c].isna().any(), f"NaN in {c}"


# ---------------------------------------------------------------------------
# compute_sample_summary
# ---------------------------------------------------------------------------

class TestComputeSampleSummary:
    """Test per-sample aggregation."""

    @pytest.fixture
    def scores_df(self):
        rng = np.random.RandomState(42)
        return pd.DataFrame({
            "tumor_serine_capacity_score": rng.randn(200),
            "nk_sm_balance_score": rng.randn(200),
            "nk_topology_permissive_score": rng.randn(200),
            "sst_axis_score": rng.randn(200),
            "sample_id": [f"S{i % 5}" for i in range(200)],
            "tissue": "tumor",
            "condition": ["tumor"] * 200,
        })

    def test_returns_dataframe(self, scores_df):
        summary = compute_sample_summary(scores_df)
        assert isinstance(summary, pd.DataFrame)

    def test_groups_by_sample(self, scores_df):
        summary = compute_sample_summary(scores_df)
        assert summary.index.nlevels == 3 or len(summary) == 5

    def test_handles_missing_group_cols(self, scores_df):
        df = scores_df.drop(columns=["tissue"])
        summary = compute_sample_summary(df)
        assert isinstance(summary, pd.DataFrame)
