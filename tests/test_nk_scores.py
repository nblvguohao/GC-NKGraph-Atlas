"""
Tests for src/immune_scoring/nk_scores.py — NK immune state scoring.
"""

from __future__ import annotations

import tempfile
import os
import numpy as np
import pandas as pd
import pytest

from src.immune_scoring.nk_scores import (
    mean_zscore,
    compute_nk_scores,
    assign_immune_states,
    NK_MARKERS,
    NK_CYTOTOXICITY_GENES,
    NK_DYSFUNCTION_GENES,
    CAF_ECM_TGFB_GENES,
)
from src.common.log_utils import Logger


@pytest.fixture
def logger():
    """Create a temporary Logger that writes to a temp file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        log_path = os.path.join(tmpdir, "test_log.md")
        yield Logger(log_path=log_path)


class TestMeanZscore:
    """Tests for mean_zscore — the core scoring primitive."""

    def test_returns_series_with_correct_index(self):
        df = pd.DataFrame({
            "GENE1": [1.0, 2.0, 3.0],
            "GENE2": [4.0, 5.0, 6.0],
        }, index=["A", "B", "C"])
        result = mean_zscore(df, ["GENE1", "GENE2"])
        assert isinstance(result, pd.Series)
        assert list(result.index) == ["A", "B", "C"]

    def test_all_genes_missing_returns_zeros(self):
        df = pd.DataFrame({"GENE1": [1.0, 2.0]}, index=["A", "B"])
        result = mean_zscore(df, ["NONEXISTENT"])
        assert (result == 0.0).all()

    def test_single_gene_equals_its_zscore(self):
        df = pd.DataFrame({"GENE1": [0.0, 1.0, 2.0]}, index=["A", "B", "C"])
        result = mean_zscore(df, ["GENE1"])
        vals = np.array([0.0, 1.0, 2.0])
        z = (vals - vals.mean()) / vals.std(ddof=0)
        np.testing.assert_array_almost_equal(result.values, z)

    def test_mean_of_two_genes(self):
        df = pd.DataFrame({
            "G1": [0.0, 2.0],
            "G2": [2.0, 0.0],
        }, index=["A", "B"])
        result = mean_zscore(df, ["G1", "G2"])
        np.testing.assert_array_almost_equal(result.values, [0.0, 0.0])

    def test_partial_gene_coverage(self):
        df = pd.DataFrame({"G1": [0.0, 1.0]}, index=["A", "B"])
        result = mean_zscore(df, ["G1", "MISSING"])
        vals = np.array([0.0, 1.0])
        z = (vals - vals.mean()) / vals.std(ddof=0)
        np.testing.assert_array_almost_equal(result.values, z)

    def test_no_nan_in_output(self):
        df = pd.DataFrame({"GENE1": [1.0, 2.0, 3.0]}, index=["A", "B", "C"])
        result = mean_zscore(df, ["GENE1"])
        assert not result.isna().any()

    def test_nan_in_input_handled(self):
        df = pd.DataFrame({"G1": [1.0, np.nan, 3.0]}, index=["A", "B", "C"])
        result = mean_zscore(df, ["G1"])
        assert not result.isna().any()


class TestComputeNkScores:
    """Integration test for NK score computation."""

    @pytest.fixture
    def expr_df(self):
        rng = np.random.RandomState(42)
        genes = NK_MARKERS + NK_DYSFUNCTION_GENES + CAF_ECM_TGFB_GENES
        data = rng.lognormal(mean=2.0, sigma=0.5, size=(50, len(genes)))
        df = pd.DataFrame(data, columns=genes)
        df.index = [f"SAMPLE_{i:03d}" for i in range(50)]
        df.index.name = "sample_id"
        return df

    def test_returns_dataframe(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        assert isinstance(result, pd.DataFrame)

    def test_has_expected_columns(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        expected = [
            "NK_infiltration_score",
            "NK_cytotoxicity_score",
            "NK_dysfunction_score",
            "NK_exclusion_score",
        ]
        for col in expected:
            assert col in result.columns

    def test_same_index_as_input(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        assert list(result.index) == list(expr_df.index)

    def test_no_nan_in_scores(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        assert not result.isna().any().any()

    def test_dysfunction_formula(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        expected = mean_zscore(expr_df, NK_DYSFUNCTION_GENES) - mean_zscore(expr_df, NK_CYTOTOXICITY_GENES)
        assert (result["NK_dysfunction_score"] - expected).abs().max() < 1e-10

    def test_exclusion_formula(self, expr_df, logger):
        result = compute_nk_scores(expr_df, logger, "test_dataset")
        expected = mean_zscore(expr_df, CAF_ECM_TGFB_GENES) - result["NK_infiltration_score"]
        assert (result["NK_exclusion_score"] - expected).abs().max() < 1e-10


class TestAssignImmuneStates:
    """Tests for immune state label assignment. Returns (DataFrame, thresholds_dict)."""

    @pytest.fixture
    def scores_df(self):
        return pd.DataFrame({
            "NK_infiltration_score": [1.0, 1.0, -1.0, 0.0],
            "NK_cytotoxicity_score": [1.0, -1.0, 0.0, 0.0],
            "NK_dysfunction_score": [-1.0, 1.0, 0.0, 0.0],
            "NK_exclusion_score": [-1.0, -1.0, 1.0, 0.0],
        }, index=["S1", "S2", "S3", "S4"])

    def test_returns_tuple(self, scores_df, logger):
        result = assign_immune_states(scores_df, logger=logger)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_first_element_is_dataframe(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        assert isinstance(result, pd.DataFrame)

    def test_second_element_is_threshold_dict(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        assert isinstance(thresholds, dict)
        assert "infiltration_high" in thresholds

    def test_has_immune_state_column(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        assert "nk_immune_state" in result.columns

    def test_all_labels_in_valid_set(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        valid = {
            "NK-hot-cytotoxic",
            "NK-hot-dysfunctional",
            "NK-cold/excluded",
            "NK-intermediate",
        }
        for label in result["nk_immune_state"]:
            assert label in valid

    def test_returns_four_labels(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        assert len(result) == 4

    def test_custom_thresholds(self, scores_df, logger):
        thresholds = {
            "infiltration_high": 0.5,
            "cytotoxicity_high": 0.5,
            "dysfunction_high": 0.5,
            "exclusion_high": 0.5,
        }
        result, thresh = assign_immune_states(scores_df, thresholds, logger)
        assert isinstance(result, pd.DataFrame)

    def test_output_preserves_index(self, scores_df, logger):
        result, thresholds = assign_immune_states(scores_df, logger=logger)
        assert list(result.index) == ["S1", "S2", "S3", "S4"]
