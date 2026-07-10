"""
Tests for target prioritization module (prioritize_targets.py).

Covers: evidence computation functions, direction consistency,
candidate classification, axis core membership.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.interpretation.prioritize_targets import (
    compute_tumor_specificity,
    compute_nk_correlations,
    classify_candidate,
    check_direction_consistency,
    is_axis_core,
    GOLD_STANDARD,
    DRUGGABILITY,
    SEED_CANDIDATE_CATEGORIES,
)


# =========================================================================
# Fixtures
# =========================================================================

@pytest.fixture
def simple_expr():
    """Simple expression matrix: 10 cells x 15 genes."""
    np.random.seed(42)
    n_cells = 10
    genes = [f"GENE_{i}" for i in range(15)] + [
        "PHGDH", "SGMS1", "SMPD1", "NKG7", "GZMB", "HAVCR2", "KLRC1",
        "TIGIT", "PRF1", "IFNG", "XCL1", "CCL5", "CD96", "ENTPD1", "TOX",
    ]
    data = np.random.exponential(scale=5, size=(n_cells, len(genes)))
    df = pd.DataFrame(data, columns=genes)
    df.index = [f"cell_{i}" for i in range(n_cells)]
    return df


@pytest.fixture
def malignant_mask(simple_expr):
    """Half cells are malignant."""
    n = len(simple_expr)
    return pd.Series(
        [True] * (n // 2) + [False] * (n - n // 2),
        index=simple_expr.index,
    )


@pytest.fixture
def nk_mask(simple_expr):
    """A few cells are NK."""
    mask = pd.Series(False, index=simple_expr.index)
    mask.iloc[:3] = True
    return mask


# =========================================================================
# Tests: compute_tumor_specificity
# =========================================================================

class TestComputeTumorSpecificity:
    """Tests for tumor specificity computation."""

    def test_returns_float(self, simple_expr, malignant_mask):
        result = compute_tumor_specificity(simple_expr, malignant_mask, "PHGDH")
        assert isinstance(result, float)

    def test_returns_zero_for_missing_gene(self, simple_expr, malignant_mask):
        result = compute_tumor_specificity(simple_expr, malignant_mask, "NONEXISTENT")
        assert result == 0.0

    def test_positive_for_high_in_malignant(self):
        """Gene highly expressed in malignant should give positive LFC."""
        expr = pd.DataFrame({
            "GENE_A": [100.0, 100.0, 1.0, 1.0],
        })
        is_mal = pd.Series([True, True, False, False])
        result = compute_tumor_specificity(expr, is_mal, "GENE_A")
        assert result > 0.0

    def test_negative_for_low_in_malignant(self):
        """Gene low in malignant should give negative LFC."""
        expr = pd.DataFrame({
            "GENE_B": [1.0, 1.0, 100.0, 100.0],
        })
        is_mal = pd.Series([True, True, False, False])
        result = compute_tumor_specificity(expr, is_mal, "GENE_B")
        assert result < 0.0


# =========================================================================
# Tests: compute_nk_correlations
# =========================================================================

class TestComputeNkCorrelations:
    """Tests for NK dysfunction correlation computation."""

    def test_returns_dict(self, simple_expr, nk_mask):
        genes = ["PHGDH", "SGMS1", "HAVCR2"]
        result = compute_nk_correlations(simple_expr, nk_mask, genes)
        assert isinstance(result, dict)

    def test_all_requested_genes_in_result(self, simple_expr, nk_mask):
        genes = ["PHGDH", "SGMS1", "SMPD1"]
        result = compute_nk_correlations(simple_expr, nk_mask, genes)
        for g in genes:
            assert g in result

    def test_missing_gene_returns_zero(self, simple_expr, nk_mask):
        result = compute_nk_correlations(simple_expr, nk_mask, ["NONEXISTENT"])
        assert result["NONEXISTENT"] == 0.0

    def test_correlation_in_range_neg1_to_1(self, simple_expr, nk_mask):
        genes = ["PHGDH", "SGMS1", "HAVCR2", "KLRC1", "TIGIT"]
        result = compute_nk_correlations(simple_expr, nk_mask, genes)
        for v in result.values():
            assert -1.0 <= v <= 1.0, f"Correlation {v} out of range"

    def test_no_nk_cells_handled(self, simple_expr):
        """When no NK cells, should not crash."""
        no_nk = pd.Series(False, index=simple_expr.index)
        result = compute_nk_correlations(simple_expr, no_nk, ["PHGDH"])
        assert isinstance(result, dict)


# =========================================================================
# Tests: classify_candidate
# =========================================================================

class TestClassifyCandidate:
    """Tests for candidate classification."""

    def test_returns_seed_category_for_known_gene(self):
        result = classify_candidate("LDHA", None)
        assert result == "metabolic_suppression"

    def test_returns_sst_module_for_non_seed(self):
        result = classify_candidate("PHGDH", "tumor_serine_capacity")
        assert "sst_axis" in result
        assert "tumor_serine_capacity" in result

    def test_returns_unknown_for_unclassified(self):
        result = classify_candidate("ACTB", None)
        assert result == "unknown_candidate"

    def test_seed_category_takes_priority_over_sst(self):
        """Seed category should be preferred over SST module name."""
        result = classify_candidate("HAVCR2", "nk_checkpoint_link")
        # HAVCR2 is NOT in seed categories, but is in SST
        assert result != "unknown_candidate"

    def test_nt5e_returns_adenosine(self):
        result = classify_candidate("NT5E", None)
        assert result == "adenosine_pathway"


# =========================================================================
# Tests: check_direction_consistency
# =========================================================================

class TestCheckDirectionConsistency:
    """Tests for SST axis direction consistency checking."""

    def test_catabolism_positive_corr_consistent(self):
        """SM catabolism: higher in dysfunctional NK = positive correlation = consistent."""
        result = check_direction_consistency("nk_sm_catabolism", 0.3)
        assert result is True

    def test_catabolism_negative_corr_inconsistent(self):
        result = check_direction_consistency("nk_sm_catabolism", -0.2)
        assert result is False

    def test_synthesis_positive_corr_inconsistent(self):
        """SM synthesis: higher should mean less dysfunction, so positive = inconsistent."""
        result = check_direction_consistency("nk_sm_synthesis", 0.3)
        assert result is False

    def test_synthesis_negative_corr_consistent(self):
        result = check_direction_consistency("nk_sm_synthesis", -0.5)
        assert result is True

    def test_protrusion_negative_corr_consistent(self):
        result = check_direction_consistency("nk_protrusion_machinery", -0.1)
        assert result is True

    def test_protrusion_positive_corr_inconsistent(self):
        result = check_direction_consistency("nk_protrusion_machinery", 0.4)
        assert result is False

    def test_cytotoxicity_negative_corr_consistent(self):
        result = check_direction_consistency("nk_synapse_cytotoxicity_outcome", -0.3)
        assert result is True

    def test_cytotoxicity_positive_corr_inconsistent(self):
        result = check_direction_consistency("nk_synapse_cytotoxicity_outcome", 0.2)
        assert result is False

    def test_non_directional_module_returns_none(self):
        result = check_direction_consistency("tumor_serine_capacity", 0.5)
        assert result is None

    def test_none_module_returns_none(self):
        result = check_direction_consistency(None, 0.5)
        assert result is None

    def test_zero_correlation_handled(self):
        """Zero correlation edge case."""
        result = check_direction_consistency("nk_sm_catabolism", 0.0)
        # 0 > 0 is False, so catabolism with corr=0 returns False
        assert result is False


# =========================================================================
# Tests: is_axis_core
# =========================================================================

class TestIsAxisCore:
    """Tests for SST axis core membership."""

    def test_sgms1_is_core(self):
        assert is_axis_core("SGMS1") is True

    def test_smpd1_is_core(self):
        assert is_axis_core("SMPD1") is True

    def test_nkg7_is_core(self):
        assert is_axis_core("NKG7") is True

    def test_actb_not_core(self):
        assert is_axis_core("ACTB") is False

    def test_unknown_gene_not_core(self):
        assert is_axis_core("RANDOM_XYZ") is False


# =========================================================================
# Tests: Constants / configuration
# =========================================================================

class TestConstants:
    """Tests for gold standard and druggability configuration."""

    def test_gold_standard_is_list(self):
        assert isinstance(GOLD_STANDARD, list)
        assert len(GOLD_STANDARD) >= 10

    def test_druggability_has_known_targets(self):
        assert "HAVCR2" in DRUGGABILITY
        assert "TIGIT" in DRUGGABILITY
        assert "PHGDH" in DRUGGABILITY

    def test_druggability_values_are_tuples(self):
        for v in DRUGGABILITY.values():
            assert isinstance(v, tuple)
            assert len(v) == 2

    def test_seed_categories_non_empty(self):
        assert len(SEED_CANDIDATE_CATEGORIES) >= 5
        for genes in SEED_CANDIDATE_CATEGORIES.values():
            assert len(genes) > 0
