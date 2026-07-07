"""
Tests for src/common/sst_config.py — the single source of truth for SST gene modules.
"""

from __future__ import annotations

import pytest
from src.common.sst_config import (
    load_sst_modules,
    get_sst_genes,
    get_module_for_gene,
    get_genes_by_role,
    SST_MODULES,
    SST_GENES,
)


class TestLoadSstModules:
    """Verify the config loader returns correct structure."""

    def test_loads_from_default_config(self):
        """Default load returns all expected module types."""
        modules = load_sst_modules()
        assert isinstance(modules, dict)
        assert len(modules) >= 5, f"Expected >=5 modules, got {len(modules)}"

    def test_each_module_has_required_keys(self):
        """Every module must have role, cell_type_attribution, genes, expected_direction."""
        modules = load_sst_modules()
        required = {"role", "cell_type_attribution", "genes", "expected_direction"}
        for name, mod in modules.items():
            missing = required - set(mod.keys())
            assert not missing, f"Module '{name}' missing keys: {missing}"

    def test_genes_are_non_empty_lists(self):
        """No module should have empty gene lists."""
        modules = load_sst_modules()
        for name, mod in modules.items():
            assert mod["genes"], f"Module '{name}' has empty gene list"
            assert all(isinstance(g, str) for g in mod["genes"]), \
                f"Module '{name}' has non-string gene entries"

    def test_known_gene_present(self):
        """PHGDH must be in tumor_serine_capacity."""
        modules = load_sst_modules()
        assert "PHGDH" in modules["tumor_serine_capacity"]["genes"]

    def test_havcr2_in_checkpoint_link(self):
        """HAVCR2 must be the checkpoint_link gene."""
        modules = load_sst_modules()
        assert "HAVCR2" in modules["checkpoint_link"]["genes"]

    def test_caching_returns_same_object(self):
        """Repeated calls with no args return the cached modules."""
        m1 = load_sst_modules()
        m2 = load_sst_modules()
        assert m1 is m2


class TestGetSstGenes:
    """Verify flat gene set extraction."""

    def test_returns_set(self):
        genes = get_sst_genes()
        assert isinstance(genes, set)

    def test_contains_core_genes(self):
        genes = get_sst_genes()
        assert "PHGDH" in genes
        assert "SGMS1" in genes
        assert "NKG7" in genes
        assert "HAVCR2" in genes

    def test_no_duplicates(self, sample_sst_modules):
        """Genes appearing in multiple modules should not be duplicated."""
        genes = get_sst_genes(sample_sst_modules)
        assert len(genes) == sum(
            len(set(m["genes"])) for m in sample_sst_modules.values()
        )  # all unique in fixture


class TestGetModuleForGene:
    """Verify gene-to-module mapping."""

    def test_maps_phgdh_to_serine(self):
        assert get_module_for_gene("PHGDH") == "tumor_serine_capacity"

    def test_maps_nkg7_to_cytotoxicity(self):
        assert get_module_for_gene("NKG7") == "nk_synapse_cytotoxicity_outcome"

    def test_unknown_gene_returns_none(self):
        assert get_module_for_gene("NOT_A_GENE_XYZ") is None


class TestGetGenesByRole:
    """Verify role-based gene filtering."""

    def test_tumor_side_returns_serine_genes(self):
        genes = get_genes_by_role("tumor_side")
        assert "PHGDH" in genes
        assert "PSAT1" in genes

    def test_outcome_anchor_returns_cytotoxicity_genes(self):
        genes = get_genes_by_role("outcome_anchor")
        assert "NKG7" in genes
        assert "GZMB" in genes

    def test_unknown_role_returns_empty(self):
        assert get_genes_by_role("nonexistent_role") == []


class TestModuleConstants:
    """Verify the convenience SST_MODULES and SST_GENES constants."""

    def test_sst_modules_is_dict(self):
        assert isinstance(SST_MODULES, dict)
        assert len(SST_MODULES) > 0

    def test_sst_genes_is_set(self):
        assert isinstance(SST_GENES, set)
        assert len(SST_GENES) > 10
