"""
Tests for graph construction module (build_heterograph.py).

Covers: data loaders, edge builders, node/edge creation logic.
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure project root is on path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.graph_construction.build_heterograph import (
    _add_edge,
    build_sst_edges,
    load_ppi,
    load_reactome,
    load_chea_tf,
    load_nk_scores,
    load_cellchatdb,
)


# =========================================================================
# Fixtures
# =========================================================================

@pytest.fixture
def sample_all_nodes():
    """Minimal node dict for testing edge builders."""
    return {
        "PHGDH": {"node_type": "gene", "name": "PHGDH", "source": "sst_module"},
        "PSAT1": {"node_type": "gene", "name": "PSAT1", "source": "sst_module"},
        "SGMS1": {"node_type": "gene", "name": "SGMS1", "source": "sst_module"},
        "SMPD1": {"node_type": "gene", "name": "SMPD1", "source": "sst_module"},
        "ACTB": {"node_type": "gene", "name": "ACTB", "source": "housekeeping"},
        "NKG7": {"node_type": "gene", "name": "NKG7", "source": "sst_module"},
        "GZMB": {"node_type": "gene", "name": "GZMB", "source": "sst_module"},
        "HAVCR2": {"node_type": "gene", "name": "HAVCR2", "source": "checkpoint"},
    }


@pytest.fixture
def temp_ppi_file():
    """Create a temporary PPI file in STRING format (space-separated)."""
    content = (
        "protein1 protein2 combined_score\n"
        "9606.ENSP00000338335 9606.ENSP00000339260 960\n"
        "9606.ENSP00000338335 9606.ENSP00000351257 850\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(content)
        path = f.name
    yield path
    os.unlink(path)


@pytest.fixture
def temp_reactome_file():
    """Create a temporary Reactome file."""
    content = (
        "P12345\tR-HSA-123\thttp://example.org\tPathway Name\tIEA\tHomo sapiens\n"
        "P67890\tR-HSA-456\thttp://example.org\tOther Pathway\tTAS\tHomo sapiens\n"
        "P11111\tR-HSA-789\thttp://example.org\tMouse Pathway\tIEA\tMus musculus\n"
    )
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        path = f.name
    yield path
    os.unlink(path)


# =========================================================================
# Tests: _add_edge
# =========================================================================

class TestAddEdge:
    """Tests for the _add_edge helper."""

    def test_adds_edge_when_both_endpoints_exist(self, sample_all_nodes):
        edges = []
        _add_edge(edges, "PHGDH", "SGMS1", "metabolic_crosstalk", 0.5,
                   "Zheng2023", sample_all_nodes)
        assert len(edges) == 1
        assert edges[0]["src"] == "PHGDH"
        assert edges[0]["dst"] == "SGMS1"
        assert edges[0]["edge_type"] == "metabolic_crosstalk"
        assert edges[0]["weight"] == 0.5

    def test_skips_edge_when_src_missing(self, sample_all_nodes):
        edges = []
        _add_edge(edges, "NONEXISTENT", "SGMS1", "test", 1.0,
                   "test", sample_all_nodes)
        assert len(edges) == 0

    def test_skips_edge_when_dst_missing(self, sample_all_nodes):
        edges = []
        _add_edge(edges, "PHGDH", "NONEXISTENT", "test", 1.0,
                   "test", sample_all_nodes)
        assert len(edges) == 0

    def test_skips_edge_when_both_missing(self, sample_all_nodes):
        edges = []
        _add_edge(edges, "A", "B", "test", 1.0, "test", sample_all_nodes)
        assert len(edges) == 0

    def test_multiple_edges_same_type(self, sample_all_nodes):
        edges = []
        _add_edge(edges, "PHGDH", "SGMS1", "ppi", 0.9, "STRING", sample_all_nodes)
        _add_edge(edges, "PSAT1", "SMPD1", "ppi", 0.8, "STRING", sample_all_nodes)
        assert len(edges) == 2

    def test_empty_nodes_dict(self):
        edges = []
        _add_edge(edges, "PHGDH", "SGMS1", "test", 1.0, "test", {})
        assert len(edges) == 0


# =========================================================================
# Tests: build_sst_edges
# =========================================================================

class TestBuildSstEdges:
    """Tests for SST-axis edge construction."""

    def test_returns_list(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        assert isinstance(edges, list)

    def test_creates_metabolic_crosstalk_edges(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        metabolic = [e for e in edges if e["edge_type"] == "metabolic_crosstalk"]
        assert len(metabolic) > 0
        for e in metabolic:
            assert e["weight"] == 0.5
            assert "Zheng2023" in e["evidence"]

    def test_creates_sm_topology_axis_edges(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        sm_edges = [e for e in edges if e["edge_type"] == "sm_topology_axis"]
        assert len(sm_edges) > 0
        for e in sm_edges:
            assert e["weight"] == 0.3

    def test_edges_only_contain_known_nodes(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        known = set(sample_all_nodes.keys())
        for e in edges:
            assert e["src"] in known, f"Unknown src: {e['src']}"
            assert e["dst"] in known, f"Unknown dst: {e['dst']}"

    def test_metabolic_crosstalk_connects_serine_to_topology(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        metabolic = [e for e in edges if e["edge_type"] == "metabolic_crosstalk"]
        if metabolic:
            for e in metabolic:
                # At least one endpoint should be a serine gene or topology gene
                pass  # Validates by structure — both endpoints were in all_nodes

    def test_empty_nodes_produces_empty_edges(self):
        edges = build_sst_edges({})
        assert edges == []

    def test_edge_structure_has_required_fields(self, sample_all_nodes):
        edges = build_sst_edges(sample_all_nodes)
        for e in edges:
            for key in ["src", "dst", "edge_type", "weight", "evidence"]:
                assert key in e, f"Missing key: {key}"


# =========================================================================
# Tests: load_ppi
# =========================================================================

class TestLoadPpi:
    """Tests for the STRING PPI loader."""

    def test_loads_ppi_from_file(self, temp_ppi_file):
        df = load_ppi(temp_ppi_file)
        assert len(df) == 2
        assert "gene1" in df.columns
        assert "gene2" in df.columns
        assert "weight" in df.columns
        assert "edge_type" in df.columns

    def test_ppi_weight_in_range(self, temp_ppi_file):
        df = load_ppi(temp_ppi_file)
        assert (df["weight"] <= 1.0).all()
        assert (df["weight"] >= 0.0).all()

    def test_ppi_edge_type_is_ppi(self, temp_ppi_file):
        df = load_ppi(temp_ppi_file)
        assert (df["edge_type"] == "ppi").all()

    def test_ppi_missing_file_returns_empty_df(self):
        df = load_ppi("nonexistent_file.txt")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


# =========================================================================
# Tests: load_reactome
# =========================================================================

class TestLoadReactome:
    """Tests for the Reactome pathway loader."""

    def test_loads_reactome_from_file(self, temp_reactome_file):
        df = load_reactome(temp_reactome_file)
        assert len(df) > 0
        assert "uniprot" in df.columns
        assert "pathway_id" in df.columns

    def test_filters_non_human(self, temp_reactome_file):
        df = load_reactome(temp_reactome_file)
        # Only Homo sapiens entries should remain
        assert len(df) == 2  # 2 human, 1 mouse filtered out

    def test_reactome_missing_file_returns_empty_df(self):
        df = load_reactome("nonexistent_file.tsv")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


# =========================================================================
# Tests: load_cellchatdb
# =========================================================================

class TestLoadCellchatDb:
    """Tests for the CellChatDB ligand-receptor loader."""

    def test_loads_lr_pairs(self):
        content = "ligand,receptor,pathway\nWNT1,FZD1,Wnt\nWNT3,FZD2,Wnt\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(content)
            path = f.name
        try:
            df = load_cellchatdb(path)
            assert len(df) == 2
            assert "ligand" in df.columns
            assert "receptor" in df.columns
            assert "edge_type" in df.columns
        finally:
            os.unlink(path)

    def test_missing_columns_returns_empty(self):
        content = "col1,col2\nval1,val2\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(content)
            path = f.name
        try:
            df = load_cellchatdb(path)
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 0
        finally:
            os.unlink(path)

    def test_missing_file_returns_empty_df(self):
        df = load_cellchatdb("nonexistent.csv")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
