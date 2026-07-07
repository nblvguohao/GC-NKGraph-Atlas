"""
Shared test fixtures for GC-NKGraph-Atlas.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import pytest

# Ensure the project root is on sys.path
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))


# ---------------------------------------------------------------------------
# SST gene modules
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_sst_modules() -> Dict[str, Dict]:
    """A minimal SST module dict matching config schema."""
    return {
        "tumor_serine_capacity": {
            "role": "tumor_side",
            "cell_type_attribution": "malignant",
            "genes": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2",
                       "MTHFD1", "MTHFD2", "MTHFD1L", "SLC1A4", "SLC1A5"],
            "expected_direction": "NEEDS_REVIEW",
        },
        "nk_sm_synthesis": {
            "role": "immune_side",
            "cell_type_attribution": "nk",
            "genes": ["SGMS1", "SGMS2"],
            "expected_direction": "higher_is_more_topology_permissive",
        },
        "nk_sm_catabolism": {
            "role": "immune_side",
            "cell_type_attribution": "nk",
            "genes": ["SMPD1", "SMPD2", "SMPD3", "SMPD4"],
            "expected_direction": "higher_is_less_topology_permissive",
        },
        "nk_protrusion_machinery": {
            "role": "immune_side",
            "cell_type_attribution": "nk",
            "genes": ["EZR", "MSN", "RDX", "ACTR2", "ACTR3",
                       "CDC42", "RAC1", "RHOA"],
            "expected_direction": "higher_is_more_topology_permissive",
        },
        "nk_synapse_cytotoxicity_outcome": {
            "role": "outcome_anchor",
            "cell_type_attribution": "nk",
            "genes": ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG"],
            "expected_direction": "axis_positive_correlate",
        },
        "checkpoint_link": {
            "role": "therapeutic_hook",
            "cell_type_attribution": "nk",
            "genes": ["HAVCR2"],
            "expected_direction": "higher_is_less_topology_permissive",
        },
    }


# ---------------------------------------------------------------------------
# Expression data
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_expression_df() -> pd.DataFrame:
    """Small synthetic expression matrix (cells × genes)."""
    rng = np.random.RandomState(42)
    genes = [
        "PHGDH", "PSAT1", "PSPH", "SGMS1", "SGMS2", "SMPD1", "SMPD2",
        "EZR", "MSN", "RDX", "NKG7", "GNLY", "GZMB", "PRF1", "IFNG",
        "HAVCR2", "GAPDH", "ACTB",
    ]
    n_cells = 100
    data = rng.randn(n_cells, len(genes)) + np.arange(len(genes)) * 0.2
    return pd.DataFrame(data, columns=genes, index=[f"cell_{i}" for i in range(n_cells)])


# ---------------------------------------------------------------------------
# Graph fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_graph_dir(tmp_path: Path) -> Path:
    """Create a minimal graph directory with nodes.tsv and edges.tsv."""
    graph_dir = tmp_path / "graph"
    graph_dir.mkdir()

    # Nodes
    nodes = pd.DataFrame({
        "node_id": ["PHGDH", "PSAT1", "SGMS1", "SMPD1", "EZR", "NKG7", "HAVCR2"],
        "node_type": ["gene", "gene", "gene", "gene", "gene", "gene", "gene"],
        "name": ["PHGDH", "PSAT1", "SGMS1", "SMPD1", "EZR", "NKG7", "HAVCR2"],
        "source": ["sst"] * 7,
    })
    nodes.to_csv(graph_dir / "nodes.tsv", sep="\t", index=False)

    # Edges
    edges = pd.DataFrame([
        {"src": "PHGDH", "dst": "PSAT1", "edge_type": "ppi", "weight": 0.9, "evidence": "STRING"},
        {"src": "PHGDH", "dst": "SGMS1", "edge_type": "ppi", "weight": 0.8, "evidence": "STRING"},
        {"src": "SGMS1", "dst": "SMPD1", "edge_type": "ppi", "weight": 0.7, "evidence": "STRING"},
        {"src": "PHGDH", "dst": "SGMS1", "edge_type": "metabolic_crosstalk", "weight": 0.5, "evidence": "Zheng2023"},
        {"src": "EZR", "dst": "NKG7", "edge_type": "sm_topology_axis", "weight": 0.3, "evidence": "Zheng2023"},
    ])
    edges.to_csv(graph_dir / "edges.tsv", sep="\t", index=False)

    return graph_dir


# ---------------------------------------------------------------------------
# NK state labels
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_nk_labels() -> pd.DataFrame:
    """Synthetic NK state labels."""
    rng = np.random.RandomState(42)
    states = ["NK-hot-cytotoxic", "NK-hot-dysfunctional", "NK-cold/excluded", "NK-intermediate"]
    samples = [f"TCGA-{i:02d}" for i in range(80)]
    return pd.DataFrame({
        "nk_immune_state": rng.choice(states, size=len(samples)),
    }, index=samples)
