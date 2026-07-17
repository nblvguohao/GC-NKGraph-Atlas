"""Calibrate mechanism-view module coupling against node-label permutations."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Mapping, Set, Tuple

import pandas as pd

from src.common.sst_config import load_sst_modules
from src.models.multiview_fusion import (
    DEFAULT_VIEW_SPEC,
    build_graph_views,
    calibrate_view_in_context,
)


def _genes(modules: Mapping[str, Mapping], name: str) -> Set[str]:
    return {str(gene) for gene in modules.get(name, {}).get("genes", [])}


def calibration_module_pairs(
    modules: Mapping[str, Mapping],
) -> Dict[str, Tuple[Set[str], Set[str]]]:
    """Return the pre-registered source/target modules for each mechanism view."""

    topology = (
        _genes(modules, "nk_sm_synthesis")
        | _genes(modules, "nk_sm_catabolism")
        | _genes(modules, "nk_protrusion_machinery")
    )
    return {
        "metabolic_crosstalk": (_genes(modules, "tumor_serine_capacity"), topology),
        "sm_topology_axis": (topology, _genes(modules, "nk_synapse_cytotoxicity_outcome")),
    }


def run_calibration(
    graph_dir: Path,
    output_dir: Path,
    *,
    embedding_dim: int = 64,
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(graph_dir / "nodes.tsv", sep="\t")
    edges = pd.read_csv(graph_dir / "edges.tsv", sep="\t")
    views = build_graph_views(edges, nodes, DEFAULT_VIEW_SPEC, strict=True, require_nonempty=True)
    pairs = calibration_module_pairs(load_sst_modules())
    summary_rows = []
    null_rows = []
    for view_name, (source_genes, target_genes) in pairs.items():
        summary, null_values = calibrate_view_in_context(
            views,
            randomized_view_name=view_name,
            source_genes=source_genes,
            target_genes=target_genes,
            embedding_dim=embedding_dim,
            n_randomizations=n_randomizations,
            seed=seed,
        )
        summary_rows.append(
            {
                **summary,
                "source_module_genes": ";".join(sorted(source_genes)),
                "target_module_genes": ";".join(sorted(target_genes)),
            }
        )
        null_rows.extend(
            {"view": view_name, "randomization": index, "coupling": float(value)}
            for index, value in enumerate(null_values, start=1)
        )
    summary_table = pd.DataFrame(summary_rows)
    null_table = pd.DataFrame(null_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_table.to_csv(
        output_dir / "multiview_mechanism_randomization_summary.tsv", sep="\t", index=False
    )
    null_table.to_csv(
        output_dir / "multiview_mechanism_randomization_null.tsv", sep="\t", index=False
    )
    provenance = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "graph_dir": str(graph_dir),
        "embedding_dim": embedding_dim,
        "n_randomizations": n_randomizations,
        "seed": seed,
        "interpretation_boundary": (
            "Topology-specific calibration only; this is not biological or predictive validation."
        ),
    }
    (output_dir / "multiview_calibration_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    return summary_table, null_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--embedding-dim", type=int, default=64)
    parser.add_argument("--n-randomizations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    summary, _ = run_calibration(
        args.graph_dir,
        args.output_dir,
        embedding_dim=args.embedding_dim,
        n_randomizations=args.n_randomizations,
        seed=args.seed,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
