"""
GC-NKGraph-Atlas Heterogeneous Graph Construction (Phase 8).

Integrates PPI, ligand-receptor, pathway, TF-target, SST axis, and
NK state edges into a typed heterogeneous graph.

Output:
  - nodes.tsv: node_id, node_type, name, source
  - edges.tsv: src, dst, edge_type, weight, evidence
  - node_features.parquet: numerical features per node
  - edge_features.parquet: edge weights and types

Usage:
    python src/graph_construction/build_heterograph.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.log_utils import Logger  # noqa: E402
from src.common.io_utils import ensure_dir, load_table, load_config  # noqa: E402
from src.common.sst_config import load_sst_modules, get_sst_genes  # noqa: E402

logger = Logger()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --- Load SST gene modules from shared config (single source of truth) ---
SST_MODULES = load_sst_modules()


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_ppi(path: str) -> pd.DataFrame:
    """Load STRING PPI as edges. Returns DataFrame with gene1, gene2, weight, edge_type."""
    log(f"  Loading PPI from {path}...")
    try:
        f = open(path, "rt")
        header = f.read(2)
        f.close()

        is_gzip = header == "\x1f\x8b"
        if is_gzip or path.endswith(".gz"):
            import gzip
            f = gzip.open(path, "rt")
        else:
            f = open(path, "rt")

        df = pd.read_csv(f, sep=" ", comment="#")
        f.close()
        df.columns = ["protein1", "protein2", "score"]
        df["edge_type"] = "ppi"
        df["weight"] = df["score"].astype(float) / 1000.0

        # Extract gene symbols
        df["gene1"] = (
            df["protein1"].str.split(".").str[-1]
            if df["protein1"].astype(str).str.contains(r"\.", na=False).any()
            else df["protein1"]
        )
        df["gene2"] = (
            df["protein2"].str.split(".").str[-1]
            if df["protein2"].astype(str).str.contains(r"\.", na=False).any()
            else df["protein2"]
        )
        log(f"    {len(df)} raw edges")
        return df[["gene1", "gene2", "weight", "edge_type"]]
    except Exception as e:
        log(f"    PPI load failed: {e}")
        return pd.DataFrame()


def load_reactome(path: str) -> pd.DataFrame:
    """Load Reactome pathways (UniProt-to-Pathway)."""
    log(f"  Loading Reactome from {path}...")
    try:
        cols = ["uniprot", "pathway_id", "url", "pathway_name", "evidence", "species"]
        df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str, on_bad_lines="skip")
        df = df[df["species"].fillna("").str.contains("Homo sapiens", na=False)]
        log(f"    {len(df)} human gene-pathway associations")
        return df[["uniprot", "pathway_id", "pathway_name"]]
    except Exception as e:
        log(f"    Reactome load failed: {e}")
        return pd.DataFrame()


def load_chea_tf(path: str) -> pd.DataFrame:
    """Load ChEA TF-target relationships.

    Enrichr's ChEA gene-set-library export is ragged (one TF per line, with a
    variable number of target genes), so it cannot be read with
    `pd.read_csv` under a fixed column count. Parse it line by line instead:
    column 0 becomes the TF symbol (first whitespace token of the
    description field) and columns 1..N become that TF's target genes,
    NaN-padded to a rectangular frame for the caller's `iterrows()` loop.
    """
    log(f"  Loading TF-target from {path}...")
    try:
        rows = []
        with open(path, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0].strip():
                    continue
                tf = parts[0].split()[0].strip()
                targets = [t.strip() for t in parts[1:] if t.strip()]
                rows.append([tf] + targets)
        df = pd.DataFrame(rows)
        log(f"    {len(df)} TF-target edges")
        return df
    except Exception as e:
        log(f"    TF-target load failed: {e}")
        return pd.DataFrame()


def load_nk_scores(path: str) -> pd.DataFrame:
    """Load NK state scores for sample nodes."""
    try:
        df = load_table(path)
        log(f"    {len(df)} samples with NK scores")
        return df
    except Exception as e:
        log(f"    NK scores load failed: {e}")
        return pd.DataFrame()


def load_cellchatdb(path: str) -> pd.DataFrame:
    """Load CellChatDB ligand-receptor pairs."""
    log(f"  Loading CellChatDB from {path}...")
    try:
        df = pd.read_csv(path)
        if "ligand" in df.columns and "receptor" in df.columns:
            df["edge_type"] = "ligand_receptor"
            log(f"    {len(df)} LR pairs")
            return df
        return pd.DataFrame()
    except Exception as e:
        log(f"    CellChatDB failed: {e}")
        return pd.DataFrame()


def load_go_bp(path: str) -> pd.DataFrame:
    """Load GO Biological Process 2023 gene sets (Enrichr export).

    Same ragged one-term-per-line format as ChEA; reuses that parser.
    """
    log(f"  Loading GO_BP from {path}...")
    return load_chea_tf(path)


def build_geneset_membership(
    library_df: pd.DataFrame,
    panel_genes: Set[str],
) -> Dict[str, Set[str]]:
    """Invert a ragged term->genes table into gene->{terms} for panel genes only.

    `library_df` has one row per term: column 0 is the term id/name, remaining
    columns (NaN-padded) are member genes.
    """
    gene_terms: Dict[str, Set[str]] = {g: set() for g in panel_genes}
    for _, row in library_df.iterrows():
        vals = [str(v).strip() for v in row.values if str(v) != "nan"]
        if len(vals) < 2:
            continue
        term_id, genes = vals[0], vals[1:]
        for g in genes:
            gu = g.upper().strip()
            if gu in gene_terms:
                gene_terms[gu].add(term_id)
    return gene_terms


def load_msigdb_c2(npz_path: str, genelist_path: str, panel_genes: Set[str]) -> Dict[str, Set[int]]:
    """Load the MSigDB C2 gene x gene-set sparse matrix, restricted to panel genes.

    Returns gene -> set of gene-set column indices the gene belongs to.
    """
    log(f"  Loading MSigDB C2 from {npz_path}...")
    try:
        from scipy import sparse
        mat = sparse.load_npz(npz_path).tocsr()
        with open(genelist_path, encoding="utf-8") as f:
            gene_list = [line.strip() for line in f if line.strip()]
        gene_sets: Dict[str, Set[int]] = {}
        for gene, row_idx in zip(gene_list, range(mat.shape[0])):
            gu = gene.upper().strip()
            if gu in panel_genes:
                gene_sets[gu] = set(mat.indices[mat.indptr[row_idx]:mat.indptr[row_idx + 1]].tolist())
        log(f"    {len(gene_sets)}/{len(panel_genes)} panel genes found in MSigDB C2 ({mat.shape[1]} gene sets)")
        return gene_sets
    except Exception as e:
        log(f"    MSigDB C2 load failed: {e}")
        return {}


def build_jaccard_edges(
    gene_sets: Dict[str, Set],
    edge_type: str,
    evidence: str,
    all_nodes: Dict,
    jaccard_threshold: float = 0.05,
    weight: float = 0.2,
) -> List[Dict]:
    """Build similarity edges between genes sharing >= jaccard_threshold of their sets.

    Mirrors TreeNet's GO-prior edge augmentation: a fixed conservative edge
    weight (default 0.2) is used regardless of similarity magnitude, so the
    prior nudges message passing without dominating mechanism-grounded edges.
    Genes with empty sets (no annotation found in the library) are skipped.
    """
    edges: List[Dict] = []
    genes = sorted(g for g, s in gene_sets.items() if s)
    for i in range(len(genes)):
        si = gene_sets[genes[i]]
        for j in range(i + 1, len(genes)):
            sj = gene_sets[genes[j]]
            union = len(si | sj)
            if union == 0:
                continue
            jac = len(si & sj) / union
            if jac >= jaccard_threshold:
                _add_edge(edges, genes[i], genes[j], edge_type, weight,
                           f"{evidence}_jaccard{jac:.3f}", all_nodes)
    return edges


# ---------------------------------------------------------------------------
# Edge builders
# ---------------------------------------------------------------------------

def _add_edge(
    edges: List[Dict],
    src: str,
    dst: str,
    edge_type: str,
    weight: float,
    evidence: str,
    all_nodes: Dict,
) -> None:
    """Add an edge if both endpoints are in the node set."""
    if src in all_nodes and dst in all_nodes:
        edges.append({
            "src": src, "dst": dst,
            "edge_type": edge_type,
            "weight": weight,
            "evidence": evidence,
        })


def build_sst_edges(
    all_nodes: Dict[str, Dict],
    modules: Optional[Dict[str, Dict]] = None,
) -> List[Dict]:
    """Build SST-axis specific edges (metabolic_crosstalk + sm_topology_axis).

    Args:
        all_nodes: Dict of node_id -> {node_type, name, source}.
        modules: SST gene modules (default: loaded from config).

    Returns:
        List of edge dicts.
    """
    if modules is None:
        modules = load_sst_modules()

    edges: List[Dict] = []

    # metabolic_crosstalk: tumor_serine_program 鈫?nk_topology_state
    tumor_serine_genes: List[str] = modules.get("tumor_serine_capacity", {}).get("genes", [])
    nk_topology_genes: List[str] = (
        modules.get("nk_sm_synthesis", {}).get("genes", [])
        + modules.get("nk_sm_catabolism", {}).get("genes", [])
        + modules.get("nk_protrusion_machinery", {}).get("genes", [])
    )
    for tg in tumor_serine_genes:
        for ng in nk_topology_genes:
            _add_edge(edges, tg, ng, "metabolic_crosstalk", 0.5,
                       "Zheng2023_NatImmunol_SST_axis", all_nodes)

    # sm_topology_axis: within-axis co-expression edges
    axis_genes: Set[str] = set(nk_topology_genes) | set(
        modules.get("nk_synapse_cytotoxicity_outcome", {}).get("genes", [])
    )
    axis_list = list(axis_genes)
    for i in range(len(axis_list)):
        for j in range(i + 1, len(axis_list)):
            _add_edge(edges, axis_list[i], axis_list[j],
                       "sm_topology_axis", 0.3, "Zheng2023_NatImmunol", all_nodes)

    return edges


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="data/processed/graph",
                         help="Output directory (default: data/processed/graph, "
                              "the graph the manuscript reports). Use a separate "
                              "directory for exploratory variants.")
    parser.add_argument("--enable-go-prior", action="store_true",
                         help="Add GO Biological Process co-membership edges "
                              "(experimental, off by default; see Methods 2.5).")
    parser.add_argument("--enable-msigdb-prior", action="store_true",
                         help="Add MSigDB C2 co-membership edges "
                              "(experimental, off by default; see Methods 2.5).")
    args = parser.parse_args()

    log("=" * 60)
    log("HETEROGENEOUS GRAPH CONSTRUCTION (Phase 8)")
    log("=" * 60)

    out_dir = ensure_dir(args.out_dir)
    prior_dir = "data/raw/prior_networks"

    all_nodes: Dict[str, Dict] = {}
    all_edges: List[Dict] = []

    # ---- Step 1: Gene nodes from SST modules + candidate pool ----
    log("\nStep 1: Gene nodes...")
    all_sst_genes: Set[str] = set()
    for module_name, mod_data in SST_MODULES.items():
        for g in mod_data["genes"]:
            all_sst_genes.add(g)
            all_nodes[g] = {
                "node_type": "gene",
                "name": g,
                "source": f"sst_axis_{module_name}",
            }

    # Known NK receptors
    nk_receptors = [
        "KLRK1", "NCR1", "NCR2", "NCR3", "FCGR3A", "CD96", "TIGIT",
        "KLRC1", "KLRC2", "KLRD1", "CD226", "SLAMF6", "SLAMF7",
    ]
    for g in nk_receptors:
        all_nodes[g] = {"node_type": "nk_receptor", "name": g, "source": "nk_receptor_list"}
    all_sst_genes.update(nk_receptors)

    # Add genes from candidate evidence matrix
    cand_path = "results/tables/candidate_evidence_matrix.tsv"
    if os.path.exists(cand_path):
        cand = pd.read_csv(cand_path, sep="\t")
        for _, r in cand.iterrows():
            g = str(r["gene"])
            if g not in all_nodes:
                all_nodes[g] = {"node_type": "gene", "name": g, "source": "candidate_pool"}

    log(f"  Gene nodes: {len(all_nodes)}")

    # ---- Step 2: Prior network edges ----
    log("\nStep 2: Prior network edges...")

    # PPI (STRING) 鈥?high confidence only
    ppi = load_ppi(os.path.join(prior_dir, "string_ppi_physical.tsv"))
    ppi_filtered = ppi[ppi["weight"] >= 0.7] if len(ppi) > 0 else ppi
    for _, r in ppi_filtered.iterrows():
        _add_edge(all_edges, str(r["gene1"]), str(r["gene2"]),
                   "ppi", float(r["weight"]), "STRING_v12", all_nodes)
    log(f"    {len(ppi_filtered)} PPI edges (score>=0.7)")

    # Ligand-Receptor (CellChatDB)
    ccdb_path = os.path.join(prior_dir, "cellchatdb_human.csv")
    if os.path.exists(ccdb_path) and os.path.getsize(ccdb_path) > 0:
        lr = load_cellchatdb(ccdb_path)
        for _, r in lr.iterrows():
            _add_edge(all_edges, str(r.get("ligand", "")), str(r.get("receptor", "")),
                       "ligand_receptor", 0.9, "CellChatDB", all_nodes)

    # Pathway (Reactome)
    react_path = os.path.join(prior_dir, "reactome_human_genes.txt")
    if os.path.exists(react_path):
        react = load_reactome(react_path)
        top_paths = react[~react["pathway_id"].str.contains("-", na=False)]
        log(f"    {len(top_paths)} top-level pathway associations")

    # TF-target (ChEA)
    chea_path = os.path.join(prior_dir, "chea_tf_targets.txt")
    if os.path.exists(chea_path):
        chea = load_chea_tf(chea_path)
        tf_edges = 0
        for _, r in chea.iterrows():
            cols = [str(c).upper().strip() for c in r.values if str(c) != "nan"]
            if len(cols) >= 2:
                tf = cols[0]
                for target in cols[1:]:
                    if tf in all_nodes and target in all_nodes and tf != target:
                        all_edges.append({
                            "src": tf, "dst": target,
                            "edge_type": "tf_target",
                            "weight": 0.8,
                            "evidence": "ChEA_2022",
                        })
                        tf_edges += 1
        log(f"    {tf_edges} TF-target edges")

    # GO-prior / MSigDB-prior co-membership edges (experimental, opt-in) --
    # panel = all gene / nk_receptor nodes registered so far (Step 1).
    panel_genes = {
        nid for nid, props in all_nodes.items()
        if props["node_type"] in ("gene", "nk_receptor")
    }

    if args.enable_go_prior:
        go_path = os.path.join(prior_dir, "go_bp_2023.txt")
        if os.path.exists(go_path):
            go_lib = load_go_bp(go_path)
            gene_terms = build_geneset_membership(go_lib, panel_genes)
            go_edges = build_jaccard_edges(
                gene_terms, "go_prior", "GO_BP_2023", all_nodes,
                jaccard_threshold=0.05, weight=0.2,
            )
            all_edges.extend(go_edges)
            log(f"    {len(go_edges)} GO-prior edges")
        else:
            log(f"    GO_BP file not found at {go_path} (run download_prior_networks.py --include-go-msigdb)")

    if args.enable_msigdb_prior:
        npz_path = os.path.join(prior_dir, "c2_GenesetsMatrix.npz")
        genelist_path = os.path.join(prior_dir, "geneList.csv")
        if os.path.exists(npz_path) and os.path.exists(genelist_path):
            gene_sets = load_msigdb_c2(npz_path, genelist_path, panel_genes)
            msigdb_edges = build_jaccard_edges(
                gene_sets, "msigdb_prior", "MSigDB_C2", all_nodes,
                jaccard_threshold=0.05, weight=0.2,
            )
            all_edges.extend(msigdb_edges)
            log(f"    {len(msigdb_edges)} MSigDB-prior edges")
        else:
            log(f"    MSigDB C2 files not found (run download_prior_networks.py --include-go-msigdb)")

    # ---- Step 3: SST-axis edges (from shared config) ----
    log("\nStep 3: SST-axis edges...")
    sst_edges = build_sst_edges(all_nodes, SST_MODULES)
    all_edges.extend(sst_edges)
    log(f"    {len(sst_edges)} SST-axis edges")

    # ---- Step 4: NK state nodes and edges ----
    log("\nStep 4: NK state nodes...")
    nk_scores = load_nk_scores("results/tables/nk_scores_bulk.tsv")
    if len(nk_scores) > 0:
        for state in [
            "NK-hot-cytotoxic", "NK-hot-dysfunctional",
            "NK-cold/excluded", "NK-intermediate",
        ]:
            state_id = f"state_{state.lower().replace('-', '_').replace('/', '_')}"
            all_nodes[state_id] = {
                "node_type": "cell_state",
                "name": state,
                "source": "nk_scoring",
            }

        # Add dysfunction correlation edges for SST genes
        for g in all_sst_genes:
            if g in nk_scores.columns:
                dysf_col = nk_scores.get("NK_dysfunction_score", pd.Series(dtype=float))
                if len(dysf_col) > 0:
                    corr = float(nk_scores[g].corr(dysf_col)) if g in nk_scores.columns else 0.0
                    if abs(corr) > 0.1:
                        all_edges.append({
                            "src": g,
                            "dst": "state_nk_hot_dysfunctional",
                            "edge_type": "dysfunction_correlation",
                            "weight": abs(corr),
                            "evidence": f"bulk_corr_{corr:.3f}",
                        })

    # ---- Step 5: Write outputs ----
    log("\nStep 5: Writing outputs...")

    node_df = pd.DataFrame([
        {"node_id": nid, **props} for nid, props in all_nodes.items()
    ])
    node_path = os.path.join(out_dir, "nodes.tsv")
    node_df.to_csv(node_path, sep="\t", index=False)
    log(f"  Nodes: {len(node_df)} ({node_df['node_type'].value_counts().to_dict()})")

    edge_df = pd.DataFrame(all_edges)
    edge_path = os.path.join(out_dir, "edges.tsv")
    edge_df.to_csv(edge_path, sep="\t", index=False)
    log(f"  Edges: {len(edge_df)} ({edge_df['edge_type'].value_counts().to_dict()})")

    # Graph statistics
    stats = {
        "total_nodes": len(node_df),
        "total_edges": len(edge_df),
        "node_types": node_df["node_type"].value_counts().to_dict(),
        "edge_types": edge_df["edge_type"].value_counts().to_dict(),
        "density": (
            len(edge_df) / (len(node_df) * (len(node_df) - 1))
            if len(node_df) > 1 else 0.0
        ),
    }
    stats_path = os.path.join(out_dir, "graph_statistics.tsv")
    pd.DataFrame([stats]).to_csv(stats_path, sep="\t", index=False)
    log("  Graph statistics saved")

    log("\n" + "=" * 60)
    log("GRAPH CONSTRUCTION COMPLETE!")
    log("=" * 60)


if __name__ == "__main__":
    main()
