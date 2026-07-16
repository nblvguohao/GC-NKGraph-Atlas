"""
Download the prior-network data sources for the heterogeneous graph (Phase 8).

The graph in build_heterograph.py consumes three external prior networks:
  - STRING v12 protein-protein interactions (score >= 700), restricted to the
    axis-centered gene panel, via the STRING REST API.
  - CellChatDB ligand-receptor pairs, via the OmniPath REST API (which
    aggregates CellChatDB among other resources).
  - ChEA 2022 transcription-factor -> target sets, via the Enrichr gene-set
    library export.

These files are NOT committed (data/raw/ is gitignored), so this script
regenerates them from a clean clone. Run it before build_heterograph.py:

    python src/data_download/download_prior_networks.py
    python -m src.graph_construction.build_heterograph

Networks: STRING and OmniPath are queried for exactly the graph's gene panel,
so the STRING file already reflects the score>=700 threshold and the panel
restriction. ChEA is downloaded whole; build_heterograph.py restricts it to the
panel at graph-build time (and, because the panel is ~100 axis genes, no ChEA
pair currently falls entirely within it -- see Methods 2.5).

Optional (--include-go-msigdb, off by default): two additional generic-prior
sources, added as an exploratory robustness check on whether GO/MSigDB
co-membership edges can substitute for the ChEA edge type that currently
contributes 0 edges within the axis panel:
  - GO Biological Process 2023 (Enrichr gene-set-library export, same ragged
    format as ChEA) -- mirrors TreeNet's "GO-prior" edge-augmentation idea.
  - MSigDB C2 (curated pathways) gene x gene-set membership matrix, reused
    verbatim from github.com/nblvguohao/CANOPY-Router (data/msigdb/), a
    separate unpublished project by the same author -- NOT an independent
    external source, just a convenient pre-packaged MSigDB C2 snapshot.
These are experimental additions, not part of the graph reported in the
manuscript; build_heterograph.py only uses them when explicitly enabled via
--enable-go-prior / --enable-msigdb-prior.
"""
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.graph_construction.build_heterograph import SST_MODULES  # noqa: E402

OUT_DIR = Path("data/raw/prior_networks")
NK_RECEPTORS = [
    "KLRK1", "NCR1", "NCR2", "NCR3", "FCGR3A", "CD96", "TIGIT",
    "KLRC1", "KLRC2", "KLRD1", "CD226", "SLAMF6", "SLAMF7",
]
CAND_PATH = "results/tables/candidate_evidence_matrix.tsv"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def gene_panel() -> list:
    """Reconstruct the graph's gene-node universe (mirrors build_heterograph main())."""
    genes = set()
    for mod in SST_MODULES.values():
        genes.update(mod["genes"])
    genes.update(NK_RECEPTORS)
    if os.path.exists(CAND_PATH):
        genes.update(pd.read_csv(CAND_PATH, sep="\t")["gene"].astype(str))
    return sorted(genes)


def download_string(genes: list) -> None:
    """STRING v12 PPI, score>=700, restricted to the gene panel."""
    log(f"STRING: querying {len(genes)} genes (score>=700)...")
    ids = "%0d".join(genes)
    url = (
        "https://string-db.org/api/tsv/network"
        f"?identifiers={ids}&species=9606&required_score=700"
    )
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    from io import StringIO
    df = pd.read_csv(StringIO(r.text), sep="\t")
    out = df[["preferredName_A", "preferredName_B", "score"]].copy()
    out["score"] = (out["score"] * 1000).round(1)  # STRING API returns 0-1 here
    out.columns = ["protein1", "protein2", "score"]
    dest = OUT_DIR / "string_ppi_physical.tsv"
    out.to_csv(dest, sep=" ", index=False)
    log(f"STRING: wrote {len(out)} edges -> {dest}")


def download_cellchatdb(genes: list) -> None:
    """CellChatDB ligand-receptor pairs via OmniPath (kept whole; the graph
    restricts to the panel at build time)."""
    log("CellChatDB (OmniPath): querying...")
    url = (
        "https://omnipathdb.org/interactions"
        "?resources=CellChatDB&fields=sources&genesymbols=1"
    )
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    from io import StringIO
    df = pd.read_csv(StringIO(r.text), sep="\t")
    out = df[["source_genesymbol", "target_genesymbol"]].copy()
    out.columns = ["ligand", "receptor"]
    out = out.drop_duplicates()
    dest = OUT_DIR / "cellchatdb_human.csv"
    out.to_csv(dest, index=False)
    n_panel = ((out["ligand"].isin(genes)) & (out["receptor"].isin(genes))).sum()
    log(f"CellChatDB: wrote {len(out)} LR pairs ({n_panel} within panel) -> {dest}")


def download_chea() -> None:
    """ChEA 2022 TF-target sets from the Enrichr gene-set library."""
    log("ChEA 2022 (Enrichr): downloading library...")
    url = (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary"
        "?mode=text&libraryName=ChEA_2022"
    )
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    dest = OUT_DIR / "chea_tf_targets.txt"
    dest.write_text(r.text, encoding="utf-8")
    n_sets = sum(1 for line in r.text.splitlines() if line.strip())
    log(f"ChEA: wrote {n_sets} TF sets -> {dest}")


def download_go_bp() -> None:
    """GO Biological Process 2023 gene sets from the Enrichr gene-set library.

    Same ragged text format as ChEA (download_chea); build_heterograph.py
    restricts it to the gene panel and derives GO-prior similarity edges from
    term co-membership.
    """
    log("GO_Biological_Process_2023 (Enrichr): downloading library...")
    url = (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary"
        "?mode=text&libraryName=GO_Biological_Process_2023"
    )
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    dest = OUT_DIR / "go_bp_2023.txt"
    dest.write_text(r.text, encoding="utf-8")
    n_terms = sum(1 for line in r.text.splitlines() if line.strip())
    log(f"GO_BP: wrote {n_terms} terms -> {dest}")


def download_msigdb_c2() -> None:
    """MSigDB C2 (curated pathways) gene x gene-set matrix.

    Reused as-is from github.com/nblvguohao/CANOPY-Router (a separate,
    unpublished project by the same author) rather than re-derived from raw
    MSigDB, purely as a convenient pre-packaged snapshot. Treat as an internal
    convenience resource, not an independently sourced external dataset.
    """
    base = "https://raw.githubusercontent.com/nblvguohao/CANOPY-Router/main/data/msigdb"
    log("MSigDB C2 (via CANOPY-Router bundle): downloading...")
    for fname in ("c2_GenesetsMatrix.npz", "geneList.csv"):
        r = requests.get(f"{base}/{fname}", timeout=120)
        r.raise_for_status()
        dest = OUT_DIR / fname
        dest.write_bytes(r.content)
        log(f"  wrote {len(r.content)} bytes -> {dest}")


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--include-go-msigdb", action="store_true",
        help="Also fetch GO_BP_2023 and MSigDB C2 (off by default; not part "
             "of the reported manuscript graph, see module docstring).",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    genes = gene_panel()
    log(f"Gene panel: {len(genes)} genes")
    download_string(genes)
    download_cellchatdb(genes)
    download_chea()
    if args.include_go_msigdb:
        download_go_bp()
        download_msigdb_c2()
    log("Prior-network download complete.")


if __name__ == "__main__":
    main()
