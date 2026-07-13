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


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    genes = gene_panel()
    log(f"Gene panel: {len(genes)} genes")
    download_string(genes)
    download_cellchatdb(genes)
    download_chea()
    log("Prior-network download complete.")


if __name__ == "__main__":
    main()
