"""
Build gene/probe ID → HGNC symbol mappings for all platforms.

Platforms:
  - GPL570 (Affymetrix HGU133 Plus 2.0) — GSE62254
  - GPL10558 (Illumina HumanHT-12) — GSE84437
  - Ensembl gene IDs — TCGA-STAD/LIHC
"""

import os
import sys
import json
import gzip
import urllib.request
from pathlib import Path
import pandas as pd
import numpy as np

os.environ["no_proxy"] = "*"
os.environ["all_proxy"] = ""
os.environ["ALL_PROXY"] = ""

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir

ANNOT_DIR = "data/raw/annotations"
ensure_dir(ANNOT_DIR)


def build_illumina_mapping():
    """Parse GPL10558 annotation for Illumina probe→symbol mapping."""
    print("=== GPL10558 (Illumina) ===")
    annot_path = f"{ANNOT_DIR}/GPL10558.annot.gz"

    if not os.path.exists(annot_path):
        print("  Annotation file not found!")
        return {}

    try:
        with gzip.open(annot_path, "rt") as f:
            df = pd.read_csv(f, sep="\t", comment="#", low_memory=False)
        print(f"  Columns: {list(df.columns[:12])}")

        # Find symbol column
        symbol_col = None
        for c in df.columns:
            if "symbol" in c.lower():
                symbol_col = c
                break

        if not symbol_col:
            # Try common Illumina annotation column names
            for c in ["Symbol", "Gene Symbol", "GENE_SYMBOL", "gene_symbol"]:
                if c in df.columns:
                    symbol_col = c
                    break

        illumina_map = {}
        if symbol_col and "ID" in df.columns:
            for _, row in df.iterrows():
                probe = str(row["ID"]).strip()
                symbol = str(row[symbol_col]).strip()
                if symbol and symbol != "nan" and symbol != "":
                    symbol = symbol.split("///")[0].strip()
                    symbol = symbol.split("//")[0].strip()
                    if symbol:
                        illumina_map[probe] = symbol
        else:
            print(f"  WARNING: Cannot find ID/Symbol columns")
            print(f"  Available: {list(df.columns)}")

        print(f"  Mapped: {len(illumina_map)} probes")

        # Save
        with open(f"{ANNOT_DIR}/gpl10558_probe_to_symbol.tsv", "w") as f:
            f.write("probe_id\tgene_symbol\n")
            for p, s in sorted(illumina_map.items()):
                f.write(f"{p}\t{s}\n")
        return illumina_map

    except Exception as e:
        print(f"  FAILED: {e}")
        return {}


def build_ensembl_mapping():
    """Query mygene.info API for Ensembl→Symbol mapping."""
    print("\n=== Ensembl → Symbol (mygene.info) ===")

    # Get all Ensembl IDs from TCGA data
    all_ensembl = set()
    for fname in [
        "data/processed/bulk/tcga_stad_expression.tsv",
        "data/processed/bulk/tcga_lihc_expression.tsv",
    ]:
        if os.path.exists(fname):
            header = pd.read_csv(fname, sep="\t", index_col=0, nrows=0)
            for col in header.columns:
                all_ensembl.add(col.split(".")[0])  # Strip version

    print(f"  Total unique Ensembl IDs: {len(all_ensembl)}")

    if len(all_ensembl) == 0:
        print("  No Ensembl IDs found!")
        return {}

    ensembl_map = {}
    clean_list = list(all_ensembl)
    batch_size = 1000
    url = "https://mygene.info/v3/query"

    for i in range(0, len(clean_list), batch_size):
        batch = clean_list[i : i + batch_size]
        try:
            payload = json.dumps(
                {
                    "q": batch,
                    "scopes": "ensembl.gene",
                    "fields": "symbol",
                    "species": "human",
                }
            ).encode()
            req = urllib.request.Request(
                url, data=payload, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=120) as resp:
                result = json.loads(resp.read().decode())
            for hit in result:
                ensg = hit.get("query")
                sym = hit.get("symbol")
                if ensg and sym:
                    ensembl_map[ensg] = sym
            if (i // batch_size) % 10 == 0:
                print(f"    Batch {i//batch_size}: {len(ensembl_map)} mapped so far")
        except Exception as e:
            print(f"    Batch {i//batch_size}: error - {e}")
            continue

    print(f"  Total mapped: {len(ensembl_map)} / {len(all_ensembl)}")

    with open(f"{ANNOT_DIR}/ensembl_to_symbol.tsv", "w") as f:
        f.write("ensembl_id\tgene_symbol\n")
        for e, s in sorted(ensembl_map.items()):
            f.write(f"{e}\t{s}\n")
    return ensembl_map


def load_gpl570_mapping():
    """Load previously created GPL570 mapping."""
    path = f"{ANNOT_DIR}/gpl570_probe_to_symbol.tsv"
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t")
        return dict(zip(df["probe_id"], df["gene_symbol"]))
    return {}


def main():
    print("=" * 60)
    print("BUILDING GENE ID MAPPINGS")
    print("=" * 60)

    # GPL570 (already done by GEOparse)
    gpl570 = load_gpl570_mapping()
    print(f"\nGPL570 (Affymetrix): {len(gpl570)} probes already mapped")

    # GPL10558
    gpl10558 = build_illumina_mapping()

    # Ensembl
    ensembl = build_ensembl_mapping()

    print("\n" + "=" * 60)
    print("MAPPING SUMMARY")
    print("=" * 60)
    print(f"  GPL570 (Affymetrix):   {len(gpl570):,} probes")
    print(f"  GPL10558 (Illumina):   {len(gpl10558):,} probes")
    print(f"  Ensembl → Symbol:      {len(ensembl):,} genes")


if __name__ == "__main__":
    main()
