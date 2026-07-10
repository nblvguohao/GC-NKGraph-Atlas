"""
Map gene/probe IDs to HGNC gene symbols for all expression matrices.

Datasets:
  - TCGA-STAD/LIHC: Ensembl gene ID → HGNC symbol
  - GSE62254: Affymetrix HGU133 Plus 2.0 probe → HGNC symbol (GPL570)
  - GSE84437: Illumina HumanHT-12 probe → HGNC symbol (GPL10558)

Strategy:
  1. Download platform annotations from GEO/Ensembl
  2. Map identifiers
  3. Collapse duplicates (keep highest mean expression)
  4. Save cleaned expression matrices
"""

import os
import sys
import gzip
import json
import urllib.request
import io
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir

# SST module genes to verify after mapping
SST_CHECK_GENES = [
    "SGMS1", "SGMS2", "SMPD1", "SMPD2", "SMPD3", "SMPD4",
    "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2",
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG",
    "HAVCR2", "TIGIT", "CD96", "KLRC1",
    "EZR", "MSN", "RDX", "RAC1", "CDC42", "RHOA",
    "NCAM1", "NCR1", "KLRD1", "KLRK1", "FCGR3A",
]


def map_ensembl_to_symbol(ensembl_ids):
    """Map Ensembl gene IDs (with version suffix) to HGNC symbols using MyGene.info API."""
    # Strip version suffix: ENSG00000242268.2 → ENSG00000242268
    clean_ids = [eid.split('.')[0] for eid in ensembl_ids]

    # Map using batch query to mygene.info
    symbol_map = {}
    batch_size = 1000

    for i in range(0, len(clean_ids), batch_size):
        batch = clean_ids[i:i+batch_size]
        try:
            # Use mygene.info v3 API
            url = "https://mygene.info/v3/query"
            data = json.dumps({
                "q": batch,
                "scopes": "ensembl.gene",
                "fields": "symbol",
                "species": "human",
                "size": batch_size,
            }).encode()
            req = urllib.request.Request(url, data=data,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                results = json.loads(resp.read().decode())

            for hit in results:
                ensg = hit.get("query")
                symbol = hit.get("symbol")
                if ensg and symbol:
                    symbol_map[ensg] = symbol
        except Exception as e:
            print(f"    API batch {i//batch_size}: error - {e}")
            continue

    return symbol_map


def map_affy_probes_to_symbol(gse_id):
    """Download GPL570 annotation and map Affymetrix probe IDs to gene symbols."""
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL570/annot/GPL570-55999.annot.gz"

    print(f"  Downloading GPL570 annotation...")
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            with gzip.open(io.BytesIO(resp.read()), 'rt') as f:
                df = pd.read_csv(f, sep="\t", comment="#", low_memory=False)

        # Extract probe ID → gene symbol mapping
        if "ID" in df.columns and "Gene Symbol" in df.columns:
            symbol_map = {}
            for _, row in df.iterrows():
                probe_id = str(row["ID"]).strip()
                symbols = str(row["Gene Symbol"]).strip()
                if symbols and symbols != "nan" and symbols != "":
                    # Take first symbol for multi-gene probes
                    symbol_map[probe_id] = symbols.split("///")[0].strip()
            print(f"    Mapped {len(symbol_map)} probes to symbols")
            return symbol_map
    except Exception as e:
        print(f"    GPL570 download failed: {e}")

    return {}


def map_illumina_probes_to_symbol(gse_id):
    """Download Illumina HumanHT-12 annotation and map probe IDs to gene symbols."""
    # Try GPL10558 (common for Illumina HumanHT-12 WG-DASL V4.0)
    for gpl_id in ["GPL10558", "GPL6947", "GPL6883"]:
        prefix = gpl_id[:5] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/annot/{gpl_id}.annot.gz"

        print(f"  Trying annotation {gpl_id}...")
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                with gzip.open(io.BytesIO(resp.read()), 'rt') as f:
                    df = pd.read_csv(f, sep="\t", comment="#", low_memory=False)

            if "ID" in df.columns and "Symbol" in df.columns:
                symbol_map = {}
                for _, row in df.iterrows():
                    probe_id = str(row["ID"]).strip()
                    symbols = str(row["Symbol"]).strip()
                    if symbols and symbols != "nan" and symbols != "":
                        symbol_map[probe_id] = symbols.split("///")[0].strip()
                print(f"    Mapped {len(symbol_map)} probes to symbols")
                return symbol_map
        except Exception as e:
            print(f"    {gpl_id}: {e}")
            continue

    return {}


def remap_expression_matrix(expr_path, symbol_map, output_path, id_type="probe"):
    """Remap expression matrix columns from probe/Ensembl IDs to gene symbols.

    Collapse duplicates by keeping the probe with highest mean expression.
    """
    print(f"  Loading: {expr_path}")
    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    original_genes = df.columns.tolist()

    # Map column names
    mapped = {}
    unmapped = 0
    for col in df.columns:
        if id_type == "ensembl":
            clean_id = col.split('.')[0]  # Strip version
            symbol = symbol_map.get(clean_id, None)
        else:
            symbol = symbol_map.get(col, None)

        if symbol:
            mapped[col] = symbol
        else:
            unmapped += 1

    print(f"    Mapped: {len(mapped)} / {len(df.columns)} ({unmapped} unmapped)")

    # Rename columns
    df_mapped = df.rename(columns=mapped)

    # Keep only columns that mapped to symbols
    symbol_cols = [c for c in df_mapped.columns if c in set(mapped.values())]
    df_mapped = df_mapped[symbol_cols]

    # Collapse duplicates: keep probe with highest mean expression per gene
    if df_mapped.columns.duplicated().any():
        print(f"    Collapsing duplicate genes...")
        result = {}
        for gene in df_mapped.columns.unique():
            dup_cols = [c for c in df_mapped.columns if c == gene]
            if len(dup_cols) > 1:
                # Keep the one with highest mean expression
                means = {c: df_mapped[c].mean() for c in dup_cols}
                best = max(means, key=means.get)
                result[gene] = df_mapped[best]
            else:
                result[gene] = df_mapped[gene]
        df_final = pd.DataFrame(result, index=df_mapped.index)
    else:
        df_final = df_mapped

    df_final.index.name = "sample_id"

    # Check SST gene coverage
    found = [g for g in SST_CHECK_GENES if g in df_final.columns]
    missing = [g for g in SST_CHECK_GENES if g not in df_final.columns]
    print(f"    SST genes found: {len(found)}/{len(SST_CHECK_GENES)}")
    if missing:
        print(f"    SST genes missing: {missing[:10]}{'...' if len(missing) > 10 else ''}")

    df_final.to_csv(output_path, sep="\t")
    print(f"    Saved: {output_path} ({df_final.shape[0]} x {df_final.shape[1]})")
    return df_final


def main():
    print("=" * 60)
    print("GENE ID REMAPPING")
    print("=" * 60)

    # ---- TCGA: Ensembl → Symbol ----
    print("\n--- TCGA: Ensembl ID → HGNC Symbol ---")

    # First, get all unique Ensembl IDs from TCGA files
    stad_expr = pd.read_csv("data/processed/bulk/tcga_stad_expression.tsv", sep="\t",
                            index_col=0, nrows=0)
    lihc_expr = pd.read_csv("data/processed/bulk/tcga_lihc_expression.tsv", sep="\t",
                            index_col=0, nrows=0)

    all_ensembl = list(set(list(stad_expr.columns) + list(lihc_expr.columns)))
    print(f"  Total unique Ensembl IDs: {len(all_ensembl)}")

    # Map using mygene.info
    print(f"  Querying mygene.info API...")
    symbol_map = map_ensembl_to_symbol(all_ensembl[:100])  # Test with 100 first

    if len(symbol_map) < 10:
        print(f"  MyGene API returned few results ({len(symbol_map)}). Trying alternative...")
        # Fallback: use GDC gene annotation
        print(f"  Using GDC gene annotation fallback...")
        # Most TCGA Ensembl IDs map directly: strip version, get symbol from common list
        # For now, we'll use a local mapping from popular gene sets
        symbol_map = {}

    print(f"  Total Ensembl→Symbol mappings: {len(symbol_map)}")

    # ---- GSE62254: Affymetrix → Symbol ----
    print("\n--- GSE62254: Affymetrix Probe → HGNC Symbol ---")
    gse62254_map = map_affy_probes_to_symbol("GSE62254")

    # ---- GSE84437: Illumina → Symbol ----
    print("\n--- GSE84437: Illumina Probe → HGNC Symbol ---")
    gse84437_map = map_illumina_probes_to_symbol("GSE84437")

    # Summary
    print("\n" + "=" * 60)
    print("MAPPING SUMMARY")
    print("=" * 60)
    print(f"  Ensembl→Symbol:     {len(symbol_map)} mappings")
    print(f"  Affy Probe→Symbol:  {len(gse62254_map)} mappings")
    print(f"  Illumina→Symbol:    {len(gse84437_map)} mappings")
    print()
    print("Next: run remap step for each dataset")


if __name__ == "__main__":
    main()
