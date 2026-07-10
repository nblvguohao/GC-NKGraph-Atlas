"""
Remap all expression matrices from probe/Ensembl IDs to HGNC gene symbols.

After this, all expression files will use gene symbols matching the SST modules.
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

os.environ["no_proxy"] = "*"
os.environ["all_proxy"] = ""
os.environ["ALL_PROXY"] = ""

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir

ANNOT_DIR = "data/raw/annotations"

# SST checkpoint genes to verify
SST_GENES = [
    "SGMS1", "SGMS2", "SMPD1", "SMPD2", "SMPD3", "SMPD4",
    "PHGDH", "PSAT1", "PSPH",
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG",
    "HAVCR2", "TIGIT", "CD96",
    "EZR", "MSN", "RDX", "RAC1", "CDC42", "RHOA",
    "NCAM1", "NCR1", "KLRD1", "KLRK1", "FCGR3A",
]


def load_ensembl_mapping():
    """Load Ensembl→Symbol mapping."""
    path = f"{ANNOT_DIR}/ensembl_to_symbol.tsv"
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["ensembl_id"], df["gene_symbol"]))


def load_gpl570_mapping():
    """Load GPL570 probe→Symbol mapping."""
    path = f"{ANNOT_DIR}/gpl570_probe_to_symbol.tsv"
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["probe_id"], df["gene_symbol"]))


def load_illumina_mapping():
    """Build Illumina mapping from GPL10558 annotation."""
    import gzip
    path = f"{ANNOT_DIR}/GPL10558.annot.gz"
    if not os.path.exists(path):
        return {}

    symbol_map = {}
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        # Find the data table (starts after the metadata header)
        data_start = 0
        header_line = ""
        for i, line in enumerate(lines):
            if line.startswith("ID\t") or line.startswith("ID,"):
                data_start = i
                header_line = line.strip()
                break

        if data_start == 0:
            print(f"    WARNING: Could not find data table header")
            return {}

        # Parse header
        cols = header_line.split("\t")
        id_col = 0
        symbol_col = None
        for j, c in enumerate(cols):
            if c.lower() in ("symbol", "gene symbol", "gene_symbol"):
                symbol_col = j
                break

        if symbol_col is None:
            print(f"    WARNING: No symbol column found. Cols: {cols[:10]}")
            return {}

        # Parse data
        from io import StringIO
        data_text = "".join(lines[data_start:])
        # Find end of table
        end = data_text.find("\n!")
        if end > 0:
            data_text = data_text[:end]

        df = pd.read_csv(StringIO(data_text), sep="\t", low_memory=False)
        for _, row in df.iterrows():
            probe = str(row.iloc[id_col]).strip()
            symbol = str(row.iloc[symbol_col]).strip()
            if symbol and symbol != "nan" and symbol != "":
                symbol = symbol.split("///")[0].strip()
                if symbol:
                    symbol_map[probe] = symbol

    except Exception as e:
        print(f"    ERROR: {e}")
        return {}

    return symbol_map


def remap_matrix(expr_path, mapping, output_path, id_type="ensembl"):
    """Remap expression matrix columns to gene symbols.

    Collapses duplicates by keeping the row with highest mean expression.
    """
    print(f"  Loading: {expr_path} ({os.path.getsize(expr_path)/1e6:.0f} MB)")

    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    n_orig = df.shape[1]

    # Build mapping from old col names to symbols
    col_to_symbol = {}
    for col in df.columns:
        if id_type == "ensembl":
            clean = col.split(".")[0]
            sym = mapping.get(clean)
        else:
            sym = mapping.get(col)
        if sym:
            col_to_symbol[col] = sym

    # Rename and keep only mapped columns
    df_mapped = df[list(col_to_symbol.keys())].rename(columns=col_to_symbol)

    # Collapse duplicate gene symbols: keep the probe with highest mean expression
    if df_mapped.columns.duplicated().any():
        dup_count = df_mapped.columns.duplicated().sum()
        print(f"    Collapsing {dup_count} duplicate genes...")
        col_to_keep = {}
        for gene in df_mapped.columns.unique():
            dup_indices = [i for i, c in enumerate(df_mapped.columns) if c == gene]
            if len(dup_indices) > 1:
                # Keep the one with highest mean
                means = [df_mapped.iloc[:, i].mean() for i in dup_indices]
                best_idx = dup_indices[means.index(max(means))]
                col_to_keep[gene] = best_idx
            else:
                col_to_keep[gene] = dup_indices[0]
        # Build result using position indices
        result = {}
        for gene, idx in col_to_keep.items():
            result[gene] = df_mapped.iloc[:, idx].values
        df_final = pd.DataFrame(result, index=df_mapped.index)
    else:
        df_final = df_mapped

    df_final.index.name = "sample_id"

    # Check SST gene coverage
    found = [g for g in SST_GENES if g in df_final.columns]
    missing = [g for g in SST_GENES if g not in df_final.columns]
    pct_mapped = 100 * len(col_to_symbol) / n_orig

    df_final.to_csv(output_path, sep="\t")
    size_out = os.path.getsize(output_path) / 1e6
    print(f"    → {df_final.shape[0]} x {df_final.shape[1]} ({pct_mapped:.0f}% mapped, {size_out:.0f} MB)")
    print(f"    SST genes: {len(found)}/{len(SST_GENES)} found")

    if missing and len(missing) < 10:
        print(f"    Missing: {missing}")

    return df_final


def main():
    print("=" * 60)
    print("REMAP ALL EXPRESSION MATRICES TO GENE SYMBOLS")
    print("=" * 60)

    # Load mappings
    ensembl_map = load_ensembl_mapping()
    gpl570_map = load_gpl570_mapping()
    illumina_map = load_illumina_mapping()

    print(f"\nMappings loaded:")
    print(f"  Ensembl:  {len(ensembl_map):,}")
    print(f"  GPL570:   {len(gpl570_map):,}")
    print(f"  Illumina: {len(illumina_map):,}")

    # Remap each dataset
    datasets = [
        ("TCGA-STAD", "data/processed/bulk/tcga_stad_expression.tsv",
         ensembl_map, "ensembl",
         "data/processed/bulk/tcga_stad_expression.tsv"),
        ("TCGA-LIHC", "data/processed/bulk/tcga_lihc_expression.tsv",
         ensembl_map, "ensembl",
         "data/processed/bulk/tcga_lihc_expression.tsv"),
        ("GSE62254", "data/processed/bulk/gse62254_expression.tsv",
         gpl570_map, "probe",
         "data/processed/bulk/gse62254_expression.tsv"),
        ("GSE84437", "data/processed/bulk/gse84437_expression.tsv",
         illumina_map, "probe",
         "data/processed/bulk/gse84437_expression.tsv"),
    ]

    for name, path, mapping, id_type, output in datasets:
        print(f"\n--- {name} ---")
        if not os.path.exists(path):
            print(f"  SKIP: file not found")
            continue
        if not mapping:
            print(f"  SKIP: no mapping available")
            continue
        remap_matrix(path, mapping, output, id_type)

    print(f"\n{'=' * 60}")
    print("REMAP COMPLETE")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
