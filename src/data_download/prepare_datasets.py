"""
Prepare downloaded datasets for GC-NKGraph-Atlas pipeline.

1. Extract TCGA-STAD and TCGA-LIHC expression from pan-cancer RSEM TPM
2. Extract GEO series matrices to expression format
3. Extract GSE246662 scRNA RAW data

Output: data/processed/bulk/ expression matrices ready for the pipeline.
"""

import os
import sys
import gzip
import tarfile
import shutil
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir


def prepare_tcga():
    """Extract TCGA-STAD and TCGA-LIHC from pan-cancer RSEM TPM."""
    print("\n" + "=" * 60)
    print("PREPARING TCGA EXPRESSION MATRICES")
    print("=" * 60)

    rsem_path = "data/raw/bulk/tcga/tcga_RSEM_gene_tpm.gz"
    if not os.path.exists(rsem_path):
        print("ERROR: RSEM file not found!")
        return False

    # Load sample lists
    stad_samples = set()
    lihc_samples = set()
    for fname, sample_set in [
        ("data/raw/bulk/tcga/stad_samples.txt", stad_samples),
        ("data/raw/bulk/tcga/lihc_samples.txt", lihc_samples),
    ]:
        if os.path.exists(fname):
            with open(fname) as f:
                for line in f:
                    sample_set.add(line.strip())

    print(f"TCGA-STAD samples: {len(stad_samples)}")
    print(f"TCGA-LIHC samples: {len(lihc_samples)}")

    # Load RSEM data
    print("Loading pan-cancer RSEM TPM (this may take a minute)...")
    df = pd.read_csv(rsem_path, sep="\t", compression="gzip", index_col=0)
    print(f"Full dimensions: {df.shape[0]} genes x {df.shape[1]} samples")

    # Extract STAD samples - match by stripping the vial suffix (e.g. -01A -> -01)
    # Xena RSEM uses TCGA-XX-YYYY-ZZ format, GDC uses TCGA-XX-YYYY-ZZV
    def match_xena_id(gdc_id):
        """Convert GDC submitter_id to Xena sample ID format."""
        # GDC: TCGA-XX-YYYY-ZZV (V=vial: A/B/C)
        # Xena: TCGA-XX-YYYY-ZZ
        parts = gdc_id.split('-')
        if len(parts) >= 4 and len(parts[3]) > 2:
            parts[3] = parts[3][:2]  # Strip vial suffix
        return '-'.join(parts)

    stad_xena_ids = {match_xena_id(s) for s in stad_samples}
    lihc_xena_ids = {match_xena_id(s) for s in lihc_samples}
    stad_cols = [c for c in df.columns if c in stad_xena_ids]
    lihc_cols = [c for c in df.columns if c in lihc_xena_ids]

    # Save expression matrices
    ensure_dir("data/processed/bulk")

    if stad_cols:
        stad_df = df[stad_cols]
        # Transpose: genes x samples → samples x genes (project convention)
        stad_df = stad_df.T
        stad_df.index.name = "sample_id"
        output = "data/processed/bulk/tcga_stad_expression.tsv"
        stad_df.to_csv(output, sep="\t")
        print(f"\nSaved: {output} ({stad_df.shape[0]} samples x {stad_df.shape[1]} genes)")

    if lihc_cols:
        lihc_df = df[lihc_cols]
        lihc_df = lihc_df.T
        lihc_df.index.name = "sample_id"
        output = "data/processed/bulk/tcga_lihc_expression.tsv"
        lihc_df.to_csv(output, sep="\t")
        print(f"Saved: {output} ({lihc_df.shape[0]} samples x {lihc_df.shape[1]} genes)")

    # Save clinical stub (project convention expects clinical files)
    for name, samples in [("tcga_stad", stad_cols), ("tcga_lihc", lihc_cols)]:
        if samples:
            clinical = pd.DataFrame({"sample_id": list(samples)})
            clinical.set_index("sample_id", inplace=True)
            clinical.to_csv(f"data/processed/bulk/{name}_clinical.tsv", sep="\t")
            print(f"Saved: data/processed/bulk/{name}_clinical.tsv")

    return True


def prepare_geo_series(gse_id):
    """Extract GEO series matrix to expression format."""
    print(f"\n--- Preparing {gse_id} ---")

    matrix_path = f"data/raw/bulk/{gse_id.lower()}/{gse_id}_series_matrix.txt.gz"
    if not os.path.exists(matrix_path):
        print(f"  WARNING: Matrix file not found at {matrix_path}")
        return False

    try:
        # GEO series matrices have metadata lines starting with !
        # and tab-separated expression data
        with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        # Parse metadata and expression sections
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith("!series_matrix_table_begin"):
                data_start = i + 1
                break

        if data_start == 0:
            print(f"  ERROR: Could not find data start marker")
            return False

        # Read expression data
        from io import StringIO
        data_text = "".join(lines[data_start:])
        end = data_text.find("!series_matrix_table_end")
        if end > 0:
            data_text = data_text[:end]

        df = pd.read_csv(StringIO(data_text), sep="\t", index_col=0)
        # Transpose to samples x genes
        df = df.T
        df.index.name = "sample_id"

        output = f"data/processed/bulk/{gse_id.lower()}_expression.tsv"
        ensure_dir("data/processed/bulk")
        df.to_csv(output, sep="\t")
        print(f"  Saved: {output} ({df.shape[0]} samples x {df.shape[1]} probes)")

        # Clinical stub
        clinical = pd.DataFrame({"sample_id": df.index})
        clinical.set_index("sample_id", inplace=True)
        clinical.to_csv(f"data/processed/bulk/{gse_id.lower()}_clinical.tsv", sep="\t")

        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def extract_scnra_raw():
    """Extract GSE246662 RAW tar file."""
    print(f"\n--- Extracting GSE246662 RAW data ---")

    tar_path = "data/raw/scrna/gse246662/GSE246662_RAW.tar"
    extract_dir = "data/raw/scrna/gse246662/extracted"

    if not os.path.exists(tar_path):
        print(f"  WARNING: RAW tar not found at {tar_path}")
        return False

    # Check if already extracted
    existing = list(Path(extract_dir).glob("*.gz"))
    existing += list(Path(extract_dir).glob("*.mtx"))
    existing += list(Path(extract_dir).glob("*.h5"))
    if existing:
        print(f"  Already extracted: {len(existing)} files in {extract_dir}")
        return True

    ensure_dir(extract_dir)
    print(f"  Extracting {os.path.getsize(tar_path)/1e6:.0f} MB tar to {extract_dir}...")

    try:
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=extract_dir)
        files = list(Path(extract_dir).glob("*"))
        print(f"  Extracted {len(files)} files")
        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def main():
    print("=" * 60)
    print("GC-NKGraph-Atlas Data Preparation")
    print("=" * 60)

    ok = True

    # TCGA
    if not prepare_tcga():
        ok = False

    # GEO bulk
    for gse_id in ["GSE62254", "GSE84437"]:
        if not prepare_geo_series(gse_id):
            ok = False

    # scRNA
    if not extract_scnra_raw():
        ok = False

    # Summary
    print("\n" + "=" * 60)
    print("DATA PREPARATION SUMMARY")
    print("=" * 60)

    expected_files = [
        "data/processed/bulk/tcga_stad_expression.tsv",
        "data/processed/bulk/tcga_stad_clinical.tsv",
        "data/processed/bulk/tcga_lihc_expression.tsv",
        "data/processed/bulk/tcga_lihc_clinical.tsv",
        "data/processed/bulk/gse62254_expression.tsv",
        "data/processed/bulk/gse62254_clinical.tsv",
        "data/processed/bulk/gse84437_expression.tsv",
        "data/processed/bulk/gse84437_clinical.tsv",
    ]

    for f in expected_files:
        if os.path.exists(f):
            size = os.path.getsize(f)
            print(f"  [OK] {f} ({size/1e3:.0f} KB)")
        else:
            print(f"  [MISSING] {f}")

    print(f"\n{'ALL DATA READY' if ok else 'SOME DATA MISSING - CHECK ABOVE'}")


if __name__ == "__main__":
    main()
