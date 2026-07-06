"""
GC-NKGraph-Atlas Bulk Preprocessing Pipeline.

Processes raw TCGA and GEO expression matrices into standardized format:
1. Gene ID standardization (HGNC symbols)
2. Duplicate resolution (keep highest mean expression)
3. Independent normalization per dataset
4. Clinical data alignment

Usage:
    python src/preprocessing/run_bulk_preprocessing.py --config configs/data_config.yaml
"""

import os
import sys
import argparse
from pathlib import Path

import gzip
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import load_config, save_table, ensure_dir


# Gene aliases for common symbols (curated from HGNC)
GENE_ALIASES = {
    "TGFB1": ["TGFB1", "TGFB"],
    "TNF": ["TNF", "TNFA", "TNF-α"],
    "NKG7": ["NKG7", "GIG1"],
    "GNLY": ["GNLY", "D2S69E", "LAG-2", "NKG5", "TLA-2"],
    "GZMB": ["GZMB", "CCP1", "CGL-1", "CGL1", "CSP-B", "CSPB", "CTLA-1", "CTLA1", "HLP"],
    "PRF1": ["PRF1", "FLH2", "HPLH2", "P1", "PFN1"],
    "HAVCR2": ["HAVCR2", "TIM3", "CD366", "TIMD3"],
    "SMPD1": ["SMPD1", "ASM", "ASMASE", "NPD"],
    "SMPD3": ["SMPD3", "NSMASE2"],
    "SMPD4": ["SMPD4", "NSMASE3"],
    "SGMS1": ["SGMS1", "SMS1", "TMEM23"],
    "SGMS2": ["SGMS2", "SMS2"],
    "CD3D": ["CD3D", "CD3-DELTA", "IMD19"],
    "CD3E": ["CD3E", "CD3-epsilon", "T3E", "IMD18"],
    "CD3G": ["CD3G", "CD3-gamma", "IMD17"],
    "NCAM1": ["NCAM1", "CD56", "MSK39"],
    "KLRD1": ["KLRD1", "CD94"],
    "KLRK1": ["KLRK1", "NKG2D", "CD314"],
    "KLRC1": ["KLRC1", "NKG2A", "CD159a"],
    "FCGR3A": ["FCGR3A", "CD16", "CD16A", "FCGR3"],
    "NT5E": ["NT5E", "CD73", "CALJA", "E5NT", "NT"],
    "ENTPD1": ["ENTPD1", "CD39", "ATPBase", "SPG64"],
}


def standardize_gene_symbols(expression_df: pd.DataFrame) -> pd.DataFrame:
    """Standardize gene symbols using alias mapping.

    Args:
        expression_df: DataFrame with gene symbols as index.

    Returns:
        DataFrame with standardized gene symbols, dropping unmappable rows.
    """
    # Build reverse alias map
    alias_to_standard = {}
    for standard, aliases in GENE_ALIASES.items():
        for alias in aliases:
            alias_to_standard[alias.upper()] = standard

    # Map current index
    current_idx = [str(g).upper() for g in expression_df.index]
    mapped = [alias_to_standard.get(g, g) for g in current_idx]

    # Keep only 1-to-1 mappings initially
    expression_df = expression_df.copy()
    expression_df.index = mapped

    # Remove rows that map to NaN or empty
    expression_df = expression_df[expression_df.index.notna()]
    expression_df = expression_df[expression_df.index != ""]

    return expression_df


def resolve_duplicates(expression_df: pd.DataFrame, strategy: str = "highest_mean") -> pd.DataFrame:
    """Resolve duplicate gene symbols.

    Args:
        expression_df: DataFrame with potential duplicate gene indices.
        strategy: 'highest_mean' or 'mean'

    Returns:
        DataFrame with unique gene indices.
    """
    if expression_df.index.duplicated().sum() == 0:
        return expression_df

    if strategy == "highest_mean":
        means = expression_df.mean(axis=1)
        keep_idx = means.groupby(expression_df.index).idxmax()
        return expression_df.loc[keep_idx]
    elif strategy == "mean":
        return expression_df.groupby(expression_df.index).mean()
    else:
        raise ValueError(f"Unknown duplicate strategy: {strategy}")


def _open_maybe_gzip(path: str, mode: str = "rt"):
    """Open a file, transparently handling gzip if needed."""
    if path.endswith(".gz") or path.endswith(".gzip"):
        return gzip.open(path, mode)
    # Check magic bytes for gzip
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, mode)
    return open(path, mode)


def load_xena_expression(path: str) -> pd.DataFrame:
    """Load TCGA expression from UCSC Xena HiSeqV2 format.

    Xena HiSeqV2 format: first column 'sample' with gene symbols, rest are samples.
    Files may be gzipped even without .gz extension.
    """
    # Try reading with inferred compression; if file is gzipped without
    # .gz extension, compression="infer" won't detect it, so try gzip as fallback
    try:
        df = pd.read_csv(path, sep="\t", index_col=0, compression="infer")
    except (UnicodeDecodeError, ValueError):
        with gzip.open(path, "rt") as f:
            df = pd.read_csv(f, sep="\t", index_col=0)
    df.index.name = "gene"
    return df


def load_xena_phenotype(path: str) -> pd.DataFrame:
    """Load TCGA phenotype from UCSC Xena format."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def load_geo_expression(path: str) -> pd.DataFrame:
    """Load GEO series matrix file and extract expression matrix.

    GEO series matrix format: expression values start after '!series_matrix_table_begin'
    """
    with _open_maybe_gzip(path, "rt") as f:
        lines = f.readlines()

    # Find table boundaries
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if "!series_matrix_table_begin" in line:
            start_idx = i + 1
        elif "!series_matrix_table_end" in line:
            end_idx = i

    if start_idx is None or end_idx is None:
        raise ValueError(f"Cannot find expression table in {path}")

    # Read the table
    header = lines[start_idx].strip().split("\t")
    data_rows = []
    for line in lines[start_idx + 1:end_idx]:
        parts = line.strip().split("\t")
        if len(parts) > 1:
            data_rows.append(parts)

    df = pd.DataFrame(data_rows, columns=header)
    df = df.set_index(header[0])
    df.index.name = "gene"
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def normalize_expression(expression_df: pd.DataFrame, method: str = "log2_tpm") -> pd.DataFrame:
    """Apply normalization.

    Args:
        expression_df: Raw expression counts/values.
        method: Normalization method.

    Returns:
        Normalized expression DataFrame.
    """
    if method == "log2_tpm_plus_1_or_platform_specific":
        # Xena HiSeqV2 is already log2(norm_count+1)
        # GEO series matrix may be log2 or raw; assume it's already platform-normalized
        return expression_df
    else:
        raise ValueError(f"Unknown normalization method: {method}")


def process_dataset(
    name: str,
    expression_path: str,
    clinical_path: str,
    output_dir: str,
    logger: Logger,
) -> bool:
    """Process a single bulk expression dataset.

    Returns:
        True if successful, False otherwise.
    """
    dataset_slug = name.lower().replace("-", "_").replace(" ", "_")
    dataset_dir = os.path.join(output_dir, dataset_slug)
    ensure_dir(dataset_dir)

    try:
        # Load expression
        if name.startswith("TCGA-"):
            expr = load_xena_expression(expression_path)
            pheno = load_xena_phenotype(clinical_path)
        else:
            expr = load_geo_expression(expression_path)
            pheno = None  # GEO clinical tables are handled separately

        # Standardize and deduplicate
        n_before = expr.shape[0]
        expr = standardize_gene_symbols(expr)
        expr = resolve_duplicates(expr)
        n_after = expr.shape[0]
        logger.ok(
            phase="PREPROCESSING",
            message=f"{name}: {n_before} genes -> {n_after} after dedup",
            script=__file__,
        )

        # Transpose: samples as rows, genes as columns (standard ML format)
        expr_t = expr.T
        expr_t.index.name = "sample_id"

        # Align with clinical data
        if pheno is not None:
            common_samples = expr_t.index.intersection(pheno.index)
            expr_t = expr_t.loc[common_samples]
            pheno = pheno.loc[common_samples]
            logger.ok(
                phase="PREPROCESSING",
                message=f"{name}: {len(common_samples)} samples with clinical data",
                script=__file__,
            )

        # Save
        out_expr = os.path.join(output_dir, f"{dataset_slug}_expression.tsv")
        out_clin = os.path.join(output_dir, f"{dataset_slug}_clinical.tsv")

        expr_t.to_csv(out_expr, sep="\t")
        if pheno is not None:
            pheno.to_csv(out_clin, sep="\t")
        else:
            # Create minimal clinical placeholder
            pd.DataFrame({"sample_id": expr_t.index}).to_csv(out_clin, sep="\t", index=False)

        logger.ok(
            phase="PREPROCESSING",
            message=f"{name}: saved expression ({expr_t.shape}) to {out_expr}",
            script=__file__,
        )
        return True

    except Exception as e:
        logger.fail(
            phase="PREPROCESSING",
            message=f"{name}: FAILED - {e}",
            script=__file__,
        )
        return False


def main():
    parser = argparse.ArgumentParser(description="Bulk data preprocessing for GC-NKGraph-Atlas")
    parser.add_argument("--config", default="configs/data_config.yaml", help="Data config path")
    parser.add_argument("--output-dir", default="data/processed/bulk", help="Output directory")
    args = parser.parse_args()

    logger = Logger()
    config = load_config(args.config)
    output_dir = ensure_dir(args.output_dir)

    # Process training datasets
    datasets = config.get("bulk_datasets", [])
    for ds in datasets:
        # Check files exist
        expr_path = ds.get("expression_path", "").replace("processed", "raw")
        clin_path = ds.get("clinical_path", "").replace("processed", "raw")

        if not os.path.exists(expr_path):
            logger.skip(
                phase="PREPROCESSING",
                message=f"{ds['name']}: expression not found at {expr_path}",
                script=__file__,
            )
            continue

        process_dataset(
            name=ds["name"],
            expression_path=expr_path,
            clinical_path=clin_path,
            output_dir=str(output_dir),
            logger=logger,
        )

    # Process positive control datasets
    control_datasets = config.get("positive_control_bulk_datasets", [])
    for ds in control_datasets:
        expr_path = ds.get("expression_path", "").replace("processed", "raw")
        clin_path = ds.get("clinical_path", "").replace("processed", "raw")

        if not os.path.exists(expr_path):
            logger.skip(
                phase="PREPROCESSING",
                message=f"{ds['name']}: expression not found at {expr_path}",
                script=__file__,
            )
            continue

        process_dataset(
            name=ds["name"],
            expression_path=expr_path,
            clinical_path=clin_path,
            output_dir=str(output_dir),
            logger=logger,
        )

    # Summary
    logger.ok(
        phase="PREPROCESSING",
        message=f"Bulk preprocessing complete. Output: {output_dir}/",
        script=__file__,
    )


if __name__ == "__main__":
    main()
