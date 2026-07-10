"""
GC-NKGraph-Atlas Synthetic Data Generator.

Generates realistic synthetic data for end-to-end pipeline testing when
real TCGA/GEO/scRNA data is not yet downloaded. All generated data is
deterministic (fixed seed) and flagged with `_SYNTHETIC_` markers.

Usage:
    python src/common/synthetic_data.py --output-dir data/synthetic

This creates a complete minimal dataset that can drive every phase of the
pipeline — from preprocessing through model training to target prioritization.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Project imports
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# =========================================================================
# Seed candidates used across multiple generators
# =========================================================================

SST_GENES = [
    "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1", "MTHFD2",
    "SGMS1", "SGMS2", "SMPD1", "SMPD2", "SMPD3", "SMPD4",
    "EZR", "MSN", "RDX", "ACTR2", "ACTR3", "CDC42", "RAC1", "RHOA",
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "HAVCR2",
]

CANDIDATE_GENES = [
    "LDHA", "SLC16A3", "NT5E", "ENTPD1", "ADORA2A",
    "PVR", "NECTIN2", "HLA-E", "MICA", "MICB", "ADAM10", "ADAM17",
    "TGFB1", "TGFBR1", "COL1A1", "FN1", "CLDN18", "ERBB2", "FGFR2", "MET",
]

ALL_GENES = sorted(set(SST_GENES + CANDIDATE_GENES + [
    "GAPDH", "ACTB", "CD8A", "CD4", "FOXP3", "CD68", "CD163",
    "EPCAM", "KRT19", "KRT18", "CDH1", "VIM", "CD3D", "CD3E",
    "NCAM1", "KLRD1", "KLRK1", "KLRC1", "FCGR3A", "TIGIT", "CD96",
    "TOX", "XCL1", "XCL2", "CCL5", "LCP2", "LAT", "VAV1", "TLN1",
]))


# =========================================================================
# Generator functions
# =========================================================================

def generate_bulk_expression(
    n_samples: int = 300,
    n_genes: int = 100,
    seed: int = 42,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Generate synthetic bulk RNA-seq expression matrix.

    Simulates log2(TPM+1) values with realistic gene-gene correlations.

    Args:
        n_samples: Number of tumor samples.
        n_genes: Number of genes (must <= len(ALL_GENES)).
        seed: Random seed for reproducibility.
        output_dir: If set, write TSV files to this directory.

    Returns:
        Expression DataFrame (samples × genes).
    """
    rng = np.random.RandomState(seed)
    n_genes = min(n_genes, len(ALL_GENES))
    genes = ALL_GENES[:n_genes]

    # Base expression: log-normal
    base = rng.lognormal(mean=2.0, sigma=1.0, size=(n_samples, n_genes))

    # Add structure: some genes co-vary via latent factors
    n_factors = 5
    factors = rng.randn(n_samples, n_factors)
    loadings = rng.randn(n_factors, n_genes) * 0.5
    expr = base + factors @ loadings

    # Clip to realistic log2(TPM+1) range
    expr = np.clip(expr, 0, 15)

    df = pd.DataFrame(expr, columns=genes, dtype=np.float32)
    df.index = [f"TCGA-SYNTH-{i:04d}" for i in range(n_samples)]
    df.index.name = "sample_id"

    if output_dir:
        path = os.path.join(output_dir, "tcga_stad_expression_synthetic.tsv")
        df.to_csv(path, sep="\t")
        log(f"  Saved: {path}")

    return df


def generate_clinical_data(
    n_samples: int = 300,
    seed: int = 42,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Generate synthetic clinical annotations.

    Args:
        n_samples: Number of samples (must match expression).
        seed: Random seed.
        output_dir: If set, write TSV.

    Returns:
        Clinical DataFrame with sample_id index.
    """
    rng = np.random.RandomState(seed)
    sample_ids = [f"TCGA-SYNTH-{i:04d}" for i in range(n_samples)]

    # Stages
    stages = rng.choice(
        ["Stage I", "Stage II", "Stage III", "Stage IV"],
        size=n_samples, p=[0.15, 0.25, 0.40, 0.20],
    )

    # Survival (months)
    os_time = np.round(rng.exponential(scale=30, size=n_samples), 1)
    os_event = rng.binomial(1, 0.4, size=n_samples)

    # Tissue type
    tissue = rng.choice(["tumor", "normal"], size=n_samples, p=[0.85, 0.15])

    df = pd.DataFrame({
        "sample_id": sample_ids,
        "stage": stages,
        "os_time": os_time,
        "os_event": os_event,
        "tissue_type": tissue,
    }).set_index("sample_id")

    if output_dir:
        path = os.path.join(output_dir, "tcga_stad_clinical_synthetic.tsv")
        df.to_csv(path, sep="\t")
        log(f"  Saved: {path}")

    return df


def generate_nk_state_labels(
    expression_df: pd.DataFrame,
    seed: int = 42,
    output_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Generate synthetic NK immune state labels from expression data.

    Uses a subset of NK-related genes to create realistic-looking labels
    with known ground-truth separability.

    Args:
        expression_df: Bulk expression (samples × genes).
        seed: Random seed.
        output_dir: If set, write TSV.

    Returns:
        DataFrame with nk_immune_state and NK_dysfunction_score columns.
    """
    rng = np.random.RandomState(seed)

    # Build a synthetic dysfunction score from known marker genes
    cyto_genes = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG"]
    dysf_genes = ["HAVCR2", "TIGIT", "CD96", "KLRC1"]

    avail_cyto = [g for g in cyto_genes if g in expression_df.columns]
    avail_dysf = [g for g in dysf_genes if g in expression_df.columns]

    if avail_cyto and avail_dysf:
        cyto_score = expression_df[avail_cyto].mean(axis=1)
        dysf_score = expression_df[avail_dysf].mean(axis=1)
        # Add noise
        dysf_score += rng.randn(len(expression_df)) * 0.5
    else:
        dysf_score = pd.Series(rng.randn(len(expression_df)), index=expression_df.index)

    # Quantile-based state assignment
    q25 = dysf_score.quantile(0.25)
    q50 = dysf_score.quantile(0.50)
    q75 = dysf_score.quantile(0.75)

    def assign_state(s):
        if s < q25:
            return "NK-hot-cytotoxic"
        elif s < q50:
            return "NK-hot-dysfunctional"
        elif s < q75:
            return "NK-intermediate"
        else:
            return "NK-cold/excluded"

    states = dysf_score.apply(assign_state)

    df = pd.DataFrame({
        "nk_immune_state": states,
        "NK_dysfunction_score": dysf_score.values,
    }, index=expression_df.index)

    if output_dir:
        path = os.path.join(output_dir, "nk_state_labels_synthetic.tsv")
        df.to_csv(path, sep="\t")
        log(f"  Saved: {path}")

    return df


def generate_scrn_data(
    n_cells: int = 5000,
    n_genes: int = 100,
    seed: int = 42,
    output_dir: Optional[str] = None,
):
    """Generate synthetic scRNA-seq AnnData object.

    Creates cell types: NK, T, B, Myeloid, Epithelial, Stromal.
    NK cells have subtypes: cytotoxic vs dysfunctional.

    Args:
        n_cells: Total number of cells.
        n_genes: Number of genes.
        seed: Random seed.
        output_dir: If set, write .h5ad file.

    Returns:
        AnnData object (or None if anndata not installed).
    """
    try:
        import anndata
    except ImportError:
        log("  anndata not installed — skipping scRNA generation")
        return None

    rng = np.random.RandomState(seed)
    n_genes = min(n_genes, len(ALL_GENES))
    genes = ALL_GENES[:n_genes]

    # Cell type proportions
    cell_types = ["NK", "T_CD8", "T_CD4", "B", "Myeloid", "Epithelial", "Stromal"]
    proportions = [0.08, 0.15, 0.10, 0.10, 0.12, 0.30, 0.15]
    n_per_type = [max(1, int(n_cells * p)) for p in proportions]
    n_per_type[0] += n_cells - sum(n_per_type)  # adjust rounding

    # Generate expression per cell type with type-specific patterns
    all_exprs = []
    all_obs = []

    nk_dysf_genes = ["HAVCR2", "TIGIT", "CD96", "KLRC1"]
    nk_cyto_genes = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG"]

    for ct, n_ct in zip(cell_types, n_per_type):
        # Base expression
        expr = rng.lognormal(mean=1.5, sigma=0.8, size=(n_ct, n_genes))

        # Type-specific signals
        if ct == "NK":
            # NK: high cytotoxicity / receptor genes
            for g in nk_cyto_genes + nk_dysf_genes:
                if g in genes:
                    gi = genes.index(g)
                    expr[:, gi] += rng.lognormal(mean=1.5, sigma=0.3, size=n_ct)
        elif ct == "T_CD8":
            for g in ["CD8A", "GZMB", "PRF1"]:
                if g in genes:
                    expr[:, gi] += rng.lognormal(mean=1.0, sigma=0.3, size=n_ct)
        elif ct == "Epithelial":
            for g in ["EPCAM", "KRT19", "CDH1"]:
                if g in genes:
                    gi = genes.index(g)
                    expr[:, gi] += rng.lognormal(mean=2.0, sigma=0.3, size=n_ct)

        expr = np.clip(expr, 0, 12)
        all_exprs.append(expr)

        # Build observation metadata
        obs = pd.DataFrame({
            "cell_type": [ct] * n_ct,
            "condition": rng.choice(
                ["tumor", "normal"], size=n_ct, p=[0.7, 0.3]
            ),
            "sample_id": rng.choice(
                [f"S{i:02d}" for i in range(10)], size=n_ct,
            ),
        })
        # Add NK subtypes
        if ct == "NK":
            obs["nk_subtype"] = rng.choice(
                ["NK_cytotoxic", "NK_dysfunctional", "NK_tissue_resident"],
                size=n_ct, p=[0.4, 0.35, 0.25],
            )
        all_obs.append(obs)

    X = np.vstack(all_exprs).astype(np.float32)
    obs = pd.concat(all_obs, ignore_index=True)
    var = pd.DataFrame(index=genes)

    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.var_names = genes

    if output_dir:
        path = os.path.join(output_dir, "gc_integrated_synthetic.h5ad")
        adata.write(path)
        log(f"  Saved: {path} ({adata.n_obs} cells × {adata.n_vars} genes)")

        # Also write an NK subset
        nk_sub = adata[adata.obs["cell_type"] == "NK"].copy()
        path_nk = os.path.join(output_dir, "gc_nk_subset_synthetic.h5ad")
        nk_sub.write(path_nk)
        log(f"  Saved: {path_nk} ({nk_sub.n_obs} NK cells)")

    return adata


def generate_prior_networks(output_dir: str, seed: int = 42) -> None:
    """Generate minimal prior network files for graph construction.

    Creates placeholder STRING PPI, CellChatDB LR, and ChEA TF-target files
    covering the SST gene set so the graph builder can run without real data.
    """
    rng = np.random.RandomState(seed)
    genes = SST_GENES + CANDIDATE_GENES
    unique_genes = sorted(set(genes))

    # --- STRING PPI (minimal) ---
    ppi_rows = []
    for i, g1 in enumerate(unique_genes):
        for g2 in unique_genes[i + 1:i + 4]:
            if rng.rand() < 0.3:  # 30% edge density
                ppi_rows.append({
                    "protein1": g1,
                    "protein2": g2,
                    "score": round(float(rng.uniform(700, 999)), 1),
                })
    if ppi_rows:
        pd.DataFrame(ppi_rows).to_csv(
            os.path.join(output_dir, "string_ppi_physical.tsv"),
            sep=" ", index=False,
        )

    # --- CellChatDB LR ---
    lr_rows = []
    tumor_side = ["TGFB1", "MICA", "MICB", "PVR", "NECTIN2", "HLA-E"]
    nk_side = ["TGFBR1", "TGFBR2", "KLRK1", "KLRD1", "KLRC1", "CD96", "TIGIT", "HAVCR2"]
    for l in tumor_side:
        for r in nk_side:
            if rng.rand() < 0.15:
                lr_rows.append({"ligand": l, "receptor": r})
    if lr_rows:
        pd.DataFrame(lr_rows).to_csv(
            os.path.join(output_dir, "cellchatdb_human.csv"), index=False,
        )

    # --- ChEA TF-target ---
    tf_rows = []
    tfs = ["STAT1", "STAT3", "NFKB1", "RELA", "IRF1", "IRF4", "TBX21", "EOMES"]
    for tf in tfs:
        targets = [g for g in unique_genes if rng.rand() < 0.1]
        if targets:
            tf_rows.append("\t".join([tf] + targets))
    if tf_rows:
        with open(os.path.join(output_dir, "chea_tf_targets.txt"), "w") as f:
            f.write("\n".join(tf_rows))

    log(f"  Prior networks: {len(ppi_rows)} PPI, {len(lr_rows)} LR, {len(tf_rows)} TF-target")


def generate_all(
    output_dir: str = "data/synthetic",
    n_samples: int = 300,
    n_cells: int = 5000,
    n_genes: int = 100,
    seed: int = 42,
) -> Dict[str, str]:
    """Generate a complete synthetic dataset for end-to-end pipeline testing.

    Returns:
        Dict mapping dataset name to file path.
    """
    log("=" * 60)
    log("GENERATING SYNTHETIC DATASET")
    log(f"  Samples: {n_samples}, Cells: {n_cells}, Genes: {n_genes}, Seed: {seed}")
    log("=" * 60)

    d = ensure_dir(output_dir)
    prior_dir = ensure_dir(os.path.join(output_dir, "prior_networks"))
    files: Dict[str, str] = {}

    # 1. Bulk expression
    log("\n[1/6] Bulk expression...")
    expr = generate_bulk_expression(
        n_samples=n_samples, n_genes=n_genes, seed=seed, output_dir=d
    )
    files["tcga_stad_expression"] = os.path.join(d, "tcga_stad_expression_synthetic.tsv")

    # 2. Clinical data
    log("\n[2/6] Clinical data...")
    generate_clinical_data(n_samples=n_samples, seed=seed, output_dir=d)
    files["tcga_stad_clinical"] = os.path.join(d, "tcga_stad_clinical_synthetic.tsv")

    # 3. NK state labels
    log("\n[3/6] NK state labels...")
    generate_nk_state_labels(expr, seed=seed, output_dir=d)
    files["nk_state_labels"] = os.path.join(d, "nk_state_labels_synthetic.tsv")

    # 4. scRNA data
    log("\n[4/6] scRNA-seq data...")
    generate_scrn_data(n_cells=n_cells, n_genes=n_genes, seed=seed, output_dir=d)
    files["gc_integrated"] = os.path.join(d, "gc_integrated_synthetic.h5ad")
    files["gc_nk_subset"] = os.path.join(d, "gc_nk_subset_synthetic.h5ad")

    # 5. Prior networks
    log("\n[5/6] Prior networks...")
    generate_prior_networks(str(prior_dir), seed=seed)
    files["prior_networks"] = str(prior_dir)

    # 6. SST-axis scores (from scRNA)
    log("\n[6/6] Computing SST-axis scores on synthetic scRNA...")
    try:
        import scanpy as sc
        nk_path = os.path.join(d, "gc_nk_subset_synthetic.h5ad")
        adata = sc.read(nk_path)
        from src.topology.sst_axis import compute_sst_scores
        scores = compute_sst_scores(adata)

        score_dir = ensure_dir(os.path.join(d, "results", "tables"))
        scores.to_csv(
            os.path.join(score_dir, "sst_axis_scores_single_cell.tsv"),
            sep="\t", index_label="cell_id",
        )
        log("  SST scores computed and saved")
    except Exception as e:
        log(f"  SST scores skipped ({e})")

    log("\n" + "=" * 60)
    log("SYNTHETIC DATA GENERATION COMPLETE")
    for k, v in files.items():
        log(f"  {k}: {v}")
    log("=" * 60)

    return files


# =========================================================================
# CLI
# =========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate synthetic data for GC-NKGraph-Atlas pipeline testing"
    )
    parser.add_argument("--output-dir", default="data/synthetic",
                       help="Output directory (default: data/synthetic)")
    parser.add_argument("--n-samples", type=int, default=300,
                       help="Number of bulk samples (default: 300)")
    parser.add_argument("--n-cells", type=int, default=5000,
                       help="Number of scRNA cells (default: 5000)")
    parser.add_argument("--n-genes", type=int, default=100,
                       help="Number of genes (default: 100)")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed (default: 42)")
    args = parser.parse_args()

    generate_all(
        output_dir=args.output_dir,
        n_samples=args.n_samples,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
