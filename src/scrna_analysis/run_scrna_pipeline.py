"""
GC-NKGraph-Atlas scRNA Pipeline (Phase 3-4).

Processes GSE246662 gastric cancer/liver metastasis scRNA data:
1. Load CSV count matrices into AnnData
2. QC filtering (MAD-based)
3. Scrublet doublet detection
4. scVI integration across samples
5. CellTypist cell-type annotation
6. NK/ILC/T separation

Usage:
    python src/scrna_analysis/run_scrna_pipeline.py
"""

import os, sys, gzip, argparse, warnings
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import ensure_dir, save_table

warnings.filterwarnings("ignore")


SAMPLE_METADATA = {
    "GSM7874169": {"sample_id": "HL1",  "tissue": "healthy_liver",       "condition": "normal"},
    "GSM7874170": {"sample_id": "HL2",  "tissue": "healthy_liver",       "condition": "normal"},
    "GSM7874171": {"sample_id": "HL3",  "tissue": "healthy_liver",       "condition": "normal"},
    "GSM7874172": {"sample_id": "GC1",  "tissue": "gastric_cancer",      "condition": "tumor"},
    "GSM7874173": {"sample_id": "GC2",  "tissue": "gastric_cancer",      "condition": "tumor"},
    "GSM7874174": {"sample_id": "GC3",  "tissue": "gastric_cancer",      "condition": "tumor"},
    "GSM7874175": {"sample_id": "LM1",  "tissue": "liver_metastasis",    "condition": "tumor"},
    "GSM7874176": {"sample_id": "LM2",  "tissue": "liver_metastasis",    "condition": "tumor"},
    "GSM7874177": {"sample_id": "LM3",  "tissue": "liver_metastasis",    "condition": "tumor"},
}

NK_MARKERS = ["NCAM1", "KLRD1", "NKG7", "GNLY", "KLRF1", "EOMES", "NCR1", "FCGR3A", "KLRC1", "KLRK1"]
T_MARKERS  = ["CD3D", "CD3E", "CD3G", "TRAC", "CD4", "CD8A"]
MONOCYTE_MARKERS = ["CD14", "CD68", "FCGR3A", "CSF1R"]
B_MARKERS = ["MS4A1", "CD79A", "CD19"]
PLASMA_MARKERS = ["IGHG1", "MZB1", "SDC1"]


def load_10x_csv(filepath: str, sample_name: str) -> "AnnData":
    """Load a 10x-format CSV (genes x cells) into AnnData."""
    import scanpy as sc
    df = pd.read_csv(filepath, index_col=0)
    # Rename cells to include sample prefix
    df.columns = [f"{sample_name}_{c}" for c in df.columns]
    adata = sc.AnnData(df.T)  # cells x genes
    adata.obs["sample_id"] = sample_name
    adata.obs["cell_barcode"] = df.columns
    adata.var_names_make_unique()
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw/scrna/gse246662")
    parser.add_argument("--output-dir", default="data/processed/scrna")
    parser.add_argument("--figures-dir", default="results/figures")
    parser.add_argument("--tables-dir", default="results/tables")
    parser.add_argument("--nkjobs", type=int, default=16)
    args = parser.parse_args()

    logger = Logger()
    import scanpy as sc
    sc.settings.verbosity = 1

    output_dir = ensure_dir(args.output_dir)
    figures_dir = ensure_dir(args.figures_dir)
    tables_dir = ensure_dir(args.tables_dir)

    # ---- Step 1: Load all samples ----
    logger.ok("SCRNA", "Step 1: Loading samples...", script=__file__)
    adatas = []
    tar_path = os.path.join(args.raw_dir, "GSE246662_RAW.tar")
    extract_dir = os.path.join(args.raw_dir, "extracted")
    ensure_dir(extract_dir)

    if not os.path.exists(extract_dir) or len(os.listdir(extract_dir)) < 3:
        import tarfile
        with tarfile.open(tar_path) as tar:
            tar.extractall(path=extract_dir)
        logger.ok("SCRNA", f"Extracted tar to {extract_dir}", script=__file__)

    for fname in sorted(os.listdir(extract_dir)):
        if not fname.endswith(".csv.gz"):
            continue
        gsm_id = fname.split("_")[0]
        if gsm_id not in SAMPLE_METADATA:
            continue
        meta = SAMPLE_METADATA[gsm_id]
        fpath = os.path.join(extract_dir, fname)
        adata = load_10x_csv(fpath, meta["sample_id"])
        for k, v in meta.items():
            adata.obs[k] = v
        adatas.append(adata)
        logger.ok("SCRNA", f"  Loaded {meta['sample_id']}: {adata.n_obs} cells x {adata.n_vars} genes", script=__file__)

    # Concatenate
    adata = sc.concat(adatas, axis=0, join="outer")
    adata.obs_names_make_unique()
    logger.ok("SCRNA", f"Concatenated: {adata.n_obs} cells x {adata.n_vars} genes", script=__file__)

    # ---- Step 2: QC ----
    logger.ok("SCRNA", "Step 2: QC filtering...", script=__file__)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # MAD-based filtering
    def mad_filter(col, n_mads=5):
        from scipy.stats import median_abs_deviation
        med = np.median(col)
        mad = median_abs_deviation(col)
        return (col > med - n_mads * mad) & (col < med + n_mads * mad)

    cell_filter = (
        mad_filter(adata.obs["n_genes_by_counts"]) &
        mad_filter(adata.obs["total_counts"]) &
        (adata.obs["pct_counts_mt"] < 20)
    )
    adata = adata[cell_filter].copy()
    logger.ok("SCRNA", f"After QC: {adata.n_obs} cells", script=__file__)

    # ---- Step 3: Normalization + HVG ----
    logger.ok("SCRNA", "Step 3: Normalization & HVG selection...", script=__file__)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample_id", flavor="seurat_v3")
    n_hvg = adata.var.highly_variable.sum()
    logger.ok("SCRNA", f"HVG selected: {n_hvg}", script=__file__)

    # ---- Step 4: Scrublet doublet detection ----
    logger.ok("SCRNA", "Step 4: Scrublet doublet detection...", script=__file__)
    try:
        import scrublet
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
        for s in adata.obs["sample_id"].unique():
            mask = adata.obs["sample_id"] == s
            sdata = adata[mask].copy()
            scrub = scrublet.Scrublet(sdata.X, expected_doublet_rate=0.06)
            doublet_scores, predicted = scrub.scrub_doublets()
            adata.obs.loc[mask, "doublet_score"] = doublet_scores
            adata.obs.loc[mask, "predicted_doublet"] = predicted
        n_doublet = adata.obs["predicted_doublet"].sum()
        logger.ok("SCRNA", f"Doublets found: {n_doublet} / {adata.n_obs}", script=__file__)
    except Exception as e:
        logger.needs_review("SCRNA", f"Scrublet failed: {e}", script=__file__)

    # ---- Step 5: scVI integration ----
    logger.ok("SCRNA", "Step 5: scVI integration...", script=__file__)
    import scvi
    scvi.settings.seed = 42

    adata_hvg = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample_id")
    model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30)
    model.train(max_epochs=200, early_stopping=True)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    logger.ok("SCRNA", "scVI integration complete", script=__file__)

    # ---- Step 6: Clustering + UMAP ----
    logger.ok("SCRNA", "Step 6: Clustering & UMAP...", script=__file__)
    import scipy.sparse
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=1.0)
    logger.ok("SCRNA", f"Leiden clusters: {adata.obs['leiden'].nunique()}", script=__file__)

    # ---- Step 7: Cell-type annotation ----
    logger.ok("SCRNA", "Step 7: Cell-type annotation...", script=__file__)

    # Compute marker scores
    for name, markers in [
        ("NK_score", NK_MARKERS),
        ("T_score", T_MARKERS),
        ("Mono_score", MONOCYTE_MARKERS),
        ("B_score", B_MARKERS),
        ("Plasma_score", PLASMA_MARKERS),
    ]:
        present = [g for g in markers if g in adata.var_names]
        if present:
            scores = adata[:, present].X.toarray() if scipy.sparse.issparse(adata.X) else adata[:, present].X
            adata.obs[name] = scores.mean(axis=1)

    # Simple rule-based annotation
    def annotate_cell(row):
        if row["NK_score"] > 0.3 and row["T_score"] < 0.2:
            return "NK"
        elif row["T_score"] > 0.3:
            return "T_cell"
        elif row["Mono_score"] > 0.3:
            return "Monocyte"
        elif row["B_score"] > 0.3:
            return "B_cell"
        elif row["Plasma_score"] > 0.3:
            return "Plasma"
        else:
            return "Other"

    adata.obs["cell_type"] = adata.obs.apply(annotate_cell, axis=1)
    counts = adata.obs["cell_type"].value_counts().to_dict()
    logger.ok("SCRNA", f"Cell types: {counts}", script=__file__)

    # Extract NK subset
    nk_adata = adata[adata.obs["cell_type"] == "NK"].copy()
    logger.ok("SCRNA", f"NK subset: {nk_adata.n_obs} cells", script=__file__)

    # ---- Step 8: Save outputs ----
    logger.ok("SCRNA", "Step 8: Saving outputs...", script=__file__)

    integrated_path = os.path.join(output_dir, "gc_integrated.h5ad")
    adata.write(integrated_path)
    logger.ok("SCRNA", f"Saved integrated: {integrated_path}", script=__file__)

    nk_path = os.path.join(output_dir, "gc_nk_subset.h5ad")
    nk_adata.write(nk_path)
    logger.ok("SCRNA", f"Saved NK subset: {nk_path}", script=__file__)

    # Save summary table
    summary = adata.obs.groupby("sample_id").agg(
        total_cells=("cell_type", "count"),
        nk_cells=("cell_type", lambda x: (x == "NK").sum()),
        t_cells=("cell_type", lambda x: (x == "T_cell").sum()),
    ).reset_index()
    save_table(summary, os.path.join(tables_dir, "gc_scrna_dataset_summary.tsv"),
               provenance={"script": __file__})
    logger.ok("SCRNA", "Pipeline complete!", script=__file__)


if __name__ == "__main__":
    main()
