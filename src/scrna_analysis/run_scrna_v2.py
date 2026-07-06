#!/usr/bin/env python
"""
GC-NKGraph-Atlas scRNA Pipeline v2 (Phase 3-4).
Simplified, robust version with verbose logging.
"""
import os, sys, gzip, tarfile, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

LOG = None

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def load_10x_csv(fpath, sample_name):
    log(f"  Loading {fpath}...")
    df = pd.read_csv(fpath, index_col=0)
    df.columns = [f"{sample_name}_{c}" for c in df.columns]
    adata = sc.AnnData(df.T)
    adata.obs["sample_id"] = sample_name
    return adata

def checkpoint(path, adata, msg):
    """Save checkpoint and log progress."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write(path)
    log(f"  CHECKPOINT: {msg} -> {path}")

def load_checkpoint(path):
    """Load checkpoint if it exists."""
    if os.path.exists(path):
        log(f"  Found checkpoint: {path}")
        return sc.read(path)
    return None

def main():
    log("=" * 50)
    log("GC-NKGraph-Atlas scRNA Pipeline v2")
    log("CHECKPOINT-RESUME enabled")
    log("=" * 50)

    # Config
    raw_dir = "data/raw/scrna/gse246662"
    out_dir = "data/processed/scrna"
    fig_dir = "results/figures"
    tbl_dir = "results/tables"
    ckpt_dir = os.path.join(raw_dir, "checkpoints")  # checkpoints survive restarts
    os.makedirs(ckpt_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(tbl_dir, exist_ok=True)

    # Sample metadata
    meta = {
        "GSM7874169": ("HL1", "healthy_liver", "normal"),
        "GSM7874170": ("HL2", "healthy_liver", "normal"),
        "GSM7874171": ("HL3", "healthy_liver", "normal"),
        "GSM7874172": ("GC1", "gastric_cancer", "tumor"),
        "GSM7874173": ("GC2", "gastric_cancer", "tumor"),
        "GSM7874174": ("GC3", "gastric_cancer", "tumor"),
        "GSM7874175": ("LM1", "liver_metastasis", "tumor"),
        "GSM7874176": ("LM2", "liver_metastasis", "tumor"),
        "GSM7874177": ("LM3", "liver_metastasis", "tumor"),
    }

    # ---- CHECKPOINT: raw_loaded ----
    ckpt_raw = os.path.join(ckpt_dir, "01_raw_loaded.h5ad")
    adata = load_checkpoint(ckpt_raw)

    if adata is None:
        # Step 1: Extract
        log("Step 1: Extracting tar...")
        extract_dir = os.path.join(raw_dir, "extracted_v2")
        if not os.path.exists(extract_dir):
            os.makedirs(extract_dir)
            with tarfile.open(os.path.join(raw_dir, "GSE246662_RAW.tar")) as tar:
                tar.extractall(path=extract_dir)
            log(f"  Extracted to {extract_dir}")
        else:
            log(f"  Already extracted: {len(os.listdir(extract_dir))} files")

        # Step 2: Load samples
        log("Step 2: Loading samples...")
        adatas = []
        for fname in sorted(os.listdir(extract_dir)):
            if not fname.endswith(".csv.gz"):
                continue
            gsm = fname.split("_")[0]
            if gsm not in meta:
                continue
            sid, tissue, cond = meta[gsm]
            fpath = os.path.join(extract_dir, fname)
            adata = load_10x_csv(fpath, sid)
            adata.obs["tissue"] = tissue
            adata.obs["condition"] = cond
            adata.obs["gsm"] = gsm
            adatas.append(adata)
            log(f"  {sid}: {adata.n_obs} cells x {adata.n_vars} genes")

        # Concatenate (inner join: only common genes across all samples)
        log("Step 2b: Concatenating...")
        adata = sc.concat(adatas, axis=0, join="inner")
        adata.obs_names_make_unique()
        log(f"  Total: {adata.n_obs} cells x {adata.n_vars} genes")
        # Free memory
        del adatas
        checkpoint(ckpt_raw, adata, "raw data loaded")
    else:
        log(f"  [resumed] cells={adata.n_obs} genes={adata.n_vars}")

    # ---- CHECKPOINT: qc_filtered ----
    ckpt_qc = os.path.join(ckpt_dir, "02_qc_filtered.h5ad")
    ckpt_norm = os.path.join(ckpt_dir, "03_normalized.h5ad")

    if os.path.exists(ckpt_qc):
        adata = load_checkpoint(ckpt_qc)
        log(f"  [resumed from QC] cells={adata.n_obs}")
    elif os.path.exists(ckpt_raw):
        adata = load_checkpoint(ckpt_raw)
        log(f"  [resumed from raw] cells={adata.n_obs}, running QC")

    # Step 3: QC (manual calculation for robustness)
    if not os.path.exists(ckpt_qc):
        log("Step 3: QC filtering...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Manual QC metrics (works across scanpy versions)
    import scipy.sparse
    X = adata.X
    if scipy.sparse.issparse(X):
        X_dense = X.toarray() if hasattr(X, "toarray") else X.A
    else:
        X_dense = X

    adata.obs["n_genes"] = (X_dense > 0).sum(axis=1)
    adata.obs["total_counts"] = X_dense.sum(axis=1)
    mt_mask = adata.var["mt"].values
    if mt_mask.any():
        adata.obs["pct_counts_mt"] = (X_dense[:, mt_mask].sum(axis=1) / adata.obs["total_counts"]) * 100
    else:
        adata.obs["pct_counts_mt"] = 0.0

    log(f"  Before QC: {adata.n_obs} cells")
    qc_cols = [c for c in adata.obs.columns if "n_genes" in c or "total" in c or "mt" in c]
    log(f"  QC columns: {qc_cols}")
    log(f"  n_genes stats: min={adata.obs['n_genes'].min():.0f} median={adata.obs['n_genes'].median():.0f} max={adata.obs['n_genes'].max():.0f}")
    log(f"  total_counts stats: min={adata.obs['total_counts'].min():.0f} median={adata.obs['total_counts'].median():.0f} max={adata.obs['total_counts'].max():.0f}")

    # Use quantile-based filtering (more robust for sparse scRNA data)
    from scipy.stats import median_abs_deviation
    for col, name, low_q, high_q in [("n_genes", "n_genes", 0.01, 0.99),
                                       ("total_counts", "total_counts", 0.01, 0.99)]:
        lo = adata.obs[col].quantile(low_q)
        hi = adata.obs[col].quantile(high_q)
        adata.obs[f"pass_{name}"] = (adata.obs[col] >= lo) & (adata.obs[col] <= hi)
    adata.obs["pass_mt"] = adata.obs["pct_counts_mt"] < 20
    prev_n = adata.n_obs
    adata = adata[adata.obs[["pass_n_genes","pass_total_counts","pass_mt"]].all(axis=1)].copy()
    log(f"  After QC: {adata.n_obs} cells")
    checkpoint(ckpt_qc, adata, f"QC filtered: {prev_n} -> {adata.n_obs} cells")

    # ---- CHECKPOINT check for normalized data ----
    if os.path.exists(ckpt_norm):
        adata = load_checkpoint(ckpt_norm)
        log(f"  [resumed from normalized] cells={adata.n_obs} genes={adata.n_vars}")
    elif os.path.exists(ckpt_qc):
        adata = load_checkpoint(ckpt_qc)
        log(f"  [resumed from QC]")

    if not os.path.exists(ckpt_norm):
        # Step 4: Normalize + HVG
        log("Step 4: Normalization + HVG...")
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        checkpoint(ckpt_norm, adata, "normalized")

        # Simple variance-based HVG (most robust across data shapes)
        log("  Selecting HVGs by variance...")
    import scipy.sparse
    X = adata.X
    if scipy.sparse.issparse(X):
        X_d = X.toarray() if hasattr(X, "toarray") else X.A
    else:
        X_d = X
    gene_vars = np.var(X_d, axis=0)
    gene_vars = np.nan_to_num(gene_vars, nan=0.0)
    n_top = min(3000, len(gene_vars))
    top_idx = np.argsort(gene_vars)[-n_top:]
    adata.var["highly_variable"] = False
    adata.var.iloc[top_idx, adata.var.columns.get_loc("highly_variable")] = True
    log(f"  HVG: {adata.var.highly_variable.sum()} genes selected (by variance)")

    # Step 5: scVI
    ckpt_scvi = os.path.join(ckpt_dir, "04_scvi.h5ad")
    if os.path.exists(ckpt_scvi):
        adata = load_checkpoint(ckpt_scvi)
        log(f"  [resumed from scVI] cells={adata.n_obs}")
    else:
        if not os.path.exists(ckpt_norm):
            # hvg was computed on the normalized version
            pass
        log("Step 5: scVI integration...")
        import scvi
        scvi.settings.seed = 42
        scvi.settings.verbosity = 0

        adata_hvg = adata[:, adata.var.highly_variable].copy()
        scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample_id")
        model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30)
        log("  Training scVI (this takes ~10-20 min)...")
        model.train(max_epochs=200, early_stopping=True)
        adata.obsm["X_scVI"] = model.get_latent_representation()
        log(f"  scVI done: {adata.obsm['X_scVI'].shape}")
        checkpoint(ckpt_scvi, adata, "scVI integrated")

    # Step 6: UMAP + Clustering
    log("Step 6: UMAP + Leiden...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=1.0)
    log(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

    # Step 7: Cell-type annotation
    log("Step 7: Cell-type annotation...")
    nk_genes = ["NCAM1","KLRD1","NKG7","GNLY","KLRF1","EOMES","NCR1","FCGR3A"]
    t_genes  = ["CD3D","CD3E","CD3G","CD4","CD8A"]
    mono_genes = ["CD14","CD68","CSF1R"]
    b_genes = ["MS4A1","CD79A","CD19"]

    for name, genes in [("NK", nk_genes), ("T", t_genes), ("Mono", mono_genes), ("B", b_genes)]:
        present = [g for g in genes if g in adata.var_names]
        if present:
            adata.obs[f"{name}_score"] = np.array(
                adata[:, present].X.mean(axis=1)).flatten()

    def annotate(r):
        if r["NK_score"] > 0.3 and r["T_score"] < 0.2: return "NK"
        if r["T_score"] > 0.3: return "T_cell"
        if r["Mono_score"] > 0.3: return "Monocyte"
        if r["B_score"] > 0.3: return "B_cell"
        return "Other"
    adata.obs["cell_type"] = adata.obs.apply(annotate, axis=1)
    log(f"  Cell types: {dict(adata.obs['cell_type'].value_counts())}")

    # Step 8: Save
    log("Step 8: Saving outputs...")
    adata.write(os.path.join(out_dir, "gc_integrated.h5ad"))
    log(f"  Saved integrated: {os.path.join(out_dir, 'gc_integrated.h5ad')}")

    nk = adata[adata.obs["cell_type"] == "NK"].copy()
    nk.write(os.path.join(out_dir, "gc_nk_subset.h5ad"))
    log(f"  Saved NK subset: {nk.n_obs} cells")

    summary = adata.obs.groupby("sample_id").agg(
        total=("cell_type","count"),
        nk=("cell_type",lambda x: (x=="NK").sum()),
    ).reset_index()
    summary.to_csv(os.path.join(tbl_dir, "gc_scrna_dataset_summary.tsv"), sep="\t", index=False)
    log(f"  Saved summary table")

    log("=" * 50)
    log("PIPELINE COMPLETE!")
    log("=" * 50)

if __name__ == "__main__":
    main()
