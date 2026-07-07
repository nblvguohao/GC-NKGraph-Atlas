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

def sparse_sum(X, axis=0):
    """Sum of sparse matrix, handling both dense and sparse."""
    import scipy.sparse
    if scipy.sparse.issparse(X):
        return X.sum(axis=axis)
    return X.sum(axis=axis)

def sparse_nonzero_count(X, axis=1):
    """Count non-zero elements, handling sparse matrices efficiently."""
    import scipy.sparse
    if scipy.sparse.issparse(X):
        return X.getnnz(axis=axis)
    return (X > 0).sum(axis=axis)

def load_10x_csv(fpath, sample_name):
    log(f"  Loading {fpath}...")
    df = pd.read_csv(fpath, index_col=0)
    # Detect orientation: if index looks like gene symbols, data is genes x cells (needs transpose)
    # Gene symbols: short (1-15 chars), alphanumeric with dots or dashes
    # Cell barcodes: 16+ chars, start with AAAC or similar, often have -1 suffix
    first_idx = str(df.index[0])
    is_gene = len(first_idx) < 16 or first_idx.startswith(("ENSG", "AL", "AP0", "AC0"))
    needs_transpose = is_gene  # gene index = genes x cells = needs transpose
    if needs_transpose:
        df.columns = [f"{sample_name}_{c}" for c in df.columns]
        adata = sc.AnnData(df.T)
    else:
        # Already cells x genes: index = cell barcodes, columns = gene symbols
        df.index = [f"{sample_name}_{c}" for c in df.index]
        adata = sc.AnnData(df)
    adata.obs["sample_id"] = sample_name
    return adata

def checkpoint(path, adata, msg, as_sparse=True):
    """Save checkpoint and log progress. Stores as sparse CSMatrix to save space."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if as_sparse and hasattr(adata.X, "toarray"):
        pass  # already sparse
    elif as_sparse:
        import scipy.sparse
        adata.X = scipy.sparse.csr_matrix(adata.X)
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

        # Concatenate (inner join: common genes only, avoids sparse NaN issues)
        log("Step 2b: Concatenating (inner join)...")
        adata = sc.concat(adatas, axis=0, join="inner")
        adata.obs_names_make_unique()
        log(f"  Total: {adata.n_obs} cells x {adata.n_vars} genes")
        # Free memory (skip raw checkpoint - too large, 100+ GB)
        del adatas
        log(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        log(f"  [resumed] cells={adata.n_obs} genes={adata.n_vars}")

    # ---- CHECKPOINT: qc_filtered ----
    ckpt_qc = os.path.join(ckpt_dir, "02_qc_filtered.h5ad")
    ckpt_norm = os.path.join(ckpt_dir, "03_normalized.h5ad")

    if os.path.exists(ckpt_qc):
        adata = load_checkpoint(ckpt_qc)
        log(f"  [resumed from QC] cells={adata.n_obs}")
        ckpt_raw = None  # raw checkpoint no longer saved

    # Step 3: Skip QC filtering (inner-join common genes work across platforms)
    # sc.pp.calculate_qc_metrics has NaN issues with cross-platform data.
    # We keep all cells and rely on scVI to handle biological variation.
    if not os.path.exists(ckpt_qc):
        log("Step 3: QC (minimal - keep all cells, just log metrics)")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    try:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        log(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
        ngenes_col = [c for c in adata.obs.columns if "n_genes" in c and "counts" in c]
        if ngenes_col:
            log(f"  n_genes median: {adata.obs[ngenes_col[0]].median():.0f}")
    except Exception as e:
        log(f"  QC metrics logging failed (non-critical): {e}")
    checkpoint(ckpt_qc, adata, "QC logged (all cells kept)")

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
        log(f"  Normalized: {adata.n_obs} cells")
        checkpoint(ckpt_norm, adata, "normalized")

        # Robust HVG selection with multiple fallbacks
        log("  Selecting HVGs...")
    import scipy.sparse
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", n_bins=20)
        hvg_n = adata.var.highly_variable.sum()
        log(f"  HVG (seurat): {hvg_n}")
        if hvg_n == 0:
            raise ValueError("0 HVG from seurat")
    except Exception:
        log("  HVG seurat failed, using variance fallback")
        import scipy.sparse
        if scipy.sparse.issparse(adata.X):
            X_sub = adata.X.toarray() if adata.X.shape[1] < 10000 else adata.X[:, :10000].toarray()
            # For large sparse, compute variance without dense conversion: Var = E[X²] - E[X]²
            from scipy.sparse import issparse
            X_mean = np.array(adata.X.mean(axis=0)).flatten()
            X_sq_mean = np.array(adata.X.power(2).mean(axis=0)).flatten() if issparse(adata.X) else (adata.X**2).mean(axis=0)
            gene_vars = X_sq_mean - X_mean**2
            gene_vars = np.nan_to_num(gene_vars, nan=0.0, posinf=0.0, neginf=0.0)
        else:
            gene_vars = np.nanvar(adata.X, axis=0)
        gene_vars = np.nan_to_num(gene_vars, nan=0.0)
        adata.var["highly_variable"] = False
        n_genes = min(3000, len(gene_vars))
        if n_genes == 0:
            log("  ERROR: 0 genes in dataset! Using all genes.")
            adata.var["highly_variable"] = True
        elif np.max(gene_vars) == 0:
            log(f"  All {n_genes} genes have 0 variance, using first 3000")
            adata.var.iloc[:n_genes, adata.var.columns.get_loc("highly_variable")] = True
        else:
            top_idx = np.argsort(gene_vars)[-n_genes:]
            adata.var.iloc[top_idx, adata.var.columns.get_loc("highly_variable")] = True
    log(f"  HVG: {adata.var.highly_variable.sum()} genes selected")

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

    gene_sets = {"NK": nk_genes, "T": t_genes, "Mono": mono_genes, "B": b_genes}
    score_columns = []
    for name, genes in gene_sets.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            col = f"{name}_score"
            adata.obs[col] = np.array(adata[:, present].X.mean(axis=1)).flatten()
            score_columns.append(col)
            log(f"  {name}_score: {len(present)}/{len(genes)} genes found")
        else:
            log(f"  {name}_score: NO genes found (all missing after inner join)")

    # Robust annotation: check which score columns exist
    def annotate(r):
        nk = r.get("NK_score", 0)
        t  = r.get("T_score", 0)
        m  = r.get("Mono_score", 0)
        b  = r.get("B_score", 0)
        if nk > 0.3 and t < 0.2: return "NK"
        if t > 0.3: return "T_cell"
        if m > 0.3: return "Monocyte"
        if b > 0.3: return "B_cell"
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
