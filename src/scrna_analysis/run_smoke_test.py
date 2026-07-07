#!/usr/bin/env python
"""
GC-NKGraph-Atlas scRNA smoke test.
Runs end-to-end on 2 samples x 500 cells each to validate the pipeline.
"""
import os, sys, tarfile, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scanpy as sc
import pandas as pd
import numpy as np

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    log("=" * 50)
    log("SMOKE TEST: scRNA Pipeline")
    log("=" * 50)

    raw_dir = "data/raw/scrna/gse246662"
    out_dir = "data/processed/scrna/smoke_test"
    os.makedirs(out_dir, exist_ok=True)

    # Extract just 2 samples
    log("Step 1: Extracting 2 test samples...")
    extract_dir = os.path.join(raw_dir, "smoke_extract")
    if not os.path.exists(extract_dir):
        os.makedirs(extract_dir)
        with tarfile.open(os.path.join(raw_dir, "GSE246662_RAW.tar")) as tar:
            members = [m for m in tar.getmembers() if "HL1" in m.name or "GC1" in m.name]
            for m in members:
                tar.extract(m, path=extract_dir)
        log(f"  Extracted {len(os.listdir(extract_dir))} files")

    # Load 2 samples, subset 500 cells each
    log("Step 2: Loading samples (500 cells each)...")
    adatas = []
    for fname in sorted(os.listdir(extract_dir)):
        if not fname.endswith(".csv.gz"):
            continue
        sid = "HL1" if "HL1" in fname else "GC1"
        tissue = "healthy_liver" if "HL" in sid else "gastric_cancer"
        fpath = os.path.join(extract_dir, fname)
        df = pd.read_csv(fpath, index_col=0)
        df = df.iloc[:, :500]  # subset to 500 cells
        df.columns = [f"{sid}_{c}" for c in df.columns]
        adata = sc.AnnData(df.T)
        adata.obs["sample_id"] = sid
        adata.obs["tissue"] = tissue
        adatas.append(adata)
        log(f"  {sid}: {adata.n_obs} cells x {adata.n_vars} genes")

    adata = sc.concat(adatas, axis=0, join="inner")
    adata.obs_names_make_unique()
    log(f"  Total: {adata.n_obs} cells x {adata.n_vars} genes")

    # Step 3: QC
    log("Step 3: QC filtering...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_genes = (adata.X > 0).sum(axis=1)
    total_counts = adata.X.sum(axis=1)
    lo, hi = np.percentile(n_genes, [1, 99])
    adata = adata[(n_genes >= lo) & (n_genes <= hi)].copy()
    log(f"  After QC: {adata.n_obs} cells")

    # Step 4: Norm + HVG
    log("Step 4: Normalization + HVG...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", n_bins=20)
    hvg_n = adata.var.highly_variable.sum()
    log(f"  HVG: {hvg_n}")
    if hvg_n == 0:
        log("  HVG=0, using variance fallback...")
        gene_vars = np.var(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X, axis=0)
        adata.var["highly_variable"] = False
        adata.var.iloc[np.argsort(gene_vars)[-2000:], adata.var.columns.get_loc("highly_variable")] = True
        hvg_n = adata.var.highly_variable.sum()
        log(f"  HVG fallback: {hvg_n}")

    # Step 5: scVI (3 epochs only for testing)
    log("Step 5: scVI (3 epochs, smoke test)...")
    import scvi
    scvi.settings.seed = 42
    scvi.settings.verbosity = 1
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample_id")
    model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30)
    model.train(max_epochs=3)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    log(f"  scVI OK: {adata.obsm['X_scVI'].shape}")

    # Step 6: UMAP
    log("Step 6: UMAP + Leiden...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=10)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=1.0)
    log(f"  Clusters: {adata.obs['leiden'].nunique()}")

    # Step 7: Annotation (robust)
    log("Step 7: Cell-type annotation...")
    gene_sets = {
        "NK": ["NCAM1","KLRD1","NKG7","GNLY","KLRF1","EOMES"],
        "T": ["CD3D","CD3E","CD3G","CD4","CD8A"],
        "Mono": ["CD14","CD68","CSF1R"],
        "B": ["MS4A1","CD79A","CD19"],
    }
    for name, genes in gene_sets.items():
        present = [g for g in genes if g in adata.var_names]
        adata.obs[f"{name}_score"] = np.array(adata[:, present].X.mean(axis=1)).flatten() if present else 0.0
        log(f"  {name}: {len(present)}/{len(genes)} genes")

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
    log(f"  Types: {dict(adata.obs['cell_type'].value_counts())}")

    # Step 8: Save
    log("Step 8: Saving...")
    adata.write(os.path.join(out_dir, "smoke_integrated.h5ad"))
    nk = adata[adata.obs["cell_type"] == "NK"].copy()
    nk.write(os.path.join(out_dir, "smoke_nk_subset.h5ad"))
    log(f"  NK subset: {nk.n_obs} cells")

    log("=" * 50)
    log("SMOKE TEST PASSED! Full pipeline can proceed.")
    log("=" * 50)

if __name__ == "__main__":
    main()
