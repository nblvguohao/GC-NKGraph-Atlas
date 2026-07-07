"""
Standard scRNA QC filtering (methodological hardening).

PROBLEM
-------
Two failure modes existed:
  - v1 (`run_scrna_pipeline.py`) logged `After QC: 0 cells` — a QC step dropped
    100% of cells (cross-platform NaN in `calculate_qc_metrics`).
  - v2 (`run_scrna_v2.py`) over-corrected: it keeps *all* cells and does NO
    filtering (no min-genes, no mito%, no doublet removal). That is why the NK
    subset exists, but reviewers will reject "no QC" on 166k cross-platform cells.

FIX
---
Apply conventional, defensible per-cell / per-gene QC on the integrated object,
with thresholds that are logged and adjustable. This does NOT re-integrate; run
it on the concatenated raw-counts object *before* normalization/scVI, or on the
saved integrated object to produce a QC'd subset for reporting.

Standard thresholds (adjust per dataset; defaults are conservative and typical):
  - min_genes_per_cell = 200
  - max_genes_per_cell = 6000        (doublet/multiplet proxy)
  - max_pct_mito       = 20.0        (dying cells)
  - min_cells_per_gene = 3
  - optional Scrublet doublet removal (--doublets)

USAGE (run on server)
---------------------
    python src/scrna_analysis/qc_filter.py \
        --in data/processed/scrna/gc_integrated.h5ad \
        --out data/processed/scrna/gc_integrated_qc.h5ad

Emits `results/tables/scrna_qc_summary.tsv` (per-sample before/after counts).

NOTE: untested on the author's workstation (scanpy not installed here). The
thresholds are the review-relevant knobs — record whatever you finalize in
Methods §2.4.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger


def run_qc(adata, args, logger):
    import scanpy as sc

    n_start = adata.n_obs
    per_sample_before = adata.obs.get("sample_id", pd.Series(index=adata.obs_names,
                                                             data="all")).value_counts()

    # mito fraction
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    try:
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)
    except Exception as e:
        # cross-platform NaN guard (the exact bug that zeroed v1): compute manually
        logger.ok("SCRNA", f"calculate_qc_metrics fallback: {e}", script=__file__)
        X = adata.X
        counts = np.asarray(X.sum(axis=1)).flatten()
        ngenes = np.asarray((X > 0).sum(axis=1)).flatten()
        mt_counts = np.asarray(X[:, adata.var["mt"].values].sum(axis=1)).flatten()
        adata.obs["total_counts"] = counts
        adata.obs["n_genes_by_counts"] = ngenes
        adata.obs["pct_counts_mt"] = np.where(counts > 0, 100.0 * mt_counts / counts, 0.0)

    ngenes_col = next((c for c in adata.obs.columns
                       if "n_genes" in c and "counts" in c), "n_genes_by_counts")
    mito_col = next((c for c in adata.obs.columns if "pct_counts_mt" in c),
                    "pct_counts_mt")

    keep = (
        (adata.obs[ngenes_col] >= args.min_genes) &
        (adata.obs[ngenes_col] <= args.max_genes) &
        (adata.obs[mito_col] <= args.max_pct_mito)
    )
    logger.ok("SCRNA",
              f"cell QC: {int(keep.sum())}/{n_start} pass "
              f"(min_genes={args.min_genes}, max_genes={args.max_genes}, "
              f"max_mito={args.max_pct_mito}%)", script=__file__)
    adata = adata[keep.values].copy()

    # gene filter
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    logger.ok("SCRNA", f"after gene filter: {adata.n_vars} genes", script=__file__)

    # optional doublet removal
    if args.doublets:
        try:
            import scrublet as scr
            counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
            scrub = scr.Scrublet(counts)
            scores, predicted = scrub.scrub_doublets()
            adata.obs["doublet_score"] = scores
            adata.obs["predicted_doublet"] = predicted
            n_db = int(np.sum(predicted))
            adata = adata[~np.asarray(predicted)].copy()
            logger.ok("SCRNA", f"scrublet removed {n_db} doublets", script=__file__)
        except Exception as e:
            logger.fail("SCRNA", f"scrublet skipped: {e}", script=__file__)

    per_sample_after = adata.obs.get("sample_id", pd.Series(index=adata.obs_names,
                                                            data="all")).value_counts()
    summary = pd.DataFrame({
        "cells_before": per_sample_before,
        "cells_after": per_sample_after,
    }).fillna(0).astype(int)
    summary["retained_frac"] = (summary["cells_after"] /
                                summary["cells_before"].replace(0, np.nan)).round(3)
    return adata, summary


def main():
    ap = argparse.ArgumentParser(description="Standard scRNA QC filtering")
    ap.add_argument("--in", dest="in_path", required=True, help="input .h5ad")
    ap.add_argument("--out", dest="out_path", default=None, help="output .h5ad")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-genes", type=int, default=6000)
    ap.add_argument("--max-pct-mito", type=float, default=20.0)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--doublets", action="store_true", help="run Scrublet doublet removal")
    args = ap.parse_args()

    logger = Logger()
    import scanpy as sc

    adata = sc.read(args.in_path)
    logger.ok("SCRNA", f"loaded {adata.n_obs} cells x {adata.n_vars} genes from {args.in_path}",
              script=__file__)

    adata, summary = run_qc(adata, args, logger)

    out_path = args.out_path or args.in_path.replace(".h5ad", "_qc.h5ad")
    adata.write(out_path)
    tbl = "results/tables/scrna_qc_summary.tsv"
    os.makedirs("results/tables", exist_ok=True)
    summary.to_csv(tbl, sep="\t")
    logger.ok("SCRNA",
              f"QC done: {adata.n_obs} cells retained -> {out_path}; summary -> {tbl}",
              script=__file__)


if __name__ == "__main__":
    main()
