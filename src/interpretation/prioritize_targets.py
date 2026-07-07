"""
GC-NKGraph-Atlas Candidate Target Prioritization (Phase 13).

Integrates multi-evidence to rank tumor-intrinsic NK-evasion targets.

Evidence dimensions:
  1. SST axis membership (Phase 14R)
  2. Tumor cell specificity (scRNA malignant vs non-malignant)
  3. NK dysfunction / cytotoxicity correlation (scRNA)
  4. Axis direction consistency
  5. Druggability (FDA / clinical / preclinical)
  6. Literature support (gold standard genes)

Usage:
    python src/interpretation/prioritize_targets.py
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.sst_config import (  # noqa: E402
    load_sst_modules,
    get_sst_genes,
    get_module_for_gene,
)

# --- Load SST gene modules from shared config (single source of truth) ---
SST_MODULES = load_sst_modules()


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Gold standard genes & druggability — from mechanism card
# ---------------------------------------------------------------------------
GOLD_STANDARD: List[str] = [
    "SGMS1", "SMPD1", "SMPD3", "PHGDH", "PSAT1", "HAVCR2",
    "ADAM10", "ADAM17", "TIGIT", "MICA", "MICB",
]

DRUGGABILITY: Dict[str, Tuple[str, str]] = {
    "HAVCR2": ("Checkpoint inhibitor", "FDA approved (lung, liver)"),
    "TIGIT": ("Checkpoint inhibitor", "Phase 3 clinical trials"),
    "ADAM10": ("Protease inhibitor", "Phase 1/2 clinical"),
    "ADAM17": ("Protease inhibitor", "Phase 1/2 clinical"),
    "SMPD1": ("Sphingomyelinase inhibitor", "Preclinical (Zheng 2023)"),
    "SMPD3": ("Sphingomyelinase inhibitor", "Preclinical (Zheng 2023)"),
    "SGMS1": ("SM synthase", "Preclinical"),
    "SGMS2": ("SM synthase", "Preclinical"),
    "PHGDH": ("Serine metabolism", "Phase 1/2 (breast cancer)"),
    "PSAT1": ("Serine metabolism", "Preclinical"),
    "PSPH": ("Serine metabolism", "Preclinical"),
    "MICA": ("Stress ligand", "Preclinical"),
    "MICB": ("Stress ligand", "Preclinical"),
    "ULBP1": ("Stress ligand", "Preclinical"),
    "ULBP2": ("Stress ligand", "Preclinical"),
    "NT5E": ("Adenosine (CD73)", "Phase 1/2 clinical"),
    "ENTPD1": ("Adenosine (CD39)", "Phase 1/2 clinical"),
    "ADORA2A": ("Adenosine receptor", "Phase 2 clinical"),
    "TGFB1": ("TGF-beta", "FDA approved (combination)"),
    "LDHA": ("Metabolic", "Preclinical"),
}

# Category classification for seed candidates
SEED_CANDIDATE_CATEGORIES: Dict[str, List[str]] = {
    "metabolic_suppression": ["LDHA", "SLC16A3", "SLC16A1", "CA9"],
    "adenosine_pathway": ["NT5E", "ENTPD1", "ADORA2A"],
    "nk_inhibitory_ligand": ["PVR", "NECTIN2", "HLA-E"],
    "stress_ligand_shedding": ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "ADAM10", "ADAM17"],
    "caf_ecm_exclusion": ["TGFB1", "TGFBR1", "TGFBR2", "COL1A1", "COL1A2", "FN1", "FAP"],
    "gastric_cancer_target": ["CLDN18", "ERBB2", "FGFR2", "MET"],
}


# ---------------------------------------------------------------------------
# Evidence computation
# ---------------------------------------------------------------------------

def compute_tumor_specificity(
    expr: pd.DataFrame,
    is_malignant: pd.Series,
    gene: str,
) -> float:
    """Compute log2 fold-change: malignant vs non-malignant.

    Returns 0.0 if gene is not in expression matrix.
    """
    if gene not in expr.columns:
        return 0.0
    mal = np.log1p(expr.loc[is_malignant, gene]).mean()
    non = np.log1p(expr.loc[~is_malignant, gene]).mean()
    return float(mal - non) if non > 0 else 0.0


def compute_nk_correlations(
    expr_full: pd.DataFrame,
    is_nk: pd.Series,
    gene_list: List[str],
) -> Dict[str, float]:
    """Correlate each gene with NK dysfunction score across NK cells.

    Returns a dict mapping gene -> Pearson correlation with dysfunction score.
    Genes not found in expression matrix get 0.0.
    """
    nk_expr = expr_full.loc[is_nk]

    # NK dysfunction score definition (Phase 5)
    dysf_genes = ["KLRC1", "TIGIT", "CD96", "HAVCR2", "TOX", "ENTPD1"]
    cyto_genes = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "XCL1", "XCL2", "CCL5"]

    avail_dysf = [g for g in dysf_genes if g in nk_expr.columns]
    avail_cyto = [g for g in cyto_genes if g in nk_expr.columns]

    dysf = nk_expr[avail_dysf].mean(axis=1) if avail_dysf else pd.Series(0.0, index=nk_expr.index)
    cyto = nk_expr[avail_cyto].mean(axis=1) if avail_cyto else pd.Series(0.0, index=nk_expr.index)
    dysf_score: pd.Series = dysf - cyto

    results: Dict[str, float] = {}
    for gene in gene_list:
        if gene not in nk_expr.columns:
            results[gene] = 0.0
        else:
            results[gene] = float(nk_expr[gene].corr(dysf_score))
    return results


def classify_candidate(
    gene: str,
    sst_module: Optional[str],
    modules: Optional[Dict[str, Dict]] = None,
) -> str:
    """Return the target category label for a gene.

    Checks seed candidate categories first, then falls back to SST module name.
    """
    for category, genes in SEED_CANDIDATE_CATEGORIES.items():
        if gene in genes:
            return category

    if sst_module:
        return f"sst_axis_{sst_module.lower().replace(' ', '_')}"

    return "unknown_candidate"


def check_direction_consistency(
    sst_module: Optional[str],
    dysf_corr: float,
) -> Optional[bool]:
    """Check if gene's dysfunction correlation matches its SST module's expected direction.

    Returns:
        True: consistent with expected direction
        False: opposite of expected direction
        None: no directional prediction available
    """
    if sst_module is None:
        return None

    if sst_module == "nk_sm_catabolism":
        return dysf_corr > 0  # more catabolism → more dysfunction
    elif sst_module in [
        "nk_sm_synthesis",
        "nk_protrusion_machinery",
        "nk_synapse_cytotoxicity_outcome",
    ]:
        return dysf_corr < 0  # more of these → less dysfunction
    else:
        return None


def is_axis_core(gene: str, modules: Optional[Dict[str, Dict]] = None) -> bool:
    """Check if gene belongs to core SST axis modules (protrusion→cytotoxicity)."""
    if modules is None:
        modules = load_sst_modules()

    core_modules = [
        "nk_protrusion_machinery",
        "nk_synapse_cytotoxicity_outcome",
        "nk_sm_synthesis",
        "nk_sm_catabolism",
    ]
    for mod_name in core_modules:
        if gene in modules.get(mod_name, {}).get("genes", []):
            return True
    return False


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main() -> None:
    log("=" * 60)
    log("CANDIDATE TARGET PRIORITIZATION (Phase 13)")
    log("=" * 60)

    out_dir = "results/tables"
    integrated_path = "data/processed/scrna/gc_integrated.h5ad"
    os.makedirs(out_dir, exist_ok=True)

    # ---- Load scRNA data ----
    log("\nLoading integrated scRNA data...")
    import scanpy as sc
    import scipy.sparse

    adata = sc.read(integrated_path)
    if scipy.sparse.issparse(adata.X):
        expr = pd.DataFrame(
            adata.X.toarray(), index=adata.obs_names, columns=adata.var_names
        )
    else:
        expr = pd.DataFrame(
            adata.X, index=adata.obs_names, columns=adata.var_names
        )

    is_nk = adata.obs["cell_type"] == "NK"
    is_other_tumor = (adata.obs["cell_type"] == "Other") & (adata.obs["condition"] == "tumor")
    epi_genes = [g for g in ["EPCAM", "KRT19", "KRT18", "KRT8", "CDH1"] if g in expr.columns]
    is_epithelial = pd.Series(False, index=expr.index)
    if epi_genes:
        is_epithelial = expr[epi_genes].mean(axis=1) > expr[epi_genes].mean(axis=1).quantile(0.75)
    is_malignant = is_other_tumor | is_epithelial
    log(f"  NK cells: {is_nk.sum()}, Malignant: {is_malignant.sum()}")

    # ---- Build candidate pool ----
    log("\nBuilding candidate pool...")

    # SST genes from shared config
    sst_genes: Set[str] = get_sst_genes(SST_MODULES)

    # Seed candidates from design doc
    seed_candidates: List[str] = [
        "LDHA", "SLC16A3", "SLC16A1", "CA9",
        "NT5E", "ENTPD1", "ADORA2A",
        "PVR", "NECTIN2", "HLA-E",
        "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "ADAM10", "ADAM17",
        "TGFB1", "TGFBR1", "TGFBR2", "COL1A1", "COL1A2", "FN1", "FAP",
        "CLDN18", "ERBB2", "FGFR2", "MET",
    ]
    candidate_pool: List[str] = list(sst_genes | set(seed_candidates))
    candidate_pool = [g for g in candidate_pool if g in expr.columns]
    log(f"  Total candidates: {len(candidate_pool)}")

    # ---- Compute evidence dimensions ----
    log("Computing evidence scores...")

    log("  1/5 Tumor specificity...")
    tumor_spec = {g: compute_tumor_specificity(expr, is_malignant, g) for g in candidate_pool}

    log("  2/5 NK dysfunction correlation...")
    nk_corr = compute_nk_correlations(expr, is_nk, candidate_pool)

    log("  3/5 NK cytotoxicity correlation...")
    nk_nk = expr.loc[is_nk]
    cyto_genes = ["NKG7", "GNLY", "GZMB", "PRF1", "IFNG"]
    avail_cyto = [g for g in cyto_genes if g in nk_nk.columns]
    cyto_score = (
        nk_nk[avail_cyto].mean(axis=1)
        if avail_cyto
        else pd.Series(0.0, index=nk_nk.index)
    )
    nk_cytotox_corr: Dict[str, float] = {}
    for g in candidate_pool:
        nk_cytotox_corr[g] = float(nk_nk[g].corr(cyto_score)) if g in nk_nk.columns else 0.0

    log("  4/5 SST axis membership...")
    gene_to_module: Dict[str, Optional[str]] = {}
    for g in candidate_pool:
        gene_to_module[g] = get_module_for_gene(g, SST_MODULES)

    log("  5/5 Literature & druggability...")

    # ---- Assemble evidence matrix ----
    log("\nAssembling evidence matrix...")
    rows: List[Dict[str, Any]] = []
    for g in candidate_pool:
        sst_module = gene_to_module.get(g)
        sst_member = 1 if sst_module else 0
        dysf_corr = nk_corr.get(g, 0.0)
        dir_consistent = check_direction_consistency(sst_module, dysf_corr)
        drug_info = DRUGGABILITY.get(g, ("", ""))
        drug_category, drug_stage = drug_info
        is_gold = 1 if g in GOLD_STANDARD else 0
        category = classify_candidate(g, sst_module)

        rows.append({
            "gene": g,
            "target_category": category,
            "sst_axis_membership": sst_module or "",
            "in_sst_axis": sst_member,
            "in_axis_core": 1 if is_axis_core(g, SST_MODULES) else 0,
            "tumor_specificity_log2": round(tumor_spec.get(g, 0.0), 4),
            "nk_dysfunction_correlation": round(dysf_corr, 4),
            "nk_cytotoxicity_correlation": round(nk_cytotox_corr.get(g, 0.0), 4),
            "axis_direction_consistent": (
                dir_consistent if dir_consistent is not None else ""
            ),
            "druggability_category": drug_category,
            "druggability_stage": drug_stage,
            "gold_standard": is_gold,
        })

    ev_df = pd.DataFrame(rows)

    # ---- Compute composite target score ----
    log("Computing composite target scores...")
    ts = ev_df["tumor_specificity_log2"].abs()
    ev_df["tumor_specificity_score"] = (
        (ts - ts.min()) / max(ts.max() - ts.min(), 1e-10)
    )

    ev_df["target_score"] = (
        0.30 * ev_df["tumor_specificity_score"]
        + 0.20 * ev_df["nk_dysfunction_correlation"].abs()
        + 0.30 * ev_df["in_sst_axis"]
        + 0.10 * ev_df["in_axis_core"]
        + 0.10 * ev_df["gold_standard"]
    )

    # Rank
    ev_df = ev_df.sort_values("target_score", ascending=False).reset_index(drop=True)
    ev_df["rank"] = range(1, len(ev_df) + 1)

    # ---- Output ----
    log("\nTop 30 candidates:")
    top_cols = [
        "rank", "gene", "target_category", "target_score",
        "tumor_specificity_log2", "nk_dysfunction_correlation",
        "axis_direction_consistent", "druggability_stage",
    ]
    print(ev_df[top_cols].head(30).to_string(index=False))

    ev_path = os.path.join(out_dir, "candidate_evidence_matrix.tsv")
    ev_df.to_csv(ev_path, sep="\t", index=False)
    log(f"\nSaved evidence matrix: {ev_path}")

    top_path = os.path.join(out_dir, "top_candidate_targets.tsv")
    ev_df.head(50).to_csv(top_path, sep="\t", index=False)
    log(f"Saved top 50: {top_path}")

    # ---- Summary ----
    log("\n--- PRIORITIZATION SUMMARY ---")
    log(f"  Total candidates evaluated: {len(ev_df)}")
    log(f"  SST-axis genes: {int(ev_df['in_sst_axis'].sum())}")
    log(f"  Top candidate: {ev_df.iloc[0]['gene']} (score={ev_df.iloc[0]['target_score']:.4f})")
    gold_in_top = int(ev_df.head(20)["gold_standard"].sum())
    log(f"  Gold standards in top 20: {gold_in_top}")

    # ---- Recommended validation assays ----
    log("\n--- RECOMMENDED ASSAYS (Top 10) ---")
    assay_map = {
        "metabolic_suppression": "NK-tumor co-culture + serine deprivation",
        "adenosine_pathway": "Co-culture with A2AR inhibitor",
        "stress_ligand_shedding": "ADAM activity assay + NK cytotoxicity",
        "caf_ecm_exclusion": "TGFb blockade + NK infiltration assay",
        "gastric_cancer_target": "qPCR/IHC in GC tissue + NK co-culture validation",
    }
    for _, r in ev_df.head(10).iterrows():
        cat = str(r["target_category"])
        assay = None
        for prefix, a in assay_map.items():
            if cat.startswith(prefix):
                assay = a
                break
        if assay is None:
            if "serine" in cat or "metabolic" in cat:
                assay = "NK-tumor co-culture + serine deprivation"
            elif "sm_" in cat or "sphingolipid" in cat:
                assay = "Sphingomyelinase-inhibitor rescue +/- Tim3 blockade"
            elif "protrusion" in cat:
                assay = "SEM / super-resolution membrane-protrusion imaging"
            else:
                assay = "qPCR/IHC in GC tissue + NK co-culture validation"
        log(f"  {r['gene']:<8} → {assay}")

    log("\nPhase 13 CANDIDATE PRIORITIZATION COMPLETE!")


if __name__ == "__main__":
    main()
