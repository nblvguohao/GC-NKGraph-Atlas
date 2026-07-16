"""
Cross-reference tumor-intrinsic candidate targets against IntOGen driver genes
(exploratory, inspired by github.com/nblvguohao/CANOPY-Router's own IntOGen
validation step, but using the real IntOGen 2023-05-31 compendium instead of
that project's hardcoded 10-gene-per-cancer proxy list).

Purpose: this project's candidate pool is deliberately curated around the
serine-SM-topology axis plus a handful of seeded known gastric-cancer targets
(ERBB2/FGFR2/MET), not drawn at random from the transcriptome. So a
hypergeometric enrichment test against the whole genome would be testing the
wrong null (any axis-curated gene list trivially "enriches" for cancer
relevance vs a genome background) and isn't reported here. Instead this
produces a descriptive breakdown by target_category:
  - seeded "gastric_cancer_target" genes are expected to overlap IntOGen
    (they were chosen because they're known drivers) -- a positive control.
  - SST-axis / NK-ligand genes are the axis-discovery's actual contribution;
    whether they overlap IntOGen or not, both directions are informative:
    overlap = converges with driver-gene consensus; no overlap = the axis
    finds immune-evasion-specific genes a generic driver detector would miss.

Data: data/raw/intogen/Compendium_Cancer_Genes.tsv (IntOGen 2023-05-31 release,
CANCER_TYPE codes STAD and HCC -- HCC is IntOGen's code for the liver cancer
cohort matching this project's TCGA-LIHC arm).

Output: results/tables/intogen_crossref.tsv

Run: python src/interpretation/intogen_crossref.py
"""
import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir

INTOGEN_PATH = "data/raw/intogen/Compendium_Cancer_Genes.tsv"
CAND_POOL_PATH = "results/tables/candidate_evidence_matrix.tsv"
CAND_FINAL_PATH = "results/tables/tumor_intrinsic_candidates.tsv"
OUT_PATH = "results/tables/intogen_crossref.tsv"

# CANOPY-Router's own hardcoded proxy lists (scripts/external_validation_intogen.py),
# kept only as a secondary, much smaller cross-check -- NOT a substitute for the
# real compendium above.
CANOPY_ROUTER_INTOGEN_PROXY = {
    "STAD": ["TP53", "CDH1", "ARID1A", "PIK3CA", "RHOA", "KRAS", "SMAD4", "CTNNB1", "RNF43", "ERBB2"],
    "LIHC": ["TP53", "CTNNB1", "AXIN1", "TERT", "ALB", "ARID1A", "RPS6KA3", "VEGFA", "MET", "CCNE1"],
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_intogen_genes(cancer_code: str) -> set:
    df = pd.read_csv(INTOGEN_PATH, sep="\t")
    return set(df.loc[df["CANCER_TYPE"] == cancer_code, "SYMBOL"].astype(str))


def main() -> None:
    log("Loading IntOGen 2023-05-31 compendium...")
    stad_genes = load_intogen_genes("STAD")
    hcc_genes = load_intogen_genes("HCC")  # HCC = IntOGen's liver-cancer code
    union_genes = stad_genes | hcc_genes
    log(f"  IntOGen STAD driver genes: {len(stad_genes)}")
    log(f"  IntOGen HCC (liver) driver genes: {len(hcc_genes)}")
    log(f"  Union: {len(union_genes)}")

    proxy_union = set(CANOPY_ROUTER_INTOGEN_PROXY["STAD"]) | set(CANOPY_ROUTER_INTOGEN_PROXY["LIHC"])

    log("\nLoading candidate tables...")
    pool = pd.read_csv(CAND_POOL_PATH, sep="\t")
    final = pd.read_csv(CAND_FINAL_PATH, sep="\t")
    log(f"  Candidate pool: {len(pool)} genes")
    log(f"  Final ranked candidates: {len(final)} genes")

    rows = []
    for _, r in pool.iterrows():
        gene = str(r["gene"])
        category = str(r.get("target_category", ""))
        in_final = gene in set(final["gene"])
        rows.append({
            "gene": gene,
            "target_category": category,
            "in_final_37": in_final,
            "in_intogen_stad": gene in stad_genes,
            "in_intogen_hcc_liver": gene in hcc_genes,
            "in_intogen_either": gene in union_genes,
            "in_canopy_router_proxy_list": gene in proxy_union,
        })
    out = pd.DataFrame(rows).sort_values(["in_final_37", "target_category"], ascending=[False, True])
    ensure_dir("results/tables")
    out.to_csv(OUT_PATH, sep="\t", index=False)
    log(f"\nWritten {len(out)} rows to {OUT_PATH}")

    # ---- Summary by category ----
    log("\n" + "=" * 70)
    log("SUMMARY: overlap with real IntOGen driver genes, by target_category")
    log("=" * 70)
    for scope_name, df in [("Full candidate pool (n=%d)" % len(pool), out),
                            ("Final ranked 37", out[out["in_final_37"]])]:
        log(f"\n-- {scope_name} --")
        for cat, sub in df.groupby("target_category"):
            n = len(sub)
            n_hit = sub["in_intogen_either"].sum()
            hit_genes = sorted(sub.loc[sub["in_intogen_either"], "gene"])
            log(f"  {cat:35s} {n_hit}/{n} in IntOGen" + (f"  -> {hit_genes}" if hit_genes else ""))

    log("\n" + "=" * 70)
    log("DONE (exploratory cross-reference; no hypergeometric test reported --")
    log("see module docstring for why a genome-background enrichment test")
    log("would test the wrong null given how the candidate pool was curated).")
    log("=" * 70)


if __name__ == "__main__":
    main()
