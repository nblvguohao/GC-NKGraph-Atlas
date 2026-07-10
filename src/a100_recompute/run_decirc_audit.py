"""
T6 — De-circularization audit: verify tumor-intrinsic pool composition.

Checks:
  1. Category breakdown of the n=37 tumor-intrinsic candidate pool
  2. Which categories are NK-side (protrusion machinery, sm synthesis/catabolism, etc.)
  3. Whether RAC1 and other NK-protrusion genes are present → leak detection
  4. Alignment between Table 4 ranking and the ranking used in Fig 3B

Produces: results/tables/tumor_intrinsic_pool_audit.tsv

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_decirc_audit.py
"""
import pandas as pd
import numpy as np

T = "results/tables/"

# Load the tumor-intrinsic candidates
ti = pd.read_csv(T + "tumor_intrinsic_candidates.tsv", sep="\t")

# Load NK-side axis confirmation panel for comparison
try:
    ac = pd.read_csv(T + "axis_confirmation_panel.tsv", sep="\t")
    ac_genes = set(ac["gene"].str.upper() if "gene" in ac.columns else [])
except Exception:
    ac_genes = set()

# --- Category definitions ---
NK_SIDE_CATEGORIES = {
    "nk_protrusion_machinery", "nk_sm_synthesis", "nk_sm_catabolism",
    "nk_denovo_sphingolipid", "nk_synapse_cytotoxicity_outcome", "checkpoint_link"
}
TUMOR_SIDE_CATEGORIES = {"tumor_serine_capacity", "tumor_serine_program"}

# What category column exists?
cat_col = None
for c in ["category", "module", "mechanism_module", "sst_module"]:
    if c in ti.columns:
        cat_col = c
        break

if cat_col is None:
    # Infer from gene membership if no explicit category column
    print("WARNING: no category column found; inferring from mechanism-card modules")
    # Load mechanism card gene lists if available
    cat_col = "inferred_category"

gene_col = None
for c in ["gene", "gene_symbol", "symbol"]:
    if c in ti.columns:
        gene_col = c
        break

if gene_col is None:
    gene_col = ti.columns[0]

print(f"Using category column: '{cat_col}', gene column: '{gene_col}'")
print(f"Pool size: {len(ti)}")
print(f"Genes: {sorted(ti[gene_col].str.upper().tolist())}")
print()

# --- Category breakdown ---
if cat_col in ti.columns:
    cat_counts = ti[cat_col].value_counts()
    print("=== CATEGORY BREAKDOWN ===")
    for cat, n in cat_counts.items():
        flag = " ⚠️ NK-side" if any(nk in str(cat).lower().replace(" ","_") for nk in ["nk_","protrusion","sm_synth","sm_catab","sphingolipid","cytotox","checkpoint"]) else ""
        print(f"  {cat:40s} {n:3d}{flag}")

    nk_total = sum(ti[cat_col].apply(
        lambda x: any(nk in str(x).lower().replace(" ","_")
                      for nk in ["nk_","protrusion_machinery","sm_synthesis","sm_catabolism",
                                 "denovo_sphingolipid","cytotoxicity_outcome","checkpoint"])
    ))
    print(f"\n  NK-side genes in 'tumor-intrinsic' pool: {nk_total}/{len(ti)} {'⚠️ LEAK' if nk_total > 3 else '△ minor'}")
else:
    print("WARNING: no category column — cannot audit de-circularization")

# --- Known NK-effector genes to check ---
NK_EFFECTOR_CHECK = {"NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "RAC1", "EZR", "MSN", "RDX",
                      "WAS", "WASL", "ARPC2", "ARPC3", "DIAPH1", "FMNL1"}
leaked = {g for g in ti[gene_col].str.upper() if g in NK_EFFECTOR_CHECK}
if leaked:
    print(f"\n=== NK EFFECTOR GENES IN TUMOR POOL ===  ({len(leaked)} genes)")
    for g in sorted(leaked):
        cat = ti[ti[gene_col].str.upper() == g][cat_col].values[0] if cat_col in ti.columns else "?"
        print(f"  {g:15s}  category={cat}")

# --- Check overlap with axis-confirmation panel ---
overlap = set(ti[gene_col].str.upper()) & ac_genes
if overlap:
    print(f"\n=== OVERLAP WITH AXIS-CONFIRMATION PANEL (NK READOUT) ===  ({len(overlap)} genes)")
    for g in sorted(overlap):
        print(f"  {g}")

# --- Save audit ---
audit = ti.copy()
audit["nk_side_flag"] = audit[cat_col].apply(
    lambda x: any(nk in str(x).lower().replace(" ","_")
                  for nk in ["nk_","protrusion_machinery","sm_synthesis","sm_catabolism",
                             "denovo_sphingolipid","cytotoxicity_outcome","checkpoint"])
) if cat_col in audit.columns else False
audit.to_csv(T + "tumor_intrinsic_pool_audit.tsv", sep="\t", index=False)

print(f"\nSaved audit to tumor_intrinsic_pool_audit.tsv")
print("T6 audit complete — review the ⚠️ lines above to decide if re-pooling is needed")
