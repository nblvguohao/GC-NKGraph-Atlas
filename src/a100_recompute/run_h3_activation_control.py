"""
T14 — H3 activation-control: test whether protrusion~cytotoxicity survives
after removing a generic NK-activation component.

Rationale (from peer review):
  The "cytotoxicity-output" module (NKG7/GNLY/GZMB/PRF1/IFNG + LAT/VAV1/
  TLN1/ITGAL/ITGB2/LCP2) and the "protrusion-machinery" module share an
  NK-activation/synapse program. A positive H3 could just be co-activation
  rather than a specific recovery of the protrusion→cytotoxicity link.

Approach:
  1. Define a generic NK-activation signature from literature genes
     (e.g. CD69, TNF, XCL1, XCL2, CCL3, CCL4, CCL5, CSF2, IFNG — classic
     NK activation markers NOT in protrusion or cytotoxicity modules)
  2. Compute per-cell activation score in scNK data
  3. Compute partial correlation: protrusion ~ cytotoxicity | activation
  4. Report: raw r, partial r, and variance explained by activation vs residual

Produces: results/tables/h3_activation_control.tsv

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_h3_activation_control.py
"""
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

T = "results/tables/"

# --- 1. Define a generic NK-activation signature ---
# Classic NK activation markers NOT overlapping with protrusion or cytotoxicity modules.
# Sources: Vivier 2008, Chiossone 2018.
NK_ACTIVATION_GENES = [
    "CD69",      # early activation marker
    "TNF",       # cytokine, not in cytotoxicity module
    "XCL1", "XCL2",  # chemokines
    "CCL3", "CCL4", "CCL5",  # chemokines
    "CSF2",      # GM-CSF
    "IL2RA",     # CD25
    "ICOS",      # co-stimulatory
    "TNFSF10",  # TRAIL
    "FASLG",     # Fas ligand
    "CD38",      # activation ectoenzyme
    "HLA-DRA", "HLA-DRB1",  # MHC-II (activation)
    "MKI67",     # proliferation
]

# --- 2. Load scRNA data ---
sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")
if "cell_type" in sc.columns:
    nk = sc[sc["cell_type"] == "NK"].copy()
else:
    # Table already NK-only (produced by sst_axis.py on a pre-subset NK h5ad)
    nk = sc.copy()
print(f"scNK cells: {len(nk)}")

# Identify columns for protrusion and cytotoxicity scores
prot_col = None
cyto_col = None
for c in nk.columns:
    if "protrusion" in c.lower() and "score" in c.lower():
        prot_col = c
    if ("cytotox" in c.lower() or "synapse" in c.lower()) and "score" in c.lower():
        cyto_col = c

if prot_col is None:
    # Try alternative column names
    for c in nk.columns:
        if "protrusion" in c.lower():
            prot_col = c
            break
if cyto_col is None:
    for c in nk.columns:
        if "cytotox" in c.lower():
            cyto_col = c
            break

print(f"Protrusion column: {prot_col}")
print(f"Cytotoxicity column: {cyto_col}")

if prot_col is None or cyto_col is None:
    print("ERROR: cannot identify protrusion/cytotoxicity columns")
    print("Available columns:", list(nk.columns[:30]))
    import sys; sys.exit(1)

# --- 3. Compute generic activation score ---
# Find which activation genes are in the dataset
act_genes_present = [g for g in NK_ACTIVATION_GENES if g in nk.columns]
print(f"Activation genes found: {len(act_genes_present)}/{len(NK_ACTIVATION_GENES)}")
print(f"  {act_genes_present}")

if len(act_genes_present) < 3:
    # Fallback: use IFNG + TNF + CCL4 as minimal activation signature
    act_genes_present = [g for g in ["IFNG", "TNF", "CCL4", "CCL5", "XCL1"]
                         if g in nk.columns]
    print(f"Fallback activation genes: {act_genes_present}")

if len(act_genes_present) < 2:
    print("ERROR: too few activation genes available")
    import sys; sys.exit(1)

# Compute mean z-score activation score
act_expr = nk[act_genes_present].copy()
# Z-score per gene
act_z = (act_expr - act_expr.mean()) / (act_expr.std() + 1e-10)
activation_score = act_z.mean(axis=1)
nk["_activation_score"] = activation_score

# --- 4. Compute correlations ---
x = nk[prot_col].values
y = nk[cyto_col].values
a = activation_score.values

valid = ~(np.isnan(x) | np.isnan(y) | np.isnan(a))
x, y, a = x[valid], y[valid], a[valid]
n_valid = len(x)
print(f"\nValid cells: {n_valid}")

# Raw H3
r_raw, p_raw = stats.pearsonr(x, y)
print(f"\n=== RAW H3 ===")
print(f"protrusion ~ cytotoxicity:  r = {r_raw:.4f}  p = {p_raw:.2e}  r² = {r_raw**2:.4f}")

# Activation ~ protrusion
r_ap, p_ap = stats.pearsonr(a, x)
print(f"activation ~ protrusion:    r = {r_ap:.4f}  p = {p_ap:.2e}")

# Activation ~ cytotoxicity
r_ac, p_ac = stats.pearsonr(a, y)
print(f"activation ~ cytotoxicity:  r = {r_ac:.4f}  p = {p_ac:.2e}")

# Partial correlation: protrusion ~ cytotoxicity | activation
# Residualize both on activation, then correlate residuals
from scipy.stats import linregress
# protrusion ~ activation
slope_xa, intercept_xa, _, _, _ = linregress(a, x)
x_resid = x - (slope_xa * a + intercept_xa)
# cytotoxicity ~ activation
slope_ya, intercept_ya, _, _, _ = linregress(a, y)
y_resid = y - (slope_ya * a + intercept_ya)

r_partial, p_partial = stats.pearsonr(x_resid, y_resid)
print(f"\n=== PARTIAL (|activation) ===")
print(f"protrusion ~ cytotoxicity | activation:  r_partial = {r_partial:.4f}  p = {p_partial:.2e}  r² = {r_partial**2:.4f}")

# --- 5. Same for bulk ---
print(f"\n=== BULK CONTROL ===")
try:
    liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
    bulk_prot = liver.get("nk_protrusion_machinery_score", None)
    bulk_cyto = liver.get("nk_synapse_cytotoxicity_outcome_score", None)
    if bulk_prot is not None and bulk_cyto is not None:
        r_b, p_b = stats.pearsonr(bulk_prot.dropna(), bulk_cyto.dropna())
        print(f"Bulk H3 raw: r={r_b:.4f} p={p_b:.2e}")
        # Activation genes in bulk
        bulk_act_genes = [g for g in act_genes_present if g in liver.columns]
        if len(bulk_act_genes) >= 2:
            idx = bulk_prot.dropna().index.intersection(bulk_cyto.dropna().index)
            bulk_act = liver.loc[idx, bulk_act_genes]
            bulk_act_z = (bulk_act - bulk_act.mean()) / (bulk_act.std() + 1e-10)
            bulk_act_score = bulk_act_z.mean(axis=1)
            bp, bc, ba = bulk_prot.loc[idx].values, bulk_cyto.loc[idx].values, bulk_act_score.values
            s_xa, _, _, _, _ = linregress(ba, bp)
            s_ya, _, _, _, _ = linregress(ba, bc)
            r_bp, p_bp = stats.pearsonr(bp - (s_xa*ba + 0), bc - (s_ya*ba + 0))
            print(f"Bulk H3 partial|activation: r={r_bp:.4f} p={p_bp:.2e}")
except Exception as e:
    print(f"Bulk check skipped: {e}")

# --- 6. Verdict ---
print(f"\n{'='*60}")
r2_raw = r_raw**2
r2_partial = r_partial**2
r2_explained_by_activation = r2_raw - r2_partial

print(f"Variance explained: raw={r2_raw:.4f}  partial={r2_partial:.4f}  Δ(activation)={r2_explained_by_activation:.4f}")

if r_partial > 0.15 and p_partial < 0.001:
    print("VERDICT: H3 ROBUST to activation control")
    print("→ protrusion→cytotoxicity is NOT just co-activation")
    print("→ the 'effector arm recovers' claim survives")
elif r_partial > 0.05 and p_partial < 0.05:
    print("VERDICT: H3 PARTIALLY survives activation control")
    print("→ protrusion→cytotoxicity has a residual relationship after controlling for activation,")
    print("  but a substantial portion is shared with the activation program")
    print("→ recommend: write 'partially independent of a generic NK-activation signature'")
else:
    print("VERDICT: H3 does NOT survive activation control")
    print("→ protrusion~cytotoxicity is largely explained by co-activation")
    print("→ the 'effector arm recovers' claim should be DOWNGRADED to 'protrusion-machinery and cytotoxicity")
    print("  transcripts co-vary with NK activation state, consistent with but not specific to the SST axis'")

# --- Save ---
results = {
    "resolution": ["single-cell NK", "single-cell NK", "single-cell NK"],
    "test": ["raw protrusion~cytotox", "partial(|activation)", "activation_component"],
    "r": [r_raw, r_partial, np.sqrt(r2_explained_by_activation)],
    "r2": [r2_raw, r2_partial, r2_explained_by_activation],
    "p": [p_raw, p_partial, np.nan],
    "n": [n_valid, n_valid, n_valid],
    "activation_genes_used": [len(act_genes_present)] * 3,
}
pd.DataFrame(results).to_csv(T + "h3_activation_control.tsv", sep="\t", index=False)
print(f"\nSaved to {T}h3_activation_control.tsv")
print("T14 PASS")
