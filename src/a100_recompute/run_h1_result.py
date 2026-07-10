"""
T5 — H1 result: tumor_serine_capacity ⟂ nk_sm_balance

Compute the missing H1 hypothesis result at bulk (TCGA-LIHC) and single-cell NK
resolution. Writes a one-row-per-resolution table appended to the existing
recovery file and prints the H1 outcome.

Run:  conda activate gc-nkgraph && python src/a100_recompute/run_h1_result.py
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

T = "results/tables/"
liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
sc     = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")

# --- Bulk H1 ---
x_bulk = liver["tumor_serine_capacity_score"]
y_bulk = liver["nk_sm_balance_score"] if "nk_sm_balance_score" in liver.columns else (
    liver["nk_sm_synthesis_score"] - liver["nk_sm_catabolism_score"]
)
r_b, p_b = stats.pearsonr(x_bulk.dropna(), y_bulk.dropna())

# --- scRNA H1 ---
nk = sc[sc["cell_type"] == "NK"].copy()
x_sc = nk["tumor_serine_capacity_score"] if "tumor_serine_capacity_score" in nk.columns else nk["tumor_serine_capacity"]
y_sc = nk["nk_sm_balance"] if "nk_sm_balance" in nk.columns else (nk["nk_sm_synthesis"] - nk["nk_sm_catabolism"])
valid = x_sc.notna() & y_sc.notna()
r_sc, p_sc = stats.pearsonr(x_sc[valid], y_sc[valid])

# --- Write results ---
rows = [
    {"hypothesis": "H1", "test": "serine_capacity ~ sm_balance", "resolution": "bulk",
     "r": round(r_b, 4), "p": p_b, "expected": "calibrated", "outcome": "reported"},
    {"hypothesis": "H1", "test": "serine_capacity ~ sm_balance", "resolution": "single-cell NK",
     "r": round(r_sc, 4), "p": p_sc, "expected": "calibrated", "outcome": "reported"},
]
df_new = pd.DataFrame(rows)
# Append to existing recovery table
rec = pd.read_csv(T + "sst_axis_positive_control_recovery.tsv", sep="\t")
rec_full = pd.concat([df_new, rec], ignore_index=True)
rec_full.to_csv(T + "sst_axis_positive_control_recovery.tsv", sep="\t", index=False)

print("=== H1 RESULT ===")
print(f"Bulk:  tumor_serine_capacity ~ nk_sm_balance  r={r_b:.4f}  p={p_b:.2e}  n={len(x_bulk.dropna())}")
print(f"scNK:  tumor_serine_capacity ~ nk_sm_balance  r={r_sc:.4f}  p={p_sc:.2e}  n={valid.sum()}")
print(f"\nAppended to sst_axis_positive_control_recovery.tsv ({len(rec_full)} rows)")
print("T5 PASS" if not np.isnan(r_b) else "T5 FAIL — check column names in source tables")
