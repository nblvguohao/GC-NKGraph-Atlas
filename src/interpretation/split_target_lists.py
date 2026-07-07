"""
De-circularize the candidate target list (fixes R5 circularity).

PROBLEM
-------
`prioritize_targets.py` computes `target_score` with 40% weight on axis membership
(0.30 in_sst_axis + 0.10 in_axis_core) AND uses abs(tumor_specificity_log2). NK
effector markers (NKG7, PRF1, GZMB...) are strongly DEPLETED in malignant cells
(large negative log2FC -> large abs), so they score high on both terms and
dominate the ranking. The result is a list of NK readout genes, not the promised
tumor-intrinsic druggable targets.

FIX
---
Split the single ranking into two honest tables, from the existing evidence
matrix (no re-run of scRNA needed):

  1. tumor_intrinsic_candidates.tsv
       - REQUIRE tumor_specificity_log2 > 0 (gene must be expressed by tumor cells)
       - re-score with SIGNED tumor specificity + mechanistic relevance
         (tumor-serine module or axis-druggable enzyme) + druggability + NK
         association. Axis membership is annotation, not a score driver.
  2. axis_confirmation_panel.tsv
       - the NK-side axis markers (in_sst_axis, negative tumor specificity),
         labeled explicitly as the AXIS READOUT, not as targets.

USAGE
-----
    python src/interpretation/split_target_lists.py \
        --in results/tables/candidate_evidence_matrix.tsv
"""

import argparse
import os

import numpy as np
import pandas as pd

# tumor-side serine program + druggable SM enzymes = mechanistically privileged
MECHANISTIC_TUMOR_SIDE = {
    "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1", "MTHFD1L",
    "MTHFD2", "SLC1A4", "SLC1A5",
}
AXIS_DRUGGABLE_ENZYME = {"SMPD1", "SMPD2", "SMPD3", "SMPD4", "SGMS1", "SGMS2"}

ASSAY_BY_CATEGORY = {
    "metabolic_suppression": "NK-tumor co-culture + serine-pathway inhibition (PHGDH inhibitor)",
    "adenosine_pathway": "Co-culture with CD73/A2AR inhibitor + NK cytotoxicity",
    "stress_ligand_shedding": "ADAM protease-activity assay + NK cytotoxicity",
    "caf_ecm_exclusion": "TGF-beta blockade + NK infiltration/cytotoxicity",
    "gastric_cancer_target": "IHC/qPCR in GC tissue + NK co-culture cytotoxicity",
    "sst_axis_tumor_serine_capacity": "NK-tumor co-culture + serine deprivation",
    "sst_axis_nk_sm_catabolism": "Sphingomyelinase-inhibitor rescue +/- Tim3 (HAVCR2) blockade",
    "sst_axis_nk_sm_synthesis": "SM-synthase modulation + NK membrane-SM readout",
}


def assay_for(cat):
    for k, v in ASSAY_BY_CATEGORY.items():
        if str(cat).startswith(k):
            return v
    return "IHC/qPCR in GC tissue + NK co-culture cytotoxicity"


def minmax(s):
    s = s.astype(float)
    rng = s.max() - s.min()
    return (s - s.min()) / rng if rng > 1e-12 else s * 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path",
                    default="results/tables/candidate_evidence_matrix.tsv")
    ap.add_argument("--out-dir", default="results/tables")
    args = ap.parse_args()

    ev = pd.read_csv(args.in_path, sep="\t")
    os.makedirs(args.out_dir, exist_ok=True)

    # ---------- Table 1: tumor-intrinsic candidates ----------
    ti = ev[ev["tumor_specificity_log2"] > 0].copy()

    ti["mechanistic_bonus"] = ti["gene"].apply(
        lambda g: 1.0 if g in MECHANISTIC_TUMOR_SIDE
        else (0.7 if g in AXIS_DRUGGABLE_ENZYME else 0.0))
    ti["druggable"] = ti["druggability_stage"].notna() & (ti["druggability_stage"] != "")

    # signed tumor-specificity (higher = more tumor-restricted), 0..1
    ti["tumor_spec_signed"] = minmax(ti["tumor_specificity_log2"])
    ti["nk_assoc"] = minmax(ti["nk_dysfunction_correlation"].abs())

    # transparent re-weighted score: tumor-intrinsic evidence leads
    ti["target_score_v2"] = (
        0.40 * ti["tumor_spec_signed"]
        + 0.25 * ti["mechanistic_bonus"]
        + 0.20 * ti["druggable"].astype(float)
        + 0.15 * ti["nk_assoc"]
    ).round(4)
    ti["recommended_assay"] = ti["target_category"].apply(assay_for)

    ti = ti.drop(columns=[c for c in ("rank",) if c in ti.columns])
    ti = ti.sort_values("target_score_v2", ascending=False).reset_index(drop=True)
    ti.insert(0, "rank", range(1, len(ti) + 1))
    cols = ["rank", "gene", "target_category", "target_score_v2",
            "tumor_specificity_log2", "mechanistic_bonus", "druggability_stage",
            "gold_standard", "nk_dysfunction_correlation", "recommended_assay"]
    out1 = os.path.join(args.out_dir, "tumor_intrinsic_candidates.tsv")
    ti[cols].to_csv(out1, sep="\t", index=False)

    # ---------- Table 2: axis-confirmation panel (readout, not targets) ----------
    panel = ev[(ev["in_sst_axis"] == 1) & (ev["tumor_specificity_log2"] <= 0)].copy()
    panel["role"] = "AXIS_READOUT_not_target"
    panel = panel.sort_values("nk_cytotoxicity_correlation", ascending=False)
    out2 = os.path.join(args.out_dir, "axis_confirmation_panel.tsv")
    panel[["gene", "sst_axis_membership", "role", "tumor_specificity_log2",
           "nk_cytotoxicity_correlation", "nk_dysfunction_correlation"]].to_csv(
        out2, sep="\t", index=False)

    print(f"[tumor-intrinsic candidates: {len(ti)}] -> {out1}")
    print(ti[cols].head(15).to_string(index=False))
    print(f"\n[axis-confirmation panel: {len(panel)}] -> {out2}")
    print(panel[["gene", "sst_axis_membership", "nk_cytotoxicity_correlation"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
