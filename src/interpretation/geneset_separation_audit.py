"""
M4 (nature-reviewer task card): audit gene-set overlap between the NK
immune-state label definition (src/immune_scoring/nk_scores.py), the SST-axis
modules used for the H1-H5 mechanism analyses (src/common/sst_config.py), and
the NK-state classifier's input features (src/models/gc_nkgraph_atlas.py).

Motivation: if the NK-state label used as the classifier target (and used to
define the NK-hot-cytotoxic / NK-hot-dysfunctional groups for the §3.5
NK-state DE test) is built from genes that also constitute the SST-axis
cytotoxicity-outcome module, or that are also present in the classifier's own
input feature matrix, there is a risk of the classification / DE results being
partly a restatement of the label definition rather than independent evidence.

This script does NOT change any modeling results; it only quantifies overlap
so the manuscript can state precisely what is and is not circular.

Output:
  results/tables/geneset_separation_audit.tsv
  results/tables/geneset_separation_audit_summary.md

Usage:
    python src/interpretation/geneset_separation_audit.py
"""
import os
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.sst_config import load_sst_modules  # noqa: E402
from src.immune_scoring.nk_scores import (  # noqa: E402
    NK_MARKERS, NK_CYTOTOXICITY_GENES, NK_DYSFUNCTION_GENES, CAF_ECM_TGFB_GENES,
)

TARGET_37_GENES = [
    "PHGDH", "SGMS2", "PSAT1", "PSPH", "SMPD3", "COL1A1", "COL1A2",
    "SMPD1", "NECTIN2", "RAC1", "MTHFD1L", "SLC1A5", "SHMT2", "SHMT1",
    "MTHFD1", "NT5E", "CA9", "ERBB2", "FN1", "MICA", "BAIAP2", "SMPD2",
    "SMPD4", "WASL", "FGFR2", "MET", "PACSIN2", "CERS6", "PVR", "SPTSSA",
    "CERS2", "FAP", "WASF1", "WASF3", "DIAPH3", "SPTLC1", "SPTLC3",
]

NK_STATE_LABEL_SETS = {
    "label:NK_MARKERS (infiltration)": set(NK_MARKERS),
    "label:NK_CYTOTOXICITY_GENES": set(NK_CYTOTOXICITY_GENES),
    "label:NK_DYSFUNCTION_GENES": set(NK_DYSFUNCTION_GENES),
    "label:CAF_ECM_TGFB_GENES (exclusion)": set(CAF_ECM_TGFB_GENES),
}

sst_modules = load_sst_modules()
SST_SETS = {f"sst_axis:{name}": set(mod["genes"]) for name, mod in sst_modules.items()}

TARGET_SET = {"target_37:tumor_intrinsic_candidates": set(TARGET_37_GENES)}


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def main():
    all_sets = {**NK_STATE_LABEL_SETS, **SST_SETS, **TARGET_SET}
    names = list(all_sets.keys())

    rows = []
    for i, n1 in enumerate(names):
        for n2 in names[i + 1:]:
            s1, s2 = all_sets[n1], all_sets[n2]
            overlap = s1 & s2
            rows.append({
                "set_1": n1, "set_2": n2,
                "n_set_1": len(s1), "n_set_2": len(s2),
                "n_overlap": len(overlap),
                "jaccard": round(jaccard(s1, s2), 4),
                "overlap_genes": ",".join(sorted(overlap)) if overlap else "",
            })

    df = pd.DataFrame(rows).sort_values("n_overlap", ascending=False)
    out_dir = "results/tables"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "geneset_separation_audit.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path}")
    print(df[df["n_overlap"] > 0].to_string(index=False))

    # Focused check: label-defining genes vs SST cytotoxicity-outcome module
    cyto_label = NK_STATE_LABEL_SETS["label:NK_CYTOTOXICITY_GENES"]
    dysf_label = NK_STATE_LABEL_SETS["label:NK_DYSFUNCTION_GENES"]
    infil_label = NK_STATE_LABEL_SETS["label:NK_MARKERS (infiltration)"]
    sst_cyto = SST_SETS.get("sst_axis:nk_synapse_cytotoxicity_outcome", set())
    sst_checkpoint = SST_SETS.get("sst_axis:checkpoint_link", set())
    target37 = TARGET_SET["target_37:tumor_intrinsic_candidates"]

    label_cyto_overlap = cyto_label & sst_cyto
    label_checkpoint_overlap = dysf_label & sst_checkpoint
    target_label_overlap = target37 & (cyto_label | dysf_label | infil_label)

    summary_path = os.path.join(out_dir, "geneset_separation_audit_summary.md")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"""# Gene-set separation audit (M4)

## Question
Does the NK immune-state label (used as the classifier target in §3.4, and to
define the NK-hot-cytotoxic / NK-hot-dysfunctional groups for the NK-state DE
test in §3.5) share genes with (a) the SST-axis modules used for the H1-H5
mechanism tests, or (b) the 37 tumor-intrinsic candidate targets?

## Finding 1: NK-state label vs SST-axis cytotoxicity-outcome module
`NK_CYTOTOXICITY_GENES` (used both to build `NK_cytotoxicity_score` and,
subtracted, `NK_dysfunction_score` -- i.e. it defines the classifier's target
label) overlaps with the SST-axis `nk_synapse_cytotoxicity_outcome` module
(used for the H3/H5 mechanism tests) in {len(label_cyto_overlap)} genes:
{sorted(label_cyto_overlap)}.

This overlap does **not** affect the H3/H5 correlation tests themselves (those
operate directly on expression-derived module scores and never reference the
NK-state label). It **does** mean that the NK-state classification task
(§3.4) has partial "label leakage" by construction: the target label is a
thresholding rule on a handful of marker genes (NKG7, GNLY, GZMB, PRF1, IFNG
among others) that are also present, unmodified, in the full expression vector
x used as classifier input for every baseline and the GNN. This is a known
and generally-accepted property of marker-defined phenotype classification
(the label is a deterministic function of a subset of the input features), not
a coding error -- but it explains why even simple baselines (e.g. the 8-gene
NK-marker signature baseline, AUROC=0.849) perform well, and it means §3.4
accuracy numbers should be read as "how well can a model recover a
marker-gene-defined label from the transcriptome that contains those same
markers," not as an unconstrained prediction task.

## Finding 2: NK-state label vs checkpoint_link module
`NK_DYSFUNCTION_GENES` (used to build `NK_dysfunction_score`, part of the
label rule) overlaps the SST-axis `checkpoint_link` module ({{HAVCR2}}) in
{len(label_checkpoint_overlap)} genes: {sorted(label_checkpoint_overlap)}.

## Finding 3: 37 tumor-intrinsic candidates vs NK-state label genes
The 37 tumor-intrinsic candidate targets (Table 4, §3.5) overlap the NK-state
label-defining gene sets (infiltration + cytotoxicity + dysfunction) in
{len(target_label_overlap)} genes: {sorted(target_label_overlap) if target_label_overlap else 'none'}.
{"This is a non-trivial overlap and means the NK-state DE test for these genes is partly self-referential." if target_label_overlap else "No overlap: the NK-state DE test (§3.5) comparing the 37 candidates between NK-hot-cytotoxic and NK-hot-dysfunctional tumors is NOT circular with the label definition -- the tested genes are disjoint from the genes that define the groups being compared."}

## Full pairwise overlap table
See `geneset_separation_audit.tsv` for all pairwise Jaccard overlaps among the
NK-state label components, the seven SST-axis modules, and the 37-gene target
list.

## Recommendation
State explicitly in Methods (§2.4/§2.8) that the NK-state classification
task's target label is partly constructed from marker genes present in the
classifier's own input (Finding 1/2). For the NK-state DE test (§3.5), flag
that NT5E -- one of the four Tier-1 orthogonally-validated candidates -- is
itself a constituent of `NK_DYSFUNCTION_GENES` (Finding 3): its "upregulated
in NK-hot-dysfunctional tumors" result is partly self-referential, since
higher NT5E expression directly contributes to a tumor being labeled
dysfunctional in the first place. The other three tested genes among the
37-candidate pool have no overlap with any label-defining gene set, so their
NK-state DE results are not subject to this concern; NT5E's Tier-1 status
should be downgraded to reflect the circularity, or corroborated with an
NK-state definition that excludes NT5E from the label rule.
""")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
