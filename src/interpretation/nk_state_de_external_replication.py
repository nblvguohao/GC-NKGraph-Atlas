"""
S3 (nature-reviewer task card, optional external replication): replicate the
NK-state differential-expression test (Sec 3.5 / Table 4) in two INDEPENDENT
external gastric bulk cohorts that are already processed locally for the
effector-arm external validation (Sec 3.3) -- GSE62254 (ACRG, n=300) and
GSE84437 (n=483) -- rather than acquiring a new dataset.

The original TCGA-STAD test (n=134 hot-cytotoxic vs n=20 hot-dysfunctional,
underpowered per Limitations item 14) compared the 37 tumor-intrinsic
candidate genes between NK-hot-cytotoxic and NK-hot-dysfunctional tumors.
Here we independently derive NK-hot-cytotoxic / NK-hot-dysfunctional labels
in each external cohort using the same scoring rule
(src/immune_scoring/nk_scores.py), then re-run the same DE test.

Output:
  results/tables/nk_state_de_external_replication.tsv
  results/tables/nk_state_de_external_replication_summary.md

Usage:
    python src/interpretation/nk_state_de_external_replication.py
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.immune_scoring.nk_scores import compute_nk_scores, assign_immune_states  # noqa: E402
from src.common.log_utils import Logger  # noqa: E402

GENES_37 = [
    "PHGDH", "SGMS2", "PSAT1", "PSPH", "SMPD3", "COL1A1", "COL1A2",
    "SMPD1", "NECTIN2", "RAC1", "MTHFD1L", "SLC1A5", "SHMT2", "SHMT1",
    "MTHFD1", "NT5E", "CA9", "ERBB2", "FN1", "MICA", "BAIAP2", "SMPD2",
    "SMPD4", "WASL", "FGFR2", "MET", "PACSIN2", "CERS6", "PVR", "SPTSSA",
    "CERS2", "FAP", "WASF1", "WASF3", "DIAPH3", "SPTLC1", "SPTLC3",
]

# Genes flagged in the TCGA-STAD internal test (Table 4 / Sec 3.5), for
# direct directional comparison against these external cohorts.
TCGA_STAD_DIRECTION = {
    "SGMS2": "UP", "NT5E": "UP", "PSPH": "DOWN", "SHMT1": "DOWN",
    "SHMT2": "DOWN", "MTHFD1L": "DOWN", "MTHFD1": "DOWN", "RAC1": "DOWN",
}

COHORTS = {
    "GSE62254": "data/processed/bulk/gse62254_expression.tsv",
    "GSE84437": "data/processed/bulk/gse84437_expression.tsv",
}


def log(msg):
    print(msg, flush=True)


def run_cohort(name, expr_path, logger):
    log(f"\n{'='*70}\n{name}\n{'='*70}")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    log(f"  Expression: {expr.shape[0]} samples x {expr.shape[1]} genes")

    scores = compute_nk_scores(expr, logger, dataset_name=name)
    labeled, thresholds = assign_immune_states(scores, logger=logger)
    counts = labeled["nk_immune_state"].value_counts()
    log(f"  NK-state distribution: {counts.to_dict()}")

    n_cyto = int(counts.get("NK-hot-cytotoxic", 0))
    n_dysf = int(counts.get("NK-hot-dysfunctional", 0))
    if n_cyto < 5 or n_dysf < 5:
        log(f"  SKIP: too few samples in one group (cytotoxic={n_cyto}, dysfunctional={n_dysf})")
        return None

    cyto_samples = labeled[labeled["nk_immune_state"] == "NK-hot-cytotoxic"].index
    dysf_samples = labeled[labeled["nk_immune_state"] == "NK-hot-dysfunctional"].index

    rows = []
    for gene in GENES_37:
        if gene not in expr.columns:
            continue
        a = expr.loc[cyto_samples, gene].dropna()
        b = expr.loc[dysf_samples, gene].dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        log2fc = b.mean() - a.mean()  # dysfunctional - cytotoxic (log-scale already for microarray)
        t, p = stats.ttest_ind(b, a, equal_var=False)
        direction = "UP" if log2fc > 0 else "DOWN"
        rows.append({
            "cohort": name, "gene": gene, "n_cytotoxic": len(a), "n_dysfunctional": len(b),
            "log2FC_dysf_vs_cyto": round(float(log2fc), 4), "p": p,
            "direction_in_dysf": direction,
        })
    df = pd.DataFrame(rows)
    log(f"  Tested {len(df)}/{len(GENES_37)} genes with coverage")
    return df, n_cyto, n_dysf


def main():
    logger = Logger()
    all_results = []
    cohort_ns = {}
    for name, path in COHORTS.items():
        if not os.path.exists(path):
            log(f"SKIP {name}: {path} not found")
            continue
        result = run_cohort(name, path, logger)
        if result is None:
            continue
        df, n_cyto, n_dysf = result
        all_results.append(df)
        cohort_ns[name] = (n_cyto, n_dysf)

    if not all_results:
        log("No cohorts produced usable results.")
        sys.exit(1)

    combined = pd.concat(all_results, ignore_index=True)
    out_dir = "results/tables"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "nk_state_de_external_replication.tsv")
    combined.to_csv(out_path, sep="\t", index=False)
    log(f"\nSaved: {out_path}")

    # ── Directional concordance with TCGA-STAD for the flagged genes ──
    log(f"\n{'='*70}\nDIRECTIONAL CONCORDANCE WITH TCGA-STAD (internal test)\n{'='*70}")
    concord_rows = []
    for gene, tcga_dir in TCGA_STAD_DIRECTION.items():
        sub = combined[combined["gene"] == gene]
        if sub.empty:
            continue
        for _, r in sub.iterrows():
            match = "concordant" if r["direction_in_dysf"] == tcga_dir else "discordant"
            concord_rows.append({
                "gene": gene, "cohort": r["cohort"], "tcga_stad_direction": tcga_dir,
                "external_direction": r["direction_in_dysf"],
                "external_log2FC": r["log2FC_dysf_vs_cyto"], "external_p": r["p"],
                "concordance": match,
            })
            log(f"  {gene:<10} {r['cohort']:<10} TCGA-STAD={tcga_dir:<5} "
                f"external={r['direction_in_dysf']:<5} ({match})")
    concord_df = pd.DataFrame(concord_rows)
    concord_path = os.path.join(out_dir, "nk_state_de_external_concordance.tsv")
    concord_df.to_csv(concord_path, sep="\t", index=False)
    log(f"\nSaved: {concord_path}")

    n_concordant = (concord_df["concordance"] == "concordant").sum() if not concord_df.empty else 0
    n_total = len(concord_df)
    log(f"\nConcordance rate: {n_concordant}/{n_total}")

    # ── Summary ──
    summary_path = os.path.join(out_dir, "nk_state_de_external_replication_summary.md")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"""# NK-state DE external replication (S3 optional item)

## Method
The TCGA-STAD NK-state DE test (Sec 3.5, Table 4) compares the 37
tumor-intrinsic candidate genes between NK-hot-cytotoxic (n=134) and
NK-hot-dysfunctional (n=20, underpowered per Limitations) tumors. Rather
than acquiring a new external dataset, we re-ran the identical labeling rule
(`src/immune_scoring/nk_scores.py`) on the two bulk gastric cohorts already
used for the effector-arm external validation (Sec 3.3): GSE62254 (ACRG,
n=300, GPL570) and GSE84437 (n=483, GPL6947).

## Cohort NK-state distributions
""")
        for name, (n_cyto, n_dysf) in cohort_ns.items():
            f.write(f"- **{name}**: NK-hot-cytotoxic n={n_cyto}, NK-hot-dysfunctional n={n_dysf}\n")
        f.write(f"""
## Directional concordance with the TCGA-STAD internal test
For the genes flagged in Table 4 (Tier 1 UP genes SGMS2/NT5E; Tier 3/X DOWN
genes PSPH/SHMT1/SHMT2/MTHFD1L/MTHFD1/RAC1), we compare the direction of the
dysfunctional-vs-cytotoxic difference in each external cohort against the
TCGA-STAD direction. Concordance rate: **{n_concordant}/{n_total}**
({100*n_concordant/n_total:.0f}% if n_total>0 else 'n/a').

See `nk_state_de_external_concordance.tsv` for the full per-gene, per-cohort
comparison, and `nk_state_de_external_replication.tsv` for DE results on all
37 genes in both cohorts.

## Caveats
- These external cohorts are microarray (not RNA-seq), and NK-state labels
  are derived independently per cohort from cohort-specific score
  distributions (median-based thresholds), not transferred from TCGA-STAD --
  so this is a test of directional consistency, not an exact replication of
  the same patients or platform.
- Sample sizes in the dysfunctional group vary by cohort and may still be
  small; see the per-cohort n's above.
""")
    log(f"Saved: {summary_path}")
    log("\nDONE.")


if __name__ == "__main__":
    main()
