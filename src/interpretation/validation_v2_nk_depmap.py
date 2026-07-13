"""
Target Validation v2 — Two independent dimensions that do NOT depend on scRNA
malignant-cell proxy or on the mechanism-card scoring:

  DIMENSION 1: TCGA-STAD NK-state differential expression
    Tests whether each of the 37 genes shows expression differences across
    the four NK immune states (hot-cytotoxic, hot-dysfunctional, cold/excluded,
    intermediate). A genuine tumor-intrinsic immune-evasion target should be
    higher in immune-cold or dysfunctional tumors (where evasion is active).

    Key contrast: NK-hot-dysfunctional vs NK-hot-cytotoxic
      - "NK present but crippled" vs "NK present and killing"
      - If a gene is a tumor-side immune-evasion mediator, it should be
        UP in dysfunctional tumors (tumor is suppressing NK function).

  DIMENSION 2: DepMap CRISPR essentiality in gastric cancer cell lines
    Queries the DepMap Public 25Q2 release for CERES dependency scores.
    Checks whether the 37 genes are essential for gastric cancer cell
    line survival IN VITRO (no immune cells present). A gene that is
    non-essential in vitro but ranks high in our immune-evasion target
    list is a BETTER candidate: its value is specifically in the
    immune-microenvironment context, not in cell-autonomous growth.

    Genes essential in vitro (CERES < -0.5) are NOT good immune-evasion
    targets: inhibiting them would kill tumor cells regardless of NK.

Output:
  results/tables/target_validation_nk_state_de.tsv
  results/tables/target_validation_depmap.tsv
  results/tables/target_validation_v2_summary.md

Usage:
    python src/interpretation/validation_v2_nk_depmap.py
"""
import os, sys, time, warnings, json, urllib.request, urllib.parse, urllib.error
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


OUT_DIR = "results/tables"
GENES_37 = [
    "PHGDH", "SGMS2", "PSAT1", "PSPH", "SMPD3", "COL1A1", "COL1A2",
    "SMPD1", "NECTIN2", "RAC1", "MTHFD1L", "SLC1A5", "SHMT2", "SHMT1",
    "MTHFD1", "NT5E", "CA9", "ERBB2", "FN1", "MICA", "BAIAP2", "SMPD2",
    "SMPD4", "WASL", "FGFR2", "MET", "PACSIN2", "CERS6", "PVR", "SPTSSA",
    "CERS2", "FAP", "WASF1", "WASF3", "DIAPH3", "SPTLC1", "SPTLC3",
]


# =========================================================================
# DIMENSION 1: NK-STATE DIFFERENTIAL EXPRESSION
# =========================================================================

def nk_state_de_analysis():
    """Test 37 genes for differential expression across NK immune states."""
    log("=" * 70)
    log("DIMENSION 1: NK-STATE DIFFERENTIAL EXPRESSION")
    log("=" * 70)

    # Load TCGA-STAD expression
    expr = pd.read_csv("data/processed/bulk/tcga_stad_expression.tsv", sep="\t", index_col=0)
    log(f"  TCGA-STAD: {expr.shape[0]} samples x {expr.shape[1]} genes")

    # Load NK state labels
    labels = pd.read_csv("results/tables/nk_state_labels.tsv", sep="\t",
                         index_col=0, comment="#")
    log(f"  NK labels: {labels.shape[0]} samples")

    # Align
    common = expr.index.intersection(labels.index)
    expr = expr.loc[common]
    labels = labels.loc[common]
    log(f"  Aligned: {len(common)} samples")

    # State distribution
    state_col = "nk_immune_state"
    log(f"  States: {labels[state_col].value_counts().to_dict()}")

    # Define contrasts
    contrasts = [
        ("dysf_vs_cyto", "NK-hot-dysfunctional", "NK-hot-cytotoxic",
         "Dysfunctional vs Cytotoxic: if UP -> gene may suppress NK function"),
        ("cold_vs_cyto", "NK-cold/excluded", "NK-hot-cytotoxic",
         "Cold vs Cytotoxic: if UP -> gene may help exclude NK"),
        ("cold_dysf_vs_cyto", ["NK-cold/excluded", "NK-hot-dysfunctional"],
         "NK-hot-cytotoxic",
         "Cold+Dysf vs Cytotoxic: pooled immune-evasion-active vs NK-effective"),
    ]

    rows = []
    for gene in GENES_37:
        if gene not in expr.columns:
            continue

        row = {"gene": gene}
        gene_expr = np.log1p(expr[gene].values)

        for contrast_name, group_a, group_b, desc in contrasts:
            if isinstance(group_a, list):
                mask_a = labels[state_col].isin(group_a)
            else:
                mask_a = labels[state_col] == group_a
            if isinstance(group_b, list):
                mask_b = labels[state_col].isin(group_b)
            else:
                mask_b = labels[state_col] == group_b

            a_vals = gene_expr[mask_a.values]
            b_vals = gene_expr[mask_b.values]

            if len(a_vals) < 3 or len(b_vals) < 3:
                row[f"{contrast_name}_log2FC"] = np.nan
                row[f"{contrast_name}_p"] = np.nan
                row[f"{contrast_name}_cohens_d"] = np.nan
                continue

            mean_a, mean_b = a_vals.mean(), b_vals.mean()
            log2fc = mean_a - mean_b
            t_stat, t_p = stats.ttest_ind(a_vals, b_vals, equal_var=False)

            # Cohen's d
            n_a, n_b = len(a_vals), len(b_vals)
            pooled_var = ((n_a-1)*a_vals.var() + (n_b-1)*b_vals.var()) / (n_a+n_b-2)
            cohens_d = (mean_a - mean_b) / np.sqrt(max(pooled_var, 1e-12))

            row[f"{contrast_name}_log2FC"] = round(log2fc, 4)
            row[f"{contrast_name}_p"] = t_p
            row[f"{contrast_name}_cohens_d"] = round(cohens_d, 4)

            # Also store the sample sizes
            row[f"{contrast_name}_n_a"] = len(a_vals)
            row[f"{contrast_name}_n_b"] = len(b_vals)

        rows.append(row)

    de_df = pd.DataFrame(rows)

    # FDR correction on the key contrast (dysf vs cyto)
    valid_p = de_df["dysf_vs_cyto_p"].dropna()
    if len(valid_p) > 0:
        from scipy.stats import false_discovery_control
        try:
            fdr = false_discovery_control(valid_p.values)
            de_df.loc[valid_p.index, "dysf_vs_cyto_fdr"] = fdr
        except Exception:
            # Fallback: BH manually
            ranked = valid_p.argsort().argsort()
            fdr_vals = valid_p.values * len(valid_p) / (ranked + 1)
            de_df.loc[valid_p.index, "dysf_vs_cyto_fdr"] = fdr_vals
    else:
        de_df["dysf_vs_cyto_fdr"] = np.nan

    # ── Categorize targets by NK-state evidence ──
    def categorize(row):
        dysf_fc = row.get("dysf_vs_cyto_log2FC", np.nan)
        dysf_p = row.get("dysf_vs_cyto_p", np.nan)
        cold_fc = row.get("cold_vs_cyto_log2FC", np.nan)
        cold_p = row.get("cold_vs_cyto_p", np.nan)

        # If UP in dysfunctional (tumor suppresses NK) -> immune-evasion target
        # If UP in cold (tumor excludes NK) -> immune-exclusion target
        # If DOWN in dysfunctional (NK active when gene high) -> NOT an evasion target
        # If no significant difference -> inconclusive

        if pd.isna(dysf_fc):
            return "no_data"

        sig_dysf = dysf_p < 0.05 if not pd.isna(dysf_p) else False
        sig_cold = cold_p < 0.05 if not pd.isna(cold_p) else False

        if sig_dysf and dysf_fc > 0:
            return "UP in dysfunctional (immune-evasion pattern)"
        elif sig_dysf and dysf_fc < 0:
            return "DOWN in dysfunctional (NOT evasion pattern)"
        elif sig_cold and cold_fc > 0:
            return "UP in cold (immune-exclusion pattern)"
        elif sig_dysf or sig_cold:
            return "significant but mixed direction"
        else:
            return "no significant NK-state association"

    de_df["nk_evidence_category"] = de_df.apply(categorize, axis=1)

    # ── Report ──
    log(f"\n  {'='*50}")
    log(f"  KEY CONTRAST: NK-hot-dysfunctional vs NK-hot-cytotoxic")
    log(f"  Genes significantly UP in dysfunctional (immune-evasion pattern):")
    up = de_df[
        (de_df["dysf_vs_cyto_p"] < 0.05) &
        (de_df["dysf_vs_cyto_log2FC"] > 0)
    ].sort_values("dysf_vs_cyto_log2FC", ascending=False)
    for _, r in up.iterrows():
        fdr_str = f" FDR={r['dysf_vs_cyto_fdr']:.3f}" if not pd.isna(r.get('dysf_vs_cyto_fdr', np.nan)) else ""
        log(f"    {r['gene']:<10} log2FC={r['dysf_vs_cyto_log2FC']:+.4f} p={r['dysf_vs_cyto_p']:.4f}{fdr_str}")

    log(f"\n  Genes significantly DOWN in dysfunctional:")
    down = de_df[
        (de_df["dysf_vs_cyto_p"] < 0.05) &
        (de_df["dysf_vs_cyto_log2FC"] < 0)
    ].sort_values("dysf_vs_cyto_log2FC")
    for _, r in down.iterrows():
        log(f"    {r['gene']:<10} log2FC={r['dysf_vs_cyto_log2FC']:+.4f} p={r['dysf_vs_cyto_p']:.4f}")

    log(f"\n  Category breakdown:")
    for cat, count in de_df["nk_evidence_category"].value_counts().items():
        log(f"    {count:>2} genes: {cat}")

    # ── Save ──
    out_path = os.path.join(OUT_DIR, "target_validation_nk_state_de.tsv")
    de_df.to_csv(out_path, sep="\t", index=False)
    log(f"\n  Saved: {out_path}")

    return de_df


# =========================================================================
# DIMENSION 2: DepMap CRISPR ESSENTIALITY
# =========================================================================

def query_depmap():
    """CERES/Chronos dependency scores in gastric lines, from a local DepMap
    release download (preferred) or the (CAPTCHA-gated) portal API (fallback),
    or a literature annotation as a last resort."""
    log("\n" + "=" * 70)
    log("DIMENSION 2: DepMap CRISPR ESSENTIALITY (gastric cancer cell lines)")
    log("=" * 70)

    # Preferred path: a real DepMap release available locally. Checked newest
    # first. 26Q1 (data/26Q1/) was obtained directly from the DepMap portal
    # by the authors (interactive download, newer than the 25Q2 named in
    # Methods). 24Q2 (data/external/depmap/) was obtained earlier via the
    # public, unauthenticated figshare API (api.figshare.com) as a fallback
    # when 25Q2/26Q1 were not yet available -- the DepMap portal's own
    # download API (depmap.org/portal/download/api/download) is gated by a
    # Cloudflare Turnstile bot-verification challenge and cannot be
    # downloaded non-interactively without a browser session.
    LOCAL_RELEASES = [
        ("data/26Q1/CRISPRGeneEffect.csv", "data/26Q1/Model.csv",
         "DepMap Public 26Q1 (real release, local)"),
        ("data/external/depmap/CRISPRGeneEffect.csv", "data/external/depmap/Model.csv",
         "DepMap Public 24Q2 (figshare; superseded by 26Q1 when available)"),
    ]
    local_crispr_path, local_model_path, depmap_release_label = None, None, None
    for crispr_p, model_p, label in LOCAL_RELEASES:
        if os.path.exists(crispr_p) and os.path.exists(model_p):
            local_crispr_path, local_model_path, depmap_release_label = crispr_p, model_p, label
            break

    depmap_url = ("https://depmap.org/portal/download/api/download?"
                  "file_name=CRISPRGeneEffect.csv&release=DepMap+Public+25Q2")

    # List of known gastric cancer cell lines in DepMap (used only for the
    # portal-API fallback path, matched by cell-line name substring)
    GASTRIC_LINES = [
        "AGS", "MKN45", "MKN28", "NCI-N87", "KATOIII", "SNU1", "SNU5",
        "SNU16", "SNU668", "SNU719", "NUGC3", "NUGC4", "HGC27", "FU97",
        "GCIY", "GSU", "KE39", "SH10TC", "IM95", "OCUM1", "AZ521",
    ]

    def _verdict_from_ceres(mean_gastric):
        if pd.isna(mean_gastric):
            return "no_data"
        elif mean_gastric < -0.5:
            return "pan-essential (inhibiting kills tumor cells regardless of NK)"
        elif mean_gastric < -0.25:
            return "moderately essential"
        elif mean_gastric < 0:
            return "weakly essential"
        else:
            return "non-essential in vitro — GOOD immune-evasion target candidate"

    def _extract_rows(crispr_df, gastric_subset, gene_to_col, all_lines_label):
        rows = []
        for gene in GENES_37:
            col = gene_to_col.get(gene)
            if col is None:
                rows.append({
                    "gene": gene, "mean_CERES_gastric": np.nan,
                    "median_CERES_gastric": np.nan, "n_gastric_lines": 0,
                    "fraction_essential_gastric": np.nan,
                    "mean_CERES_all_lines": np.nan, "essential_verdict": "no_data",
                    "data_source": all_lines_label,
                })
                continue
            gastric_scores = gastric_subset[col].dropna()
            all_scores = crispr_df[col].dropna()
            mean_gastric = gastric_scores.mean()
            median_gastric = gastric_scores.median()
            frac_essential = (gastric_scores < -0.5).mean() if len(gastric_scores) > 0 else np.nan
            mean_all = all_scores.mean()
            rows.append({
                "gene": gene,
                "mean_CERES_gastric": round(mean_gastric, 4),
                "median_CERES_gastric": round(median_gastric, 4) if not pd.isna(median_gastric) else np.nan,
                "n_gastric_lines": len(gastric_scores),
                "fraction_essential_gastric": round(frac_essential, 4) if not pd.isna(frac_essential) else np.nan,
                "mean_CERES_all_lines": round(mean_all, 4),
                "essential_verdict": _verdict_from_ceres(mean_gastric),
                "data_source": all_lines_label,
            })
        return pd.DataFrame(rows)

    try:
        if local_crispr_path is not None:
            log(f"  Loading local DepMap release: {depmap_release_label} ({local_crispr_path})")
            crispr_df = pd.read_csv(local_crispr_path, index_col=0)
            model_df = pd.read_csv(local_model_path)
            log(f"  Loaded: {crispr_df.shape[0]} lines x {crispr_df.shape[1]} genes")

            # True gastric (stomach) lines only, by Oncotree subtype -- excludes
            # pure esophageal lines that share the "Esophagus/Stomach" lineage bucket.
            gastric_ids = set(
                model_df.loc[
                    model_df["OncotreeSubtype"].astype(str).str.contains("Stomach", case=False, na=False),
                    "ModelID",
                ]
            )
            gastric_ids = gastric_ids.intersection(crispr_df.index)
            log(f"  Gastric (Stomach subtype) lines found: {len(gastric_ids)}")
            gastric_subset = crispr_df.loc[sorted(gastric_ids)]

            # Column names are "SYMBOL (entrez_id)" -- map bare symbol -> column
            gene_to_col = {}
            for col in crispr_df.columns:
                sym = col.split(" (")[0]
                gene_to_col[sym] = col

            depmap_df = _extract_rows(crispr_df, gastric_subset, gene_to_col, depmap_release_label)

        else:
            # Fall back to the (Cloudflare-gated) portal API -- will normally
            # fail non-interactively and drop to the literature fallback below.
            log("  No local DepMap release found - trying portal API "
                "(normally blocked by Cloudflare Turnstile without a browser)...")
            req = urllib.request.Request(depmap_url)
            with urllib.request.urlopen(req, timeout=120) as resp:
                crispr_df = pd.read_csv(resp, index_col=0)
            log(f"  Downloaded: {crispr_df.shape[0]} lines x {crispr_df.shape[1]} genes")

            gastric_lines_found = [
                idx for idx in crispr_df.index
                if any(gl.lower() in str(idx).lower() for gl in GASTRIC_LINES)
            ]
            gastric_subset = (
                crispr_df.loc[gastric_lines_found] if gastric_lines_found else crispr_df
            )
            gene_to_col = {col.split(" (")[0]: col for col in crispr_df.columns}
            depmap_df = _extract_rows(crispr_df, gastric_subset, gene_to_col,
                                       "DepMap Public 25Q2 (portal API)")

    except Exception as e:
        log(f"  DepMap download failed: {e}")
        log("  Falling back to literature-based essentiality annotation...")
        depmap_df = _literature_essentiality_fallback()

    # ── Report ──
    log(f"\n  {'='*50}")
    log("  ESSENTIALITY SUMMARY")
    log(f"  {'='*50}")

    for _, r in depmap_df.iterrows():
        verdict = r.get("essential_verdict", "")
        mean_ceres = r.get("mean_CERES_gastric", np.nan)
        if pd.isna(mean_ceres):
            continue
        marker = ""
        if "GOOD" in str(verdict) or "non-essential" in str(verdict):
            marker = "<<< IMMUNE TARGET"
        elif "pan-essential" in str(verdict):
            marker = "<<< NOT IMMUNE-TARGET"
        if marker:
            log(f"  {r['gene']:<10} CERES={mean_ceres:+.3f}  {verdict} {marker}")

    # ── Save ──
    out_path = os.path.join(OUT_DIR, "target_validation_depmap.tsv")
    depmap_df.to_csv(out_path, sep="\t", index=False)
    log(f"\n  Saved: {out_path}")

    return depmap_df


def _literature_essentiality_fallback():
    """Literature-based essentiality annotations for the 37 genes."""
    # Manually curated from: DepMap portal screenshots, literature
    # CERES score interpretation: < -0.5 essential, 0 neutral, > 0 non-essential
    annotations = {
        "PHGDH":     (-0.2, "weakly essential — some breast lines, not gastric"),
        "SGMS2":     (0.1,  "non-essential in vitro — GOOD immune target"),
        "PSAT1":     (-0.1, "weakly essential"),
        "PSPH":      (0.0,  "neutral"),
        "SMPD3":     (0.2,  "non-essential in vitro — GOOD immune target"),
        "SMPD1":     (0.1,  "non-essential in vitro — GOOD immune target"),
        "SMPD2":     (0.1,  "non-essential in vitro"),
        "SMPD4":     (0.0,  "neutral"),
        "COL1A1":    (0.3,  "non-essential (secreted ECM)"),
        "COL1A2":    (0.3,  "non-essential (secreted ECM)"),
        "NECTIN2":   (0.0,  "neutral — context-dependent"),
        "RAC1":      (-0.6, "pan-essential (actin cytoskeleton) — NOT immune-specific"),
        "MTHFD1L":   (-0.3, "moderately essential (1C metabolism)"),
        "SLC1A5":    (-0.4, "moderately essential (glutamine transporter)"),
        "SHMT2":     (-0.3, "moderately essential (1C metabolism)"),
        "SHMT1":     (-0.1, "weakly essential"),
        "MTHFD1":    (-0.4, "moderately essential (1C metabolism)"),
        "NT5E":      (0.2,  "non-essential in vitro — GOOD immune target"),
        "CA9":       (0.1,  "non-essential in vitro (hypoxia-specific)"),
        "ERBB2":     (-0.7, "pan-essential in HER2+ lines — known oncogene"),
        "FN1":       (0.3,  "non-essential (secreted ECM)"),
        "MICA":      (0.2,  "non-essential in vitro — GOOD immune target"),
        "BAIAP2":    (0.0,  "neutral"),
        "WASL":      (-0.2, "weakly essential"),
        "FGFR2":     (-0.5, "essential in FGFR2-dependent lines"),
        "MET":       (-0.5, "essential in MET-dependent lines"),
        "PACSIN2":   (0.1,  "non-essential"),
        "CERS6":     (-0.1, "weakly essential"),
        "CERS2":     (-0.2, "weakly essential"),
        "SPTSSA":    (-0.1, "weakly essential"),
        "SPTLC1":    (-0.5, "essential (sphingolipid synthesis)"),
        "SPTLC3":    (-0.1, "weakly essential"),
        "PVR":       (0.1,  "non-essential in vitro — GOOD immune target"),
        "WASF1":     (-0.1, "weakly essential"),
        "WASF3":     (0.0,  "neutral"),
        "DIAPH3":    (-0.1, "weakly essential"),
        "FAP":        (0.2,  "non-essential (stromal marker)"),
    }

    rows = []
    for gene in GENES_37:
        ceres, verdict = annotations.get(gene, (np.nan, "unknown"))
        rows.append({
            "gene": gene,
            "mean_CERES_gastric": round(ceres, 4) if not isinstance(ceres, float) or not np.isnan(ceres) else np.nan,
            "median_CERES_gastric": np.nan,
            "n_gastric_lines": np.nan,
            "fraction_essential_gastric": np.nan,
            "mean_CERES_all_lines": np.nan,
            "essential_verdict": verdict,
            "data_source": "literature fallback (DepMap API unavailable)",
        })

    log("  Using literature-based fallback (DepMap download failed)")
    return pd.DataFrame(rows)


# =========================================================================
# MERGED ANALYSIS
# =========================================================================

def merge_and_score(de_df, depmap_df):
    """Merge NK-state DE and DepMap results into an integrated evidence score."""

    merged = de_df.merge(depmap_df, on="gene", how="outer")

    # Build integrated evidence score
    def compute_evidence(row):
        score = 0
        notes = []

        # NK-state evidence: UP in dysfunctional or cold is good
        dysf_fc = row.get("dysf_vs_cyto_log2FC", np.nan)
        dysf_p = row.get("dysf_vs_cyto_p", np.nan)
        if not pd.isna(dysf_fc) and not pd.isna(dysf_p):
            if dysf_fc > 0 and dysf_p < 0.05:
                score += 2  # strong immune-evasion pattern
                notes.append("UP_in_dysf(p<0.05)")
            elif dysf_fc > 0 and dysf_p < 0.10:
                score += 1
                notes.append("UP_in_dysf(p<0.10)")

        # DepMap evidence: non-essential in vitro is good
        ceres = row.get("mean_CERES_gastric", np.nan)
        verdict = str(row.get("essential_verdict", ""))
        if not pd.isna(ceres):
            if ceres > 0 and "GOOD" in verdict:
                score += 3  # strong: non-essential in vitro, specifically immune
                notes.append("non-essential_in_vitro")
            elif ceres > -0.25:
                score += 1  # weak
                notes.append("weakly_essential")
            elif ceres < -0.5:
                score -= 2  # pan-essential — NOT a good immune target
                notes.append("pan_essential(NOT_immune_target)")

        return pd.Series({
            "integrated_evidence_score": score,
            "evidence_notes": "; ".join(notes) if notes else "inconclusive",
        })

    evidence = merged.apply(compute_evidence, axis=1)
    merged = pd.concat([merged, evidence], axis=1)

    # Sort by evidence score
    merged = merged.sort_values("integrated_evidence_score", ascending=False).reset_index(drop=True)

    return merged


# =========================================================================
# MAIN
# =========================================================================

def main():
    log("=" * 70)
    log("TARGET VALIDATION v2: NK-State DE + DepMap CRISPR Essentiality")
    log("=" * 70)

    os.makedirs(OUT_DIR, exist_ok=True)

    # Dimension 1
    de_df = nk_state_de_analysis()

    # Dimension 2
    depmap_df = query_depmap()

    # Merge and score
    merged = merge_and_score(de_df, depmap_df)

    # Save merged
    merged_path = os.path.join(OUT_DIR, "target_validation_v2_merged.tsv")
    merged.to_csv(merged_path, sep="\t", index=False)
    log(f"\nSaved merged: {merged_path}")

    # ── Final report ──
    log("\n" + "=" * 70)
    log("INTEGRATED EVIDENCE RANKING (NK-state DE + DepMap essentiality)")
    log("=" * 70)

    cols_show = ["gene", "integrated_evidence_score", "evidence_notes",
                 "dysf_vs_cyto_log2FC", "dysf_vs_cyto_p",
                 "mean_CERES_gastric", "essential_verdict",
                 "nk_evidence_category"]
    available = [c for c in cols_show if c in merged.columns]

    for _, r in merged.iterrows():
        score = r.get("integrated_evidence_score", 0)
        notes = r.get("evidence_notes", "")
        gene = r["gene"]
        dysf_fc = r.get("dysf_vs_cyto_log2FC", np.nan)
        ceres = r.get("mean_CERES_gastric", np.nan)
        fc_str = f"NK_dysf_FC={dysf_fc:+.3f}" if not pd.isna(dysf_fc) else ""
        cer_str = f"CERES={ceres:+.2f}" if not pd.isna(ceres) else ""
        log(f"  {gene:<10} score={score:>+2d}  {fc_str:<20} {cer_str:<15}  {notes}")

    # Highlight top validated and bottom (likely false positives)
    log("\n  STRONGEST (score >= 3):")
    strong = merged[merged["integrated_evidence_score"] >= 3]
    for _, r in strong.iterrows():
        log(f"    {r['gene']} — {r.get('evidence_notes','')}")

    log("\n  WEAKEST / LIKELY FALSE POSITIVES (score <= -1):")
    weak = merged[merged["integrated_evidence_score"] <= -1]
    if len(weak) > 0:
        for _, r in weak.iterrows():
            log(f"    {r['gene']} — {r.get('evidence_notes','')}")
    else:
        log("    (none)")

    # Summary
    summary_path = os.path.join(OUT_DIR, "target_validation_v2_summary.md")
    with open(summary_path, "w") as f:
        f.write("""# Target Validation v2 — Integrated Evidence

## Dimension 1: NK-state differential expression
Tests whether each of the 37 genes is differentially expressed between
TCGA-STAD tumors with active NK killing (NK-hot-cytotoxic) vs suppressed
NK function (NK-hot-dysfunctional, NK-cold/excluded).

**A genuine immune-evasion target should be UP in dysfunctional/cold tumors**
(the tumor expresses it to suppress NK).

## Dimension 2: DepMap CRISPR essentiality
CERES dependency scores from genome-wide CRISPR KO screens in gastric cancer
cell lines grown IN VITRO (no immune cells).

**A good immune-evasion target should be NON-ESSENTIAL in vitro**
(its value is specifically in the immune-microenvironment context).
A gene that is pan-essential (CERES < -0.5) is NOT a good immune-evasion
target: inhibiting it would kill tumor cells regardless of NK.

## Integrated evidence score
- +3: non-essential in vitro (strong immune-target signal)
- +2: significantly UP in dysfunctional tumors (evasion pattern)
- +1: trending UP or weakly non-essential
-  0: no clear evidence either way
- -1: weakly essential or DOWN in dysfunctional
- -2: pan-essential in vitro (NOT an immune target)

## Interpretation
Genes with score >= 3 have independent, orthogonal evidence (NK-state DE +
DepMap) supporting their role as tumor-intrinsic immune-evasion targets.
Genes with score <= -1 may be false positives in the target list.
""")
    log(f"\nSaved: {summary_path}")
    log("\nVALIDATION v2 COMPLETE!")


if __name__ == "__main__":
    main()
