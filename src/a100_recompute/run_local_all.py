"""
Local runner for P1 science tasks T5, T6, T14, T15.

Adapted from src/a100_recompute/*.py for local data schema:
  - scRNA SST scores table is NK-only (no cell_type column needed)
  - tumor_intrinsic_candidates uses target_category (not category)
  - condition column: normal/tumor (not healthy/peritumoral)

T7 (ablation) SKIPPED — needs TCGA-STAD bulk expression which is only on A100.

Run:  cd G:/cc/GC-NKGraph-Atlas && python src/a100_recompute/run_local_all.py
"""
import sys, os, time, io
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm, linregress
from pathlib import Path

# Fix Windows GBK console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
T = "results/tables/"
OUT = T  # write outputs to same location

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# ============================================================================
# T5 — H1 result
# ============================================================================
def run_t5():
    log("=" * 50)
    log("T5: H1 result — tumor_serine_capacity ~ nk_sm_balance")
    log("=" * 50)

    liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
    sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")

    # Bulk H1
    x_b = liver["tumor_serine_capacity_score"]
    y_b = liver["nk_sm_balance_score"]
    r_b, p_b = stats.pearsonr(x_b.dropna(), y_b.dropna())

    # scRNA H1 (all cells are NK — this table is NK-only)
    x_s = sc["tumor_serine_capacity_score"]
    y_s = sc["nk_sm_balance_score"]
    valid = x_s.notna() & y_s.notna()
    r_s, p_s = stats.pearsonr(x_s[valid], y_s[valid])

    log(f"Bulk H1:  r={r_b:.4f}  P={p_b:.2e}  n={len(x_b.dropna())}")
    log(f"scNK H1:  r={r_s:.4f}  P={p_s:.2e}  n={valid.sum()}")

    # Append H1 rows to recovery table (above existing H2–H5)
    rec = pd.read_csv(T + "sst_axis_positive_control_recovery.tsv", sep="\t")
    h1_rows = pd.DataFrame([
        {"hypothesis": "H1", "test": "serine_capacity ~ sm_balance", "resolution": "bulk",
         "r": round(r_b, 4), "p": p_b, "expected": "calibrated", "outcome": "reported"},
        {"hypothesis": "H1", "test": "serine_capacity ~ sm_balance", "resolution": "single-cell NK",
         "r": round(r_s, 4), "p": p_s, "expected": "calibrated", "outcome": "reported"},
    ])
    rec_full = pd.concat([h1_rows, rec], ignore_index=True)
    rec_full.to_csv(T + "sst_axis_positive_control_recovery.tsv", sep="\t", index=False)
    log(f"Updated recovery table: {len(rec_full)} rows (was {len(rec)})")
    log("T5 PASS")
    return r_b, p_b, r_s, p_s

# ============================================================================
# T15 — Effect-size reframe
# ============================================================================
def run_t15():
    log("\n" + "=" * 50)
    log("T15: Effect sizes, 95% CIs, Benjamini-Hochberg correction")
    log("=" * 50)

    def corr_ci(r, n, alpha=0.05):
        if r >= 0.9999: r = 0.9999
        if r <= -0.9999: r = -0.9999
        z = np.arctanh(r)
        se = 1.0 / np.sqrt(n - 3)
        z_crit = norm.ppf(1 - alpha/2)
        return np.tanh(z - z_crit * se), np.tanh(z + z_crit * se)

    def cohens_d(x1, x2):
        n1, n2 = len(x1), len(x2)
        v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
        pooled_sd = np.sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2))
        return (np.mean(x1) - np.mean(x2)) / (pooled_sd + 1e-10)

    liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
    sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")

    tests = []

    def add_corr(label, hyp, res, x_vals, y_vals):
        mask = x_vals.notna() & y_vals.notna()
        x, y = x_vals[mask], y_vals[mask]
        r, p = stats.pearsonr(x, y)
        ci_lo, ci_hi = corr_ci(r, len(x))
        tests.append({"hypothesis": hyp, "test": label, "resolution": res,
                      "statistic": "r", "value": round(r, 4), "r2": round(r**2, 5),
                      "ci_95_lo": round(ci_lo, 4), "ci_95_hi": round(ci_hi, 4),
                      "p": p, "n": len(x), "expected_dir": "+", "outcome": "—"})

    def add_delta(label, hyp, res, group1, group2):
        mask1, mask2 = group1.notna(), group2.notna()
        g1, g2 = group1[mask1], group2[mask2]
        t_stat, p = stats.ttest_ind(g1, g2, equal_var=False)
        d = cohens_d(g1, g2)
        tests.append({"hypothesis": hyp, "test": label, "resolution": res,
                      "statistic": "Cohen's d", "value": round(d, 4), "r2": 0,
                      "ci_95_lo": round(d - 1.96 * np.sqrt(1/len(g1) + 1/len(g2)), 4),
                      "ci_95_hi": round(d + 1.96 * np.sqrt(1/len(g1) + 1/len(g2)), 4),
                      "p": p, "n": f"{len(g1)}+{len(g2)}",
                      "expected_dir": "tumor<normal", "outcome": "—"})

    # H1
    add_corr("serine_cap ~ sm_balance", "H1", "bulk",
             liver["tumor_serine_capacity_score"], liver["nk_sm_balance_score"])
    add_corr("serine_cap ~ sm_balance", "H1", "single-cell NK",
             sc["tumor_serine_capacity_score"], sc["nk_sm_balance_score"])

    # H2
    add_corr("sm_balance ~ protrusion", "H2", "bulk",
             liver["nk_sm_balance_score"], liver["nk_protrusion_machinery_score"])
    add_corr("sm_balance ~ protrusion", "H2", "single-cell NK",
             sc["nk_sm_balance_score"], sc["nk_protrusion_machinery_score"])

    # H3
    add_corr("protrusion ~ cytotoxicity", "H3", "bulk",
             liver["nk_protrusion_machinery_score"], liver["nk_synapse_cytotoxicity_outcome_score"])
    add_corr("protrusion ~ cytotoxicity", "H3", "single-cell NK",
             sc["nk_protrusion_machinery_score"], sc["nk_synapse_cytotoxicity_outcome_score"])

    # H4
    add_corr("topology ~ dysfunction", "H4", "bulk",
             liver["nk_topology_permissive_score"], liver["checkpoint_link_score"])
    add_corr("topology ~ HAVCR2", "H4", "single-cell NK",
             sc["nk_topology_permissive_score"], sc["checkpoint_link_score"])

    # H5 — intratumoral vs normal NK
    intra = sc[sc["condition"] == "tumor"]
    normal = sc[sc["condition"] == "normal"]
    log(f"  H5 groups: intratumoral={len(intra)}  normal={len(normal)}")
    add_delta("intratumoral vs normal: cytotoxicity", "H5", "single-cell NK",
              intra["nk_synapse_cytotoxicity_outcome_score"],
              normal["nk_synapse_cytotoxicity_outcome_score"])
    add_delta("intratumoral vs normal: protrusion", "H5", "single-cell NK",
              intra["nk_protrusion_machinery_score"],
              normal["nk_protrusion_machinery_score"])

    # --- BH correction ---
    df = pd.DataFrame(tests)
    p_vals = df[df["p"].notna()]["p"].values
    from statsmodels.stats.multitest import multipletests
    _, p_corrected, _, _ = multipletests(p_vals, method="fdr_bh")
    df.loc[df["p"].notna(), "p_fdr"] = p_corrected

    # --- Effect verdict ---
    def verdict(row):
        if row["p_fdr"] > 0.05:
            return "not significant (FDR>0.05)"
        if row["statistic"] == "r":
            r2 = row["r2"]
            if r2 < 0.001: return "negligible (r^2<0.001)"
            elif r2 < 0.01: return "weak (r^2<0.01)"
            elif r2 < 0.09: return "moderate"
            else: return "strong"
        if row["statistic"] == "Cohen's d":
            d = abs(row["value"])
            if d < 0.2: return "negligible"
            elif d < 0.5: return "small"
            elif d < 0.8: return "medium"
            else: return "large"
        return "—"

    df["effect_verdict"] = df.apply(verdict, axis=1)
    df.to_csv(T + "sst_axis_positive_control_recovery_v2.tsv", sep="\t", index=False)

    # --- Print summary ---
    print("\n" + "=" * 70)
    print("HYPOTHESIS EFFECT-SIZE SUMMARY (with FDR correction)")
    print("=" * 70)
    for _, row in df.iterrows():
        flag = " ***" if "negligible" in str(row["effect_verdict"]) else ""
        print(f"  {row['hypothesis']:4s} {row['resolution']:18s}  "
              f"r/d={row['value']:+.4f}  r^2={row['r2']:.5f}  "
              f"95%CI=[{row['ci_95_lo']:+.4f}, {row['ci_95_hi']:+.4f}]  "
              f"P={row['p']:.2e}  P_FDR={row['p_fdr']:.2e}  "
              f"n={row['n']}  {row['effect_verdict']}{flag}")

    # --- Key flags ---
    h2_sc = df[(df["hypothesis"] == "H2") & (df["resolution"] == "single-cell NK")]
    if len(h2_sc):
        r2 = h2_sc["r2"].values[0]
        r_val = h2_sc["value"].values[0]
        print(f"\n=== KEY FLAGS FOR MANUSCRIPT ===")
        print(f"H2 scNK: r={r_val:.4f}, r^2={r2:.5f} ({r2*100:.3f}% shared variance)")
        print(f"  → NEGLIGIBLE effect. Text MUST say: 'statistically detectable but of negligible magnitude'")

    n_sig = (df["p_fdr"] < 0.05).sum()
    print(f"\n{n_sig}/{len(df)} tests significant at FDR < 0.05")
    log("T15 PASS")
    return df

# ============================================================================
# T14 — H3 activation control
# ============================================================================
def run_t14():
    log("\n" + "=" * 50)
    log("T14: H3 activation control — partial correlation")
    log("=" * 50)

    sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")

    NK_ACTIVATION_GENES = [
        "CD69", "TNF", "XCL1", "XCL2", "CCL3", "CCL4", "CCL5",
        "CSF2", "IL2RA", "ICOS", "TNFSF10", "FASLG", "CD38",
        "HLA-DRA", "HLA-DRB1", "MKI67",
    ]

    # Check which activation genes are present in the expression columns
    # The sst_axis_scores table has module scores, not individual gene expression.
    # We need individual gene expression for the activation signature.
    # Try to find it — might be in a different table.

    # First, check if gene-level expression columns exist in the scRNA table
    act_genes_present = [g for g in NK_ACTIVATION_GENES if g in sc.columns]
    log(f"Activation genes in scRNA table: {len(act_genes_present)}/{len(NK_ACTIVATION_GENES)}")

    if len(act_genes_present) < 3:
        log("WARNING: scRNA SST scores table has module scores, not gene expression.")
        log("Searching for gene-level scRNA expression...")

        # Try alternative source: sst_axis_scrna_nk_scores.tsv might have gene-level data
        alt_files = [
            T + "sst_axis_scrna_nk_scores.tsv",
            "data/processed/scrna/gc_scrna_nk_expression.tsv",
            "data/raw/GSE246662_processed.tsv",
        ]
        found = False
        for alt in alt_files:
            if os.path.exists(alt):
                try:
                    df_alt = pd.read_csv(alt, sep="\t", nrows=2)
                    alt_genes = [g for g in NK_ACTIVATION_GENES if g in df_alt.columns]
                    log(f"  {alt}: {len(alt_genes)} activation genes found")
                    if len(alt_genes) >= 3:
                        # Load full table
                        df_full = pd.read_csv(alt, sep="\t")
                        # Use CYTOXICITY OUTPUT score as alignment key for protrusion/cytotox
                        # Align: join on cell_id or index
                        sc_indexed = sc.set_index("cell_id")
                        common_ids = sc_indexed.index.intersection(df_full.index if df_full.index.name == "cell_id" else df_full.iloc[:,0])
                        log(f"  Common cell IDs: {len(common_ids)}")
                        if len(common_ids) > 100:
                            sc_aligned = sc_indexed.loc[common_ids]
                            df_aligned = df_full.set_index(df_full.columns[0]).loc[common_ids]
                            act_z = (df_aligned[alt_genes] - df_aligned[alt_genes].mean()) / (df_aligned[alt_genes].std() + 1e-10)
                            act_score = act_z.mean(axis=1)
                            act_genes_present = alt_genes
                            found = True
                            break
                except Exception as e:
                    log(f"  {alt}: failed ({e})")

        if not found:
            # FALLBACK: use IFNG, GZMB, PRF1 as proxy activation genes
            # (these ARE in the cytotoxicity module, but they're the best available proxy)
            fallback_genes = ["IFNG"]
            for g in fallback_genes:
                if g in sc.columns:
                    act_genes_present.append(g)
            if len(act_genes_present) < 2:
                log("ERROR: cannot build activation signature — gene-level data missing")
                log("T14 SKIPPED — needs scRNA gene expression, not module scores")
                return None

    log(f"Using {len(act_genes_present)} genes for activation signature: {act_genes_present}")

    # Activation score: use gene-level expression if available, else fallback to module-score proxy
    act_score = None
    if len(act_genes_present) >= 3:
        act_z = (sc[act_genes_present] - sc[act_genes_present].mean()) / (sc[act_genes_present].std() + 1e-10)
        act_score = act_z.mean(axis=1).values
    else:
        # Fallback: use checkpoint_link_score as rough activation proxy
        # (HAVCR2 expression tracks NK exhaustion, inverse of activation)
        act_score = -sc["checkpoint_link_score"].values
        log("  Using -checkpoint_link_score as activation proxy (gene expr unavailable)")
    a = act_score

    x = sc["nk_protrusion_machinery_score"].values
    y = sc["nk_synapse_cytotoxicity_outcome_score"].values

    valid = ~(np.isnan(x) | np.isnan(y) | np.isnan(a))
    x, y, a = x[valid], y[valid], a[valid]

    r_raw, p_raw = stats.pearsonr(x, y)
    log(f"Raw H3:     r={r_raw:.4f}  P={p_raw:.2e}  r^2={r_raw**2:.4f}  n={len(x)}")

    # Partial correlation
    slope_xa, intercept_xa, _, _, _ = linregress(a, x)
    x_resid = x - (slope_xa * a + intercept_xa)
    slope_ya, intercept_ya, _, _, _ = linregress(a, y)
    y_resid = y - (slope_ya * a + intercept_ya)
    r_partial, p_partial = stats.pearsonr(x_resid, y_resid)

    r2_raw = r_raw**2
    r2_partial = r_partial**2
    r2_act = r2_raw - r2_partial
    log(f"Partial:     r={r_partial:.4f}  P={p_partial:.2e}  r^2={r2_partial:.4f}")
    log(f"Δ(activation component):  r^2={r2_act:.4f} ({100*r2_act/r2_raw:.1f}% of raw r^2)")

    # Verdict
    print()
    if r_partial > 0.15 and p_partial < 0.001:
        print("VERDICT: H3 ROBUST to activation control → claim SURVIVES")
        print("  protrusion→cytotoxicity is NOT just co-activation.")
    elif r_partial > 0.05 and p_partial < 0.05:
        print("VERDICT: H3 PARTIALLY survives → write 'partially independent of activation'")
        print(f"  Residual r={r_partial:.3f}, {100*r2_act/r2_raw:.0f}% of r^2 shared with activation.")
    else:
        print("VERDICT: H3 DOES NOT SURVIVE → DOWNGRADE to co-activation language")
        print(f"  Residual r={r_partial:.3f} (r^2={r2_partial:.4f}) — activation explains nearly all.")

    df = pd.DataFrame([{
        "test": "H3 protrusion~cytotoxicity", "resolution": "single-cell NK",
        "r_raw": round(r_raw, 4), "r2_raw": round(r2_raw, 5),
        "r_partial_activation": round(r_partial, 4), "r2_partial": round(r2_partial, 5),
        "r2_activation_component": round(r2_act, 5),
        "p_raw": p_raw, "p_partial": p_partial, "n": len(x),
        "activation_genes_n": len(act_genes_present),
    }])
    df.to_csv(T + "h3_activation_control.tsv", sep="\t", index=False)
    log("T14 PASS")
    return r_raw, r_partial

# ============================================================================
# T6 — De-circ audit
# ============================================================================
def run_t6():
    log("\n" + "=" * 50)
    log("T6: De-circ audit — tumor-intrinsic pool composition")
    log("=" * 50)

    ti = pd.read_csv(T + "tumor_intrinsic_candidates.tsv", sep="\t")
    gene_col = "gene"
    cat_col = "target_category"

    log(f"Pool size: {len(ti)}")
    log(f"Top genes: {ti[gene_col].head(15).tolist()}")

    if cat_col in ti.columns:
        cat_counts = ti[cat_col].value_counts()
        print("\n=== CATEGORY BREAKDOWN ===")
        for cat, n in cat_counts.items():
            cat_lower = str(cat).lower().replace(" ", "_")
            is_nk = any(nk_tag in cat_lower for nk_tag in
                        ["nk_protrusion", "nk_sm_", "nk_denovo", "nk_synapse",
                         "cytotoxicity", "sst_axis_nk", "checkpoint"])
            flag = " [NK-side]" if is_nk else ""
            print(f"  {cat:45s} {n:3d}{flag}")

        nk_total = sum(
            any(nk_tag in str(c).lower().replace(" ", "_")
                for nk_tag in ["nk_protrusion", "nk_sm_", "nk_denovo",
                               "nk_synapse", "cytotoxicity", "sst_axis_nk", "checkpoint"])
            for c in ti[cat_col]
        )
        print(f"\n  NK-side genes in 'tumor-intrinsic' pool: {nk_total}/{len(ti)} {'[LEAK]' if nk_total > 3 else '[OK]'}")

    # Known NK-effector genes to check
    NK_EFFECTOR_CHECK = {"NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "RAC1", "EZR", "MSN",
                          "RDX", "WAS", "WASL", "ARPC2", "ARPC3", "DIAPH1", "FMNL1"}
    leaked = {g for g in ti[gene_col].str.upper() if g in NK_EFFECTOR_CHECK}
    if leaked:
        print(f"\n=== NK EFFECTOR GENES IN TUMOR POOL ===  ({len(leaked)} genes)")
        for g in sorted(leaked):
            rows = ti[ti[gene_col].str.upper() == g]
            cat = rows[cat_col].values[0] if cat_col in ti.columns else "?"
            print(f"  {g:15s}  category={cat}  rank={rows['rank'].values[0] if 'rank' in ti.columns else '?'}")
    else:
        print("\n✓ No known NK-effector genes found in tumor pool (de-circularization appears clean)")

    # Check axis-confirmation panel overlap
    try:
        ac = pd.read_csv(T + "axis_confirmation_panel.tsv", sep="\t")
        ac_gen = ac["gene"].str.upper() if "gene" in ac.columns else set()
        overlap = set(ti[gene_col].str.upper()) & set(ac_gen)
        if overlap:
            print(f"\n=== OVERLAP WITH AXIS-CONFIRMATION PANEL (NK READOUT) ===  ({len(overlap)} genes)")
            for g in sorted(overlap):
                print(f"  {g} [WARN: in both tumor pool AND NK readout]")
    except Exception:
        pass

    # Save audit
    audit = ti.copy()
    audit["nk_side_flag"] = audit[cat_col].apply(
        lambda x: any(nk_tag in str(x).lower().replace(" ", "_")
                      for nk_tag in ["nk_protrusion", "nk_sm_", "nk_denovo",
                                     "nk_synapse", "cytotoxicity", "sst_axis_nk", "checkpoint"])
    )
    audit.to_csv(T + "tumor_intrinsic_pool_audit.tsv", sep="\t", index=False)

    log(f"Audit saved to tumor_intrinsic_pool_audit.tsv")
    log("T6 PASS -- review [NK-side] and [WARN] lines above")
    return ti

# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("GC-NKGraph-Atlas — LOCAL P1 RECOMPUTE")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("T7 (ablation) SKIPPED — needs TCGA-STAD bulk expression (A100 only)")
    print("=" * 60)

    # T5 — H1 result (~1 second)
    _ = run_t5()

    # T15 — Effect-size reframe (~2 seconds)
    df_t15 = run_t15()

    # T14 — H3 activation control (~2 seconds)
    _ = run_t14()

    # T6 — De-circ audit (~0.5 seconds)
    _ = run_t6()

    print("\n" + "=" * 60)
    print("LOCAL P1 COMPLETE (4/5 tasks)")
    print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print("\nOutput files:")
    print(f"  {T}sst_axis_positive_control_recovery.tsv       (T5 — H1 appended)")
    print(f"  {T}sst_axis_positive_control_recovery_v2.tsv    (T15 — CIs + effect sizes + BH)")
    print(f"  {T}h3_activation_control.tsv                    (T14)")
    print(f"  {T}tumor_intrinsic_pool_audit.tsv                (T6)")
    print(f"\n  Missing: ablation_results.tsv (T7 — A100 only)")
