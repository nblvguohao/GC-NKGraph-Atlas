"""
S1 (nature-reviewer task card): leave-one-sample-out sensitivity analysis for
the H3 (protrusion~cytotoxicity) pseudoreplication-corrected meta-analysis.

Motivation: the DerSimonian-Laird meta-analysis across 9 samples showed
I^2=96% (substantial between-sample heterogeneity), with per-sample r ranging
from 0.009 to 0.560. This checks whether the pooled estimate is driven by one
or two high-r samples, by re-running the meta-analysis with each sample held
out in turn.

Output:
  results/tables/h3_leave_one_sample_out.tsv

Usage:
    python src/topology/h3_leave_one_sample_out.py
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.topology.pseudoreplication_correction import (  # noqa: E402
    per_sample_correlation, dersimonian_laird, fisher_z_inv, log,
)


def main():
    log("=" * 70)
    log("S1: H3 LEAVE-ONE-SAMPLE-OUT SENSITIVITY ANALYSIS")
    log("=" * 70)

    out_dir = "results/tables"
    scores_path = os.path.join(out_dir, "sst_axis_scrna_nk_scores.tsv")
    nk = pd.read_csv(scores_path, sep="\t", index_col=0)
    log(f"Loaded {len(nk)} NK cells from {nk['sample_id'].nunique()} samples")

    z_h3, se_h3, r_h3, p_h3, n_h3 = per_sample_correlation(
        nk, "nk_protrusion_machinery_score", "nk_synapse_cytotoxicity_outcome_score")
    samples = sorted(nk["sample_id"].unique())
    log(f"Samples (k={len(samples)}): {samples}")
    log(f"Per-sample r: {np.array2string(r_h3, precision=4)}")

    meta_full = dersimonian_laird(z_h3, se_h3)
    log(f"\nFull meta-analysis (k={len(samples)}): "
        f"r={meta_full['pooled_r']:.4f}, p={meta_full['p_value']:.4e}, "
        f"I2={meta_full['I2']:.1f}%")

    rows = [{
        "excluded_sample": "none (full k=9)",
        "k_samples": len(samples),
        "pooled_r": round(meta_full["pooled_r"], 4),
        "pooled_p": meta_full["p_value"],
        "I2_pct": round(meta_full["I2"], 1),
        "ci_lower": round(fisher_z_inv(meta_full["pooled_z"] - 1.96 * meta_full["pooled_se"]), 4),
        "ci_upper": round(fisher_z_inv(meta_full["pooled_z"] + 1.96 * meta_full["pooled_se"]), 4),
    }]

    for i, excluded in enumerate(samples):
        keep = np.array([s != excluded for s in samples])
        z_loo, se_loo = z_h3[keep], se_h3[keep]
        meta_loo = dersimonian_laird(z_loo, se_loo)
        ci_lo = fisher_z_inv(meta_loo["pooled_z"] - 1.96 * meta_loo["pooled_se"])
        ci_hi = fisher_z_inv(meta_loo["pooled_z"] + 1.96 * meta_loo["pooled_se"])
        rows.append({
            "excluded_sample": excluded,
            "k_samples": len(z_loo),
            "pooled_r": round(meta_loo["pooled_r"], 4),
            "pooled_p": meta_loo["p_value"],
            "I2_pct": round(meta_loo["I2"], 1),
            "ci_lower": round(ci_lo, 4),
            "ci_upper": round(ci_hi, 4),
        })
        log(f"  Excluding {excluded} (own r={r_h3[i]:.3f}): "
            f"pooled r={meta_loo['pooled_r']:.4f}, p={meta_loo['p_value']:.2e}, "
            f"I2={meta_loo['I2']:.1f}%, 95% CI [{ci_lo:.3f}, {ci_hi:.3f}]")

    df = pd.DataFrame(rows)
    out_path = os.path.join(out_dir, "h3_leave_one_sample_out.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    log(f"\nSaved: {out_path}")

    non_full = df[df["excluded_sample"] != "none (full k=9)"]
    log(f"\nRange of pooled r across leave-one-out runs: "
        f"[{non_full['pooled_r'].min():.4f}, {non_full['pooled_r'].max():.4f}]")
    log(f"All leave-one-out CIs exclude 0: {(non_full['ci_lower'] > 0).all()}")
    log("\nS1 COMPLETE!")


if __name__ == "__main__":
    main()
