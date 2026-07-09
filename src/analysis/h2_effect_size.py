"""
T9 — H2 effect size quantification: bootstrap CI + down-sampling analysis.

The H2 correlation (SM-balance -> protrusion) is r=+0.030 at n=8,310 single NK
cells. At this sample size, virtually any non-zero correlation is significant.
This script quantifies:
  1. Bootstrap 95% CI on r (how tight is the estimate?)
  2. Variance explained (r²)
  3. Down-sampling curve: at what n does the signal become undetectable?

RUN ON SERVER (requires scRNA SST scores):
    python src/analysis/h2_effect_size.py \
        --input results/tables/sst_axis_scores_single_cell.tsv \
        --output-dir results/tables/ \
        --output-fig results/figures/figS_h2_downsample.pdf
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))

# ── Helpers ──────────────────────────────────────────────────────────────────

def pearson_r(x, y):
    """Return r and p, skipping NaN."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return np.nan, np.nan
    return stats.pearsonr(x[mask], y[mask])


def bootstrap_ci(x, y, n_bootstrap=10000, ci=95, seed=42):
    """Percentile bootstrap CI for Pearson r."""
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        r_val, _ = pearson_r(x[idx], y[idx])
        boots[i] = r_val
    alpha = (100 - ci) / 2
    lo, hi = np.percentile(boots, [alpha, 100 - alpha])
    return lo, hi, boots


def downsample_curve(x, y, min_n=50, step=200, n_repeats=50, seed=42):
    """
    Down-sample from full n to min_n in steps, repeating n_repeats times
    per level. Returns DataFrame with columns: n, r_mean, r_std, p_mean.
    """
    rng = np.random.default_rng(seed)
    n_total = len(x)
    results = []
    ns = list(range(n_total, min_n - 1, -step))
    if ns[-1] != min_n:
        ns.append(min_n)

    for n in ns:
        rs, ps = [], []
        for _ in range(n_repeats):
            idx = rng.choice(n_total, size=n, replace=False)
            r_val, p_val = pearson_r(x[idx], y[idx])
            rs.append(r_val)
            ps.append(p_val)
        results.append({
            "n": n,
            "r_mean": np.mean(rs),
            "r_std": np.std(rs),
            "r_median": np.median(rs),
            "p_mean": np.mean(ps),
            "p_median": np.median(ps),
            "frac_significant": np.mean([p < 0.05 for p in ps]),
        })

    return pd.DataFrame(results)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="T9: H2 effect size — bootstrap CI + down-sampling curve"
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to sst_axis_scores_single_cell.tsv"
    )
    parser.add_argument(
        "--output-dir", default="results/tables",
        help="Directory for output tables"
    )
    parser.add_argument(
        "--output-fig", default="results/figures/figS_h2_downsample.pdf",
        help="Path for down-sampling plot (PDF)"
    )
    parser.add_argument(
        "--n-bootstrap", type=int, default=10000,
        help="Number of bootstrap resamples"
    )
    parser.add_argument(
        "--min-n", type=int, default=100,
        help="Smallest n for down-sampling"
    )
    parser.add_argument(
        "--step", type=int, default=200,
        help="Step size for down-sampling"
    )
    args = parser.parse_args()

    # Ensure output dir
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output_fig), exist_ok=True)

    log = lambda msg: print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    log(f"Loading {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    # Detect columns
    if "nk_sm_balance_score" in df.columns:
        sm_bal = df["nk_sm_balance_score"].values
    else:
        syn = df.get("nk_sm_synthesis_score", pd.Series(0, index=df.index)).values
        cat = df.get("nk_sm_catabolism_score", pd.Series(0, index=df.index)).values
        sm_bal = syn - cat

    if "nk_protrusion_machinery_score" not in df.columns:
        log("ERROR: nk_protrusion_machinery_score column not found")
        sys.exit(1)

    prot = df["nk_protrusion_machinery_score"].values
    n_total = len(df)

    log(f"Loaded {n_total} cells")

    # ── Compute observed r ───────────────────────────────────────────────────
    r_obs, p_obs = pearson_r(sm_bal, prot)
    r2_obs = r_obs ** 2

    log(f"H2 observed: r={r_obs:.4f}, p={p_obs:.2e}, r²={r2_obs:.6f} "
        f"({r2_obs*100:.3f}% variance explained)")

    # ── Bootstrap CI ─────────────────────────────────────────────────────────
    log(f"Running bootstrap (n={args.n_bootstrap})...")
    ci_lo, ci_hi, boots = bootstrap_ci(
        sm_bal, prot, n_bootstrap=args.n_bootstrap, seed=42
    )
    log(f"95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    log(f"Bootstrap mean r = {boots.mean():.4f}")

    # ── Down-sample curve ────────────────────────────────────────────────────
    log(f"Running down-sampling curve ({n_total} -> {args.min_n}, "
        f"step={args.step})...")
    ds_df = downsample_curve(sm_bal, prot, min_n=args.min_n, step=args.step)

    # Find the n at which signal loses significance (p>=0.05)
    sig_rows = ds_df[ds_df["p_mean"] < 0.05]
    if len(sig_rows) > 0:
        n_critical = sig_rows["n"].min()
        log(f"Signal loses significance at n ≈ {n_critical}")
    else:
        n_critical = ds_df["n"].max()
        log(f"Signal significant at all tested n (down to {args.min_n})")

    # ── Save tables ──────────────────────────────────────────────────────────
    ci_df = pd.DataFrame([{
        "test": "H2_sm_balance_vs_protrusion",
        "n_cells": n_total,
        "r_observed": r_obs,
        "p_observed": p_obs,
        "r_squared": r2_obs,
        "variance_explained_pct": r2_obs * 100,
        "ci_lower_95": ci_lo,
        "ci_upper_95": ci_hi,
        "bootstrap_mean_r": boots.mean(),
        "n_critical_p05": n_critical,
    }])
    ci_path = os.path.join(args.output_dir, "h2_bootstrap_ci.tsv")
    ci_df.to_csv(ci_path, sep="\t", index=False)
    log(f"Wrote {ci_path}")

    ds_path = os.path.join(args.output_dir, "h2_downsample_curve.tsv")
    ds_df.to_csv(ds_path, sep="\t", index=False)
    log(f"Wrote {ds_path}")

    # ── Summary assessment ───────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("T9 — H2 EFFECT SIZE ASSESSMENT")
    print("=" * 60)
    print(f"  Observed r      : {r_obs:.4f}")
    print(f"  Bootstrap 95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
    print(f"  Variance explained: {r2_obs*100:.3f}%")
    print(f"  CI crosses zero : {ci_lo <= 0 <= ci_hi}")
    print(f"  n for p>=0.05   : {n_critical}")
    print()
    if ci_lo > 0:
        print("  VERDICT: H2 is statistically distinguishable from zero")
        print("  but the effect is TINY (r² < 0.1%). The language in the")
        print("  manuscript must qualify this as 'weak but correctly signed'")
        print("  rather than 'recovered.'")
    else:
        print("  VERDICT: H2 CI crosses zero — the effect is not reliably")
        print("  distinguishable from zero. Manuscript should state this.")
    print("=" * 60)

    # ── Plot (matplotlib, optional) ──────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Panel A: bootstrap distribution
        ax = axes[0]
        ax.hist(boots, bins=80, color="#0072B2", edgecolor="white", alpha=0.85)
        ax.axvline(r_obs, color="#D55E00", linewidth=2, linestyle="-",
                   label=f"observed r={r_obs:.4f}")
        ax.axvline(ci_lo, color="#009E73", linewidth=1.5, linestyle="--",
                   label=f"95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
        ax.axvline(ci_hi, color="#009E73", linewidth=1.5, linestyle="--")
        ax.axvline(0, color="gray", linewidth=0.5, linestyle=":")
        ax.set_xlabel("Pearson r")
        ax.set_ylabel("Bootstrap count")
        ax.set_title("A. Bootstrap distribution of H2 (r=0.030)")
        ax.legend(fontsize=8)

        # Panel B: down-sampling curve
        ax = axes[1]
        ax.errorbar(ds_df["n"], ds_df["r_mean"], yerr=ds_df["r_std"],
                    fmt="o-", color="#0072B2", markersize=3, capsize=2,
                    label="r ± 1 SD")
        ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
        ax.axhline(r_obs, color="#D55E00", linewidth=1, linestyle="--",
                   alpha=0.5, label=f"full-sample r={r_obs:.4f}")
        ax.axvline(n_critical, color="#CC79A7", linewidth=1, linestyle=":",
                   label=f"p≥0.05 at n≈{n_critical}")
        ax.set_xlabel("Number of NK cells")
        ax.set_ylabel("Pearson r (mean ± SD)")
        ax.set_title("B. Down-sampling curve (50 repeats / level)")
        ax.invert_xaxis()
        ax.legend(fontsize=7)

        fig.suptitle(
            f"H2 effect size: SM-balance → protrusion "
            f"(r={r_obs:.4f}, r²={r2_obs:.4f}, n={n_total})",
            fontsize=11, fontweight="bold"
        )
        fig.tight_layout()
        fig.savefig(args.output_fig, dpi=300, bbox_inches="tight")
        log(f"Wrote {args.output_fig}")
        plt.close(fig)

    except ImportError:
        log("matplotlib not available; skipping figure generation")

    # ── Gate checks ──────────────────────────────────────────────────────────
    if ci_lo <= 0 and ci_hi > 0:
        log("WARNING: Bootstrap CI crosses zero — H2 is NOT reliably > 0")
    if r2_obs < 0.001:
        log("WARNING: r² < 0.1% — H2 explains negligible variance; "
            "qualify as 'weak' in manuscript")

    log("T9 complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
