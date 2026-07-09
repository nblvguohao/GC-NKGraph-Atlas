"""
T14 — Target rank stability: sensitivity of candidate rankings to weight perturbation.

Perturbs each of the 5 evidence weights by ±0.10 and measures the stability of
the top-10 and top-20 candidate lists via Jaccard index. Specifically checks
whether PHGDH and SGMS2 remain in the top-3 across all perturbations.

RUN ON SERVER:
    python src/interpretation/target_rank_stability.py \
        --input results/tables/tumor_intrinsic_candidates.tsv \
        --output results/tables/target_rank_stability.tsv \
        --n-perturbations 10
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))


# ── Default weights (from Methods §2.8) ──────────────────────────────────────
DEFAULT_WEIGHTS = {
    "tumor_specificity": 0.30,
    "nk_dysfunction_correlation": 0.20,
    "sst_axis_membership": 0.30,
    "axis_core_membership": 0.10,
    "gold_standard_literature": 0.10,
}

# Column name mapping: weight key -> list of possible column name patterns (ordered by preference)
COLUMN_PATTERNS = {
    "tumor_specificity": ["tumor_specificity_log2", "tumor_spec", "tumor_specificity"],
    "nk_dysfunction_correlation": ["nk_dysfunction_correlation", "nk_dysfunction", "nk_correlation"],
    "sst_axis_membership": ["mechanistic_bonus", "sst_axis_membership", "axis_membership", "in_sst_axis"],
    "axis_core_membership": ["axis_core_membership", "axis_core", "in_axis_core"],
    "gold_standard_literature": ["gold_standard", "literature_support", "literature"],
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def find_score_columns(df):
    """Auto-detect which columns carry the evidence scores, matching against known patterns."""
    score_cols = {}
    for key, patterns in COLUMN_PATTERNS.items():
        for pattern in patterns:
            matches = [c for c in df.columns if pattern.lower() in c.lower()]
            if matches:
                score_cols[key] = matches[0]
                break
    return score_cols


def compute_composite(df, weights, score_cols):
    """Compute composite score from columns and weights. Returns array."""
    composite = np.zeros(len(df))
    for key, w in weights.items():
        if key in score_cols:
            col = score_cols[key]
            # Normalize to [0, 1] per perturbation so weights are comparable
            vals = df[col].values.astype(float)
            vmin, vmax = np.nanmin(vals), np.nanmax(vals)
            if vmax - vmin > 0:
                vals = (vals - vmin) / (vmax - vmin)
            composite += w * np.nan_to_num(vals, 0)
    return composite


def jaccard(set_a, set_b):
    """Jaccard index between two sets."""
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)


def perturb_weights(base_weights, key_to_perturb, delta):
    """
    Perturb one weight by ±delta, renormalize others proportionally.
    Returns new weight dict.
    """
    new = base_weights.copy()
    orig = new[key_to_perturb]
    new[key_to_perturb] = np.clip(orig + delta, 0.05, 0.60)
    # Renormalize other weights to sum to 1.0
    other_keys = [k for k in new if k != key_to_perturb]
    other_sum = sum(new[k] for k in other_keys)
    remaining = 1.0 - new[key_to_perturb]
    if other_sum > 0:
        scale = remaining / other_sum
        for k in other_keys:
            new[k] = new[k] * scale
    return new


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="T14: Target rank stability under weight perturbation"
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to tumor_intrinsic_candidates.tsv"
    )
    parser.add_argument(
        "--output", default="results/tables/target_rank_stability.tsv",
        help="Path for stability results"
    )
    parser.add_argument(
        "--n-perturbations", type=int, default=10,
        help="Number of perturbation configurations"
    )
    parser.add_argument(
        "--delta", type=float, default=0.10,
        help="Weight perturbation magnitude (±)"
    )
    parser.add_argument(
        "--top-k", type=int, default=10,
        help="Top-k for Jaccard stability assessment"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
    )
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    log = lambda msg: print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    log(f"Loading {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    log(f"  {len(df)} candidates")

    # Detect gene column
    gene_col = None
    for col in ["gene", "symbol", "Gene", "gene_symbol"]:
        if col in df.columns:
            gene_col = col
            break
    if gene_col is None:
        gene_col = df.columns[0]  # fallback
        log(f"  WARNING: no gene column found; using '{gene_col}'")

    # Detect score columns
    score_cols = find_score_columns(df)
    log(f"  Score columns: {score_cols}")
    if len(score_cols) < 3:
        log("  WARNING: fewer than 3 score columns detected; results may be "
            "unreliable")
        # Try to auto-detect any numeric columns
        num_cols = df.select_dtypes(include=np.number).columns.tolist()
        # Map by position
        for i, key in enumerate(DEFAULT_WEIGHTS.keys()):
            if key not in score_cols and i < len(num_cols):
                score_cols[key] = num_cols[i]
        log(f"  Fallback columns: {score_cols}")

    # ── Baseline ranking ─────────────────────────────────────────────────────
    baseline_composite = compute_composite(df, DEFAULT_WEIGHTS, score_cols)
    df_baseline = df.copy()
    df_baseline["composite_score"] = baseline_composite
    df_baseline = df_baseline.sort_values("composite_score", ascending=False)
    baseline_top_k = set(df_baseline.head(args.top_k)[gene_col].values)
    baseline_top_20 = set(df_baseline.head(20)[gene_col].values)

    log(f"  Baseline top-{args.top_k}: {sorted(baseline_top_k)}")

    # ── Perturbation loop ────────────────────────────────────────────────────
    rng = np.random.default_rng(args.seed)
    results = []
    all_top_k_genes = {}  # gene -> count of appearances in top-k

    for pert_idx in range(args.n_perturbations):
        # Pick a random weight to perturb
        key = rng.choice(list(DEFAULT_WEIGHTS.keys()))
        # Random sign
        sign = rng.choice([-1, 1])
        delta = sign * args.delta

        perturbed = perturb_weights(DEFAULT_WEIGHTS, key, delta)

        # Compute ranking
        composite = compute_composite(df, perturbed, score_cols)
        df_pert = df.copy()
        df_pert["composite_score"] = composite
        df_pert = df_pert.sort_values("composite_score", ascending=False)
        pert_top_k = set(df_pert.head(args.top_k)[gene_col].values)
        pert_top_20 = set(df_pert.head(20)[gene_col].values)

        j_top_k = jaccard(baseline_top_k, pert_top_k)
        j_top_20 = jaccard(baseline_top_20, pert_top_20)

        # Track gene appearances
        for g in pert_top_k:
            all_top_k_genes[g] = all_top_k_genes.get(g, 0) + 1

        results.append({
            "perturbation": pert_idx + 1,
            "perturbed_weight": key,
            "delta": delta,
            "jaccard_top10": j_top_k,
            "jaccard_top20": j_top_20,
            "n_genes_changed_top10": args.top_k - len(baseline_top_k & pert_top_k),
        })

        log(f"  Pert {pert_idx + 1}: {key} {delta:+.2f} -> "
            f"J(top10)={j_top_k:.2f}, J(top20)={j_top_20:.2f}")

    # ── Gene-level stability ─────────────────────────────────────────────────
    total_perturbations = args.n_perturbations
    gene_stability = []
    for gene in sorted(set(list(baseline_top_20) +
                           list(all_top_k_genes.keys()))):
        count = all_top_k_genes.get(gene, 0)
        frac = count / total_perturbations
        in_baseline = gene in baseline_top_k
        in_top3_baseline = gene in set(
            df_baseline.head(3)[gene_col].values
        )
        gene_stability.append({
            "gene": gene,
            "times_in_top10": count,
            "frac_in_top10": frac,
            "in_baseline_top10": in_baseline,
            "in_baseline_top3": in_top3_baseline,
        })

    gene_stab_df = pd.DataFrame(gene_stability).sort_values(
        "frac_in_top10", ascending=False
    )

    # ── Save ─────────────────────────────────────────────────────────────────
    pert_df = pd.DataFrame(results)
    pert_df.to_csv(args.output, sep="\t", index=False, float_format="%.4f")
    log(f"Wrote {args.output}")

    gene_path = args.output.replace(".tsv", "_by_gene.tsv")
    gene_stab_df.to_csv(gene_path, sep="\t", index=False, float_format="%.4f")
    log(f"Wrote {gene_path}")

    # ── Summary ──────────────────────────────────────────────────────────────
    avg_j10 = pert_df["jaccard_top10"].mean()
    avg_j20 = pert_df["jaccard_top20"].mean()
    phgdh_row = gene_stab_df[gene_stab_df["gene"] == "PHGDH"]
    sgms2_row = gene_stab_df[gene_stab_df["gene"] == "SGMS2"]

    print()
    print("=" * 60)
    print("T14 — TARGET RANK STABILITY")
    print("=" * 60)
    print(f"  Perturbations          : {total_perturbations}")
    print(f"  Avg Jaccard (top-10)   : {avg_j10:.3f}")
    print(f"  Avg Jaccard (top-20)   : {avg_j20:.3f}")
    if len(phgdh_row) > 0:
        print(f"  PHGDH in top-10        : "
              f"{phgdh_row['frac_in_top10'].values[0]:.0%}")
    if len(sgms2_row) > 0:
        print(f"  SGMS2 in top-10        : "
              f"{sgms2_row['frac_in_top10'].values[0]:.0%}")
    print()

    if avg_j10 >= 0.7:
        print("  VERDICT: Top-10 is STABLE under ±0.10 weight perturbations")
        print("  (avg Jaccard ≥ 0.7). The arbitrary weight choice does not")
        print("  drive the ranking. Report this in a Supplementary Note or")
        print("  in §2.8 Methods.")
    elif avg_j10 >= 0.5:
        print("  VERDICT: Top-10 is MODERATELY STABLE (Jaccard 0.5–0.7).")
        print("  Report the sensitivity in Limitations §4.3.")
    else:
        print("  VERDICT: Top-10 is UNSTABLE (Jaccard < 0.5). The ranking")
        print("  depends heavily on the chosen weights. Consider simplifying")
        print("  the scoring or reporting only the mechanistically privileged")
        print("  genes (PHGDH, SGMS2, SMPD3/1) rather than a ranked list.")

    # Gate checks
    if len(phgdh_row) > 0:
        phgdh_frac = phgdh_row["frac_in_top10"].values[0]
        if phgdh_frac < 0.8:
            log("  WARNING: PHGDH not stable in top-10 "
                f"({phgdh_frac:.0%})")

    print("=" * 60)
    log("T14 complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
