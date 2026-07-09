"""
T11 — Gene embedding analysis: demonstrate GNN interpretability.

Produces concrete evidence that the GNN's learned gene embeddings capture
mechanism-relevant structure that tabular baselines do not provide:

  1. UMAP projection of gene embeddings, colored by SST-axis module membership
  2. Embedding distance vs scRNA co-expression correlation (Spearman)
  3. Top-5 closest gene pairs across axis modules (qualitative inspection)

RUN ON SERVER (requires trained model checkpoint):
    python src/interpretation/gene_embedding_analysis.py \
        --embeddings results/model/gene_embeddings.npy \
        --gene-list data/processed/graph/gene_nodes.tsv \
        --sst-config configs/sst_axis_config.yaml \
        --output-dir results/

Alternatively, if embeddings are stored differently (e.g., in a .pt file or
embedded in a model checkpoint), adapt the --embeddings flag accordingly. If
no checkpoint is available, the gene embeddings can be re-extracted from the
trained HGT model's gene node representations.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_PROJECT_ROOT))


# ── SST module definitions (same as sst_axis_config.yaml) ────────────────────
SST_MODULES = {
    "tumor_serine": [
        "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2",
        "MTHFD1", "MTHFD2", "MTHFD1L", "SLC1A4", "SLC1A5",
    ],
    "nk_sm_synthesis": ["SGMS1", "SGMS2"],
    "nk_sm_catabolism": ["SMPD1", "SMPD2", "SMPD3", "SMPD4"],
    "nk_denovo_sphingolipid": [
        "SPTLC1", "SPTLC2", "SPTLC3", "SPTSSA",
        "CERS2", "CERS4", "CERS5", "CERS6", "DEGS1",
    ],
    "nk_protrusion_machinery": [
        "EZR", "MSN", "RDX",
        "ACTR2", "ACTR3", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5",
        "WAS", "WASL", "WASF1", "WASF2", "WASF3", "WIPF1",
        "CDC42", "RAC1", "RHOA",
        "DIAPH1", "DIAPH3", "FMNL1",
        "BAIAP2", "PACSIN2",
    ],
    "nk_cytotoxicity": [
        "NKG7", "GNLY", "GZMB", "PRF1", "IFNG",
        "LCP2", "LAT", "VAV1", "TLN1", "ITGAL", "ITGB2",
    ],
    "checkpoint": ["HAVCR2"],
}

MODULE_COLORS = {
    "tumor_serine": "#D55E00",
    "nk_sm_synthesis": "#0072B2",
    "nk_sm_catabolism": "#009E73",
    "nk_denovo_sphingolipid": "#CC79A7",
    "nk_protrusion_machinery": "#F0E442",
    "nk_cytotoxicity": "#E69F00",
    "checkpoint": "#000000",
    "other": "#AAAAAA",
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def build_module_map(gene_list, modules):
    """Return dict: gene -> module_name (or 'other')."""
    gene_to_module = {}
    for gene in gene_list:
        assigned = "other"
        for mod_name, mod_genes in modules.items():
            if gene in mod_genes:
                assigned = mod_name
                break
        gene_to_module[gene] = assigned
    return gene_to_module


def cosine_distance_matrix(embeddings):
    """Pairwise cosine distance matrix from embedding array (n_genes × d)."""
    # Normalize to unit vectors
    norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    emb_norm = embeddings / norms
    # Cosine similarity -> distance
    sim = emb_norm @ emb_norm.T
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    return dist


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="T11: Gene embedding analysis — UMAP + distance vs coexpression"
    )
    parser.add_argument(
        "--embeddings", required=True,
        help="Path to gene_embeddings.npy (n_genes × d_dim)"
    )
    parser.add_argument(
        "--gene-list", required=True,
        help="Path to gene_nodes.tsv (must have 'gene' column)"
    )
    parser.add_argument(
        "--sst-config", default="configs/sst_axis_config.yaml",
        help="Path to sst_axis_config.yaml (for module definitions)"
    )
    parser.add_argument(
        "--scRNA-expr",
        help="Optional: scRNA expression matrix for co-expression comparison "
             "(genes × cells .tsv). If not provided, skip D2."
    )
    parser.add_argument(
        "--output-dir", default="results",
        help="Output directory"
    )
    parser.add_argument(
        "--random-seed", type=int, default=42,
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "tables"), exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "figures"), exist_ok=True)

    log = lambda msg: print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    log(f"Loading embeddings from {args.embeddings}")
    embeddings = np.load(args.embeddings)  # shape: (n_genes, d)
    n_genes, d_dim = embeddings.shape
    log(f"  {n_genes} genes, {d_dim}-dimensional embeddings")

    log(f"Loading gene list from {args.gene_list}")
    genes_df = pd.read_csv(args.gene_list, sep="\t")
    if "gene" not in genes_df.columns:
        # Try common alternatives
        for col in ["symbol", "gene_symbol", "name", "Gene"]:
            if col in genes_df.columns:
                genes_df.rename(columns={col: "gene"}, inplace=True)
                break
        else:
            genes_df["gene"] = genes_df.iloc[:, 0]  # first column as fallback
    gene_names = genes_df["gene"].tolist()

    if len(gene_names) != n_genes:
        log(f"  WARNING: gene list length ({len(gene_names)}) != "
            f"embedding rows ({n_genes}). Truncating to min.")
        min_n = min(len(gene_names), n_genes)
        gene_names = gene_names[:min_n]
        embeddings = embeddings[:min_n]

    # ── Module assignment ────────────────────────────────────────────────────
    gene_to_module = build_module_map(gene_names, SST_MODULES)
    module_labels = [gene_to_module.get(g, "other") for g in gene_names]
    n_assigned = sum(1 for m in module_labels if m != "other")
    log(f"  {n_assigned}/{len(gene_names)} genes assigned to SST modules")

    # ── D1: Embedding distance matrix ────────────────────────────────────────
    log("Computing pairwise cosine distances...")
    dist_matrix = cosine_distance_matrix(embeddings)

    # Within-module vs between-module distances
    within_dists, between_dists = [], []
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            d = dist_matrix[i, j]
            if module_labels[i] != "other" and module_labels[i] == module_labels[j]:
                within_dists.append(d)
            elif module_labels[i] != "other" and module_labels[j] != "other":
                between_dists.append(d)

    within_mean = np.mean(within_dists) if within_dists else np.nan
    between_mean = np.mean(between_dists) if between_dists else np.nan
    log(f"  Mean within-module distance: {within_mean:.4f}")
    log(f"  Mean between-module distance: {between_mean:.4f}")

    # Statistical test: are within-module pairs closer than between-module?
    if within_dists and between_dists:
        u_stat, p_val = stats.mannwhitneyu(
            within_dists, between_dists, alternative="less"
        )
        log(f"  Mann-Whitney U: within < between, p={p_val:.2e}")
    else:
        u_stat, p_val = np.nan, np.nan

    # ── D2: Embedding distance vs co-expression ──────────────────────────────
    spearman_r, spearman_p = np.nan, np.nan
    if args.scRNA_expr and os.path.exists(args.scRNA_expr):
        log(f"Loading scRNA expression from {args.scRNA_expr}")
        expr_df = pd.read_csv(args.scRNA_expr, sep="\t", index_col=0)
        # Subset to genes in embedding
        common_genes = [g for g in gene_names if g in expr_df.index]
        common_set = set(common_genes)
        log(f"  {len(common_genes)} genes in both embedding and scRNA")

        # Compute co-expression (Pearson across cells) for each gene pair
        embedding_pairs = []
        coexp_corrs = []
        for i in range(n_genes):
            gi = gene_names[i]
            if gi not in common_set:
                continue
            for j in range(i + 1, n_genes):
                gj = gene_names[j]
                if gj not in common_set:
                    continue
                if module_labels[i] == "other" and module_labels[j] == "other":
                    continue  # skip non-SST pairs for efficiency
                emb_dist = dist_matrix[i, j]
                coexp_r, _ = stats.pearsonr(
                    expr_df.loc[gi].values, expr_df.loc[gj].values
                )
                if not np.isfinite(coexp_r):
                    continue
                embedding_pairs.append(emb_dist)
                coexp_corrs.append(coexp_r)

        if len(embedding_pairs) > 10:
            spearman_r, spearman_p = stats.spearmanr(embedding_pairs, coexp_corrs)
            log(f"  Embedding dist vs coexpression: ρ={spearman_r:.4f}, "
                f"p={spearman_p:.2e}")
        else:
            log("  Too few pairs for Spearman; skipping")
    else:
        log("No scRNA expression matrix provided; skipping D2")

    # ── D3: Top-5 closest cross-module gene pairs ────────────────────────────
    log("Finding top cross-module gene pairs...")
    pairs = []
    for i in range(n_genes):
        if module_labels[i] == "other":
            continue
        for j in range(i + 1, n_genes):
            if module_labels[j] == "other":
                continue
            if module_labels[i] == module_labels[j]:
                continue  # within-module only
            pairs.append({
                "gene_A": gene_names[i],
                "module_A": module_labels[i],
                "gene_B": gene_names[j],
                "module_B": module_labels[j],
                "cosine_distance": dist_matrix[i, j],
            })

    pairs_df = pd.DataFrame(pairs).sort_values("cosine_distance")
    top_pairs = pairs_df.head(50)

    # Special interest: tumor_serine <-> nk_sm pairs
    tumor_nk = pairs_df[
        ((pairs_df["module_A"] == "tumor_serine") &
         pairs_df["module_B"].isin(["nk_sm_synthesis", "nk_sm_catabolism"])) |
        ((pairs_df["module_B"] == "tumor_serine") &
         pairs_df["module_A"].isin(["nk_sm_synthesis", "nk_sm_catabolism"]))
    ]
    log(f"  Tumor-serine — NK-SM pairs within top-100: "
        f"{len(tumor_nk.head(100))}")

    # ── UMAP projection ──────────────────────────────────────────────────────
    umap_coords = None
    try:
        import umap
        log("Running UMAP on gene embeddings...")
        reducer = umap.UMAP(
            n_components=2, metric="cosine", random_state=args.random_seed,
            n_neighbors=min(30, n_genes - 1),
        )
        umap_coords = reducer.fit_transform(embeddings)
        log(f"  UMAP done: shape {umap_coords.shape}")
    except ImportError:
        log("umap-learn not installed; skipping UMAP")
    except Exception as e:
        log(f"UMAP failed: {e}; skipping")

    # ── Save tables ──────────────────────────────────────────────────────────
    tables_dir = os.path.join(args.output_dir, "tables")

    # Module distance summary
    dist_summary = pd.DataFrame([{
        "mean_within_module_distance": within_mean,
        "mean_between_module_distance": between_mean,
        "n_within_pairs": len(within_dists),
        "n_between_pairs": len(between_dists),
        "mw_u_statistic": u_stat,
        "mw_p_value": p_val,
    }])
    dist_summary.to_csv(
        os.path.join(tables_dir, "embedding_module_distances.tsv"),
        sep="\t", index=False
    )

    # Embedding vs coexpression
    eco_df = pd.DataFrame([{
        "spearman_rho": spearman_r,
        "spearman_p": spearman_p,
        "n_gene_pairs": len(embedding_pairs) if args.scRNA_expr else 0,
    }])
    eco_df.to_csv(
        os.path.join(tables_dir, "embedding_distance_vs_coexpression.tsv"),
        sep="\t", index=False
    )

    # Top cross-module pairs
    top_pairs.to_csv(
        os.path.join(tables_dir, "embedding_top_cross_module_pairs.tsv"),
        sep="\t", index=False
    )

    # UMAP coordinates
    if umap_coords is not None:
        umap_df = pd.DataFrame({
            "gene": gene_names,
            "module": module_labels,
            "umap_1": umap_coords[:, 0],
            "umap_2": umap_coords[:, 1],
        })
        umap_df.to_csv(
            os.path.join(tables_dir, "gene_embedding_umap.tsv"),
            sep="\t", index=False
        )

    log("All tables written.")

    # ── Plot ─────────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2 if umap_coords is not None else 1,
                                 figsize=(14 if umap_coords is not None else 7, 6))
        if umap_coords is None:
            axes = [axes]

        # Panel A: module distance comparison
        ax = axes[0]
        categories = ["Within\nmodule", "Between\nmodules"]
        means = [within_mean, between_mean]
        colors = ["#0072B2", "#D55E00"]
        bars = ax.bar(categories, means, color=colors, edgecolor="white", width=0.5)
        for bar, val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                    f"{val:.4f}", ha="center", fontsize=9)
        ax.set_ylabel("Mean cosine distance")
        ax.set_title("A. Gene embedding distances\n(within vs between SST modules)")
        if not np.isnan(p_val):
            ax.text(0.5, 0.95,
                    f"Mann-Whitney p={p_val:.1e}\n(within < between)",
                    transform=ax.transAxes, ha="center", va="top",
                    fontsize=8, style="italic")

        # Panel B: UMAP
        if umap_coords is not None:
            ax = axes[1]
            for mod in sorted(set(module_labels)):
                mask = np.array([m == mod for m in module_labels])
                color = MODULE_COLORS.get(mod, "#AAAAAA")
                alpha = 0.7 if mod != "other" else 0.15
                zorder = 2 if mod != "other" else 1
                size = 8 if mod != "other" else 2
                ax.scatter(
                    umap_coords[mask, 0], umap_coords[mask, 1],
                    c=color, label=mod, alpha=alpha, s=size, zorder=zorder,
                    edgecolors="none",
                )
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.set_title("B. Gene embedding UMAP\n(colored by SST module)")
            ax.legend(fontsize=6, markerscale=1.5, loc="lower left",
                      bbox_to_anchor=(1.01, 0))

        fig.suptitle(
            f"Gene embedding interpretability (n={n_genes} genes, d={d_dim})",
            fontsize=11, fontweight="bold"
        )
        fig.tight_layout()
        fig_path = os.path.join(args.output_dir, "figures",
                                "fig6_gene_embedding_analysis.pdf")
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        log(f"Wrote {fig_path}")
        plt.close(fig)

    except ImportError:
        log("matplotlib not available; skipping plot")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("T11 — GENE EMBEDDING INTERPRETABILITY")
    print("=" * 60)
    print(f"  Genes in embedding        : {n_genes}")
    print(f"  Assigned to SST modules   : {n_assigned}")
    print(f"  Within-module distance    : {within_mean:.4f}")
    print(f"  Between-module distance   : {between_mean:.4f}")
    print(f"  Within < between (MWU)    : p={p_val:.2e}")
    if not np.isnan(spearman_r):
        print(f"  Embedding vs coexpression : ρ={spearman_r:.4f}, p={spearman_p:.2e}")
    print()
    if within_mean < between_mean and p_val < 0.05:
        print("  VERDICT: SST-related genes are closer in embedding space")
        print("  than expected by chance. The GNN embedding captures")
        print("  mechanism-relevant structure. Use this result in §4.2.")
    else:
        print("  VERDICT: No significant clustering by SST module in")
        print("  embedding space. The 'added interpretability' claim")
        print("  should be qualified accordingly in the manuscript.")
    print("=" * 60)

    log("T11 complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
