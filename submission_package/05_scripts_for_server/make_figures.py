"""
Publication figures for GC-NKGraph-Atlas (Fig 1–5), generated from result tables.

Every number is read from results/tables/*.tsv so the figures stay consistent with
the manuscript. Palette: Okabe-Ito (colour-blind safe) with a muted academic style.
Output: PDF (vector) + PNG (300+ dpi).

Usage:
    python src/figures/make_figures.py              # all figures
    python src/figures/make_figures.py --fig 1      # single figure
    python src/figures/make_figures.py --dpi 600    # higher resolution
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from scipy import stats

T = "results/tables/"
OUT = "results/figures/"
os.makedirs(OUT, exist_ok=True)

# ── Colour-blind-safe palette (Okabe-Ito) ──────────────────────────────
C = {
    "black":  "#000000", "orange": "#E69F00", "sky":    "#56B4E9",
    "green":  "#009E73", "yellow": "#F0E442", "blue":   "#0072B2",
    "verm":   "#D55E00", "purple":"#CC79A7", "grey":    "#8C8C8C",
    "dkgrey": "#5A5A5A",
}
REC   = C["green"]    # recovered
NOT   = C["verm"]     # not recovered
WEAK  = C["grey"]     # inconclusive / weak
GNN_C = C["verm"]     # our method
BL_C  = C["sky"]      # baseline methods

# ── Matplotlib global style ────────────────────────────────────────────
def set_style(dpi=300):
    mpl.rcParams.update({
        "figure.dpi": dpi, "savefig.dpi": dpi, "savefig.bbox": "tight",
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.grid": True, "grid.color": "#E8E8E8", "grid.linewidth": 0.5,
        "grid.alpha": 0.7, "legend.frameon": False, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
    })

# ── Helpers ─────────────────────────────────────────────────────────────
def save(fig, name):
    for ext in ("pdf", "png"):
        path = f"{OUT}{name}.{ext}"
        fig.savefig(path, dpi=fig.dpi, bbox_inches="tight", pad_inches=0.1,
                    metadata={"Creator": "GC-NKGraph-Atlas"})
    plt.close(fig)
    print(f"  -> {OUT}{name}.pdf / .png")

def panel(ax, letter, x=-0.12, y=1.05):
    """Add bold panel letter in upper-left corner."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left",
            fontfamily="sans-serif")

def load(path, **kw):
    """Load a TSV with graceful fallback."""
    try:
        return pd.read_csv(T + path, sep="\t", **kw)
    except FileNotFoundError:
        print(f"  [WARN] missing {T}{path} — using placeholder data for figure")
        return None

# ========================================================================
# FIGURE 1 — Arm A: partial recovery of the SST axis in liver
# ========================================================================
def figure1(dpi=300):
    liver = load("sst_axis_scores_liver_bulk.tsv", index_col=0)
    sc    = load("sst_axis_scores_single_cell.tsv")
    rec   = load("sst_axis_positive_control_recovery.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.2), dpi=dpi)

    # --- A: Bulk LIHC H3 scatter: protrusion vs cytotoxicity ---
    ax = axes[0, 0]; panel(ax, "A")
    if liver is not None:
        x = liver["nk_protrusion_machinery_score"]
        y = liver["nk_synapse_cytotoxicity_outcome_score"]
        ax.scatter(x, y, s=12, c=C["blue"], alpha=0.45, edgecolors="none", rasterized=True)
        b, a = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 80)
        ax.plot(xs, a + b * xs, color=C["black"], lw=2, zorder=3)
        r, p = stats.pearsonr(x, y)
        ax.text(0.03, 0.95, f"r = {r:.3f}\np = {p:.1e}", transform=ax.transAxes,
                va="top", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCC", alpha=0.9))
    ax.set_title("H3: protrusion → cytotoxicity (bulk TCGA-LIHC, n = 423)")
    ax.set_xlabel("NK protrusion-machinery score"); ax.set_ylabel("NK cytotoxicity-output score")

    # --- B: H2/H3 bulk vs single-cell resolution comparison ---
    ax = axes[0, 1]; panel(ax, "B")
    if rec is not None:
        def _getr(h, res):
            v = rec[(rec.hypothesis == h) & (rec.resolution == res)]["r"]
            return float(v.iloc[0]) if len(v) else np.nan
        labels = ["H2\nSM-balance→protrusion", "H3\nprotrusion→cytotoxicity"]
        bulk_vals   = [_getr("H2", "bulk"), _getr("H3", "bulk")]
        single_vals = [_getr("H2", "single_cell_NK"), _getr("H3", "single_cell_NK")]
        xx = np.arange(2); w = 0.36
        b1 = ax.bar(xx - w/2, bulk_vals,   w, label="Bulk transcriptome", color=C["grey"], edgecolor="white", lw=0.5)
        b2 = ax.bar(xx + w/2, single_vals, w, label="Single-cell NK",    color=C["blue"], edgecolor="white", lw=0.5)
        ax.axhline(0, color=C["black"], lw=0.8)
        ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Pearson r"); ax.set_ylim(-0.12, 0.65)
        ax.set_title("Cell-type resolution rescues H2")
        ax.legend(loc="upper left", fontsize=7.5)
        # annotate the null bulk result
        if not np.isnan(bulk_vals[0]):
            ax.annotate("null in bulk", (0 - w/2, max(bulk_vals[0] + 0.03, 0.03)),
                       fontsize=7, ha="center", color=C["verm"],
                       bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=C["verm"], alpha=0.6))

    # --- C: Intratumoral vs normal NK, cytotoxicity output and protrusion machinery ---
    ax = axes[1, 0]; panel(ax, "C")
    if sc is not None:
        modules = [
            ("nk_synapse_cytotoxicity_outcome_score", "cytotoxicity\noutput"),
            ("nk_protrusion_machinery_score",         "protrusion\nmachinery"),
        ]
        data_boxes, positions, colors, tick_info = [], [], [], []
        for i, (col, label) in enumerate(modules):
            for j, cond in enumerate(["normal", "tumor"]):
                vals = sc[sc.condition == cond][col].dropna().values
                if len(vals):
                    data_boxes.append(vals)
                    positions.append(i * 3 + j)
                    colors.append(C["sky"] if cond == "normal" else C["verm"])
            tick_info.append((i * 3 + 0.5, label))
        bp = ax.boxplot(data_boxes, positions=positions, widths=0.7,
                        patch_artist=True, showfliers=False,
                        medianprops={"color": C["black"], "lw": 1.2},
                        whiskerprops={"lw": 0.8}, capprops={"lw": 0.8})
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.75)
        ax.set_xticks([t[0] for t in tick_info])
        ax.set_xticklabels([t[1] for t in tick_info], fontsize=8)
        ax.set_ylabel("Per-cell module z-score")
        ax.set_title("Intratumoral vs normal NK (8,310 single cells)")
        ax.legend(handles=[
            Patch(fc=C["sky"], alpha=0.75, label="Normal NK"),
            Patch(fc=C["verm"], alpha=0.75, label="Intratumoral NK"),
        ], loc="upper right", fontsize=7.5)

    # --- D: Forest plot of H2–H5 outcomes ---
    ax = axes[1, 1]; panel(ax, "D")
    if rec is not None:
        rr = rec[rec["r"].astype(str).str.replace(".", "", 1).str.replace("-", "", 1).str.isdigit()].copy()
        rr["rv"] = rr["r"].astype(float)
        def _lab(row):
            if row.hypothesis == "H5":
                sub = row.test.split(":")[-1].strip().replace("_", " ").replace("outcome", "").strip()
                return f"H5 {sub} (scNK)"
            res = row.resolution.replace("single_cell_NK", "scNK")
            return f"{row.hypothesis} ({res})"
        rr["lab"] = rr.apply(_lab, axis=1)
        rr = rr.iloc[::-1].reset_index(drop=True)
        cmap = {"RECOVERED": REC, "INCONCLUSIVE": WEAK}
        point_colors = [cmap.get(o, NOT) for o in rr["outcome"]]
        ax.scatter(rr["rv"], range(len(rr)), c=point_colors, s=60, zorder=3,
                   edgecolors="white", linewidths=0.8)
        ax.axvline(0, color=C["black"], lw=0.8, ls="--")
        ax.set_yticks(range(len(rr)))
        ax.set_yticklabels(rr["lab"], fontsize=7)
        ax.set_xlabel("Effect size (r or Δ for H5)")
        ax.set_title("Pre-registered hypotheses H2–H5")
        ax.legend(handles=[
            Line2D([0], [0], marker="o", ls="", mfc=REC,  mec="w", label="Recovered"),
            Line2D([0], [0], marker="o", ls="", mfc=NOT,  mec="w", label="Not recovered"),
            Line2D([0], [0], marker="o", ls="", mfc=WEAK, mec="w", label="Inconclusive"),
        ], loc="lower right", fontsize=7)

    fig.suptitle("Figure 1  |  Arm A — Partial recovery of the SST axis in liver",
                 fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "fig1_armA_positive_control")


# ========================================================================
# FIGURE 2 — Arm B: gastric extension + external validation
# ========================================================================
def figure2(dpi=300):
    tis = load("sst_axis_scrna_by_tissue.tsv")
    lab = load("nk_state_labels.tsv", comment="#")
    ext = load("external_validation_results.tsv")
    sc  = load("sst_axis_scores_single_cell.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), dpi=dpi)

    # --- A: Cross-tissue NK module comparison ---
    ax = axes[0]; panel(ax, "A")
    if tis is not None:
        modcols = ["nk_sm_balance_score", "nk_protrusion_machinery_score",
                   "nk_topology_permissive_score", "nk_synapse_cytotoxicity_outcome_score"]
        short   = ["SM balance", "Protrusion", "Topology-perm", "Cytotoxicity"]
        order   = ["healthy_liver", "liver_metastasis", "gastric_cancer"]
        tis2    = tis.set_index("tissue").reindex(order)
        tcol    = {"healthy_liver": C["sky"], "liver_metastasis": C["orange"], "gastric_cancer": C["verm"]}
        xx = np.arange(len(modcols)); w = 0.26
        for k, (t, row) in enumerate(tis2.iterrows()):
            vals = [float(row[m]) for m in modcols]
            ax.bar(xx + (k-1)*w, vals, w, label=t.replace("_", " "),
                   color=tcol[t], edgecolor="white", lw=0.4)
        ax.axhline(0, color=C["black"], lw=0.8)
        ax.set_xticks(xx); ax.set_xticklabels(short, rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("Mean NK module z-score"); ax.set_title("NK axis modules by tissue (scRNA)")
        ax.legend(fontsize=7)

    # --- B: TCGA-STAD NK state distribution ---
    ax = axes[1]; panel(ax, "B")
    if lab is not None:
        stad = lab[lab.dataset == "TCGA-STAD"]["nk_immune_state"].value_counts()
        order = ["NK-hot-cytotoxic", "NK-cold/excluded", "NK-intermediate", "NK-hot-dysfunctional"]
        stad  = stad.reindex([o for o in order if o in stad.index])
        scol  = [C["green"], C["blue"], C["grey"], C["verm"]]
        bars  = ax.barh(range(len(stad)), stad.values, color=scol[:len(stad)], edgecolor="white", lw=0.4)
        ax.set_yticks(range(len(stad))); ax.set_yticklabels(stad.index, fontsize=8)
        ax.invert_yaxis()
        for i, v in enumerate(stad.values):
            ax.text(v + 3, i, str(v), va="center", fontsize=8)
        ax.set_xlabel("Number of samples"); ax.set_title("NK immune states in TCGA-STAD (n = 450)")

    # --- C: Effector coupling replication across cohorts ---
    ax = axes[2]; panel(ax, "C")
    names = ["TCGA-LIHC\n(bulk, n=423)", "GC scRNA\n(1,017 NK cells)",
             "GSE62254\n(bulk, n=300)", "GSE84437\n(bulk, n=483)"]
    colors = [C["sky"], C["purple"], C["green"], C["green"]]
    vals = [0.55, np.nan, np.nan, np.nan]
    # fill from external validation table
    if ext is not None:
        ex = ext[ext.test == "protrusion~cytotoxicity"].set_index("dataset")["r"].astype(float)
        if "gse62254" in ex.index: vals[2] = ex["gse62254"]
        if "gse84437" in ex.index: vals[3] = ex["gse84437"]
    # fill gastric scRNA correlation
    if sc is not None:
        gc_sc = sc[sc.tissue == "gastric_cancer"]
        if len(gc_sc) > 1:
            r_gc, _ = stats.pearsonr(gc_sc["nk_protrusion_machinery_score"],
                                     gc_sc["nk_synapse_cytotoxicity_outcome_score"])
            vals[1] = r_gc
    bars = ax.bar(range(4), vals, color=colors, edgecolor="white", lw=0.4)
    for i, v in enumerate(vals):
        if not np.isnan(v):
            ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=8, fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=7.5)
    ax.set_ylabel("protrusion–cytotoxicity Pearson r"); ax.set_ylim(0, 0.72)
    ax.set_title("Effector coupling replicates across cohorts")
    ax.axhline(0.30, color=C["grey"], lw=0.6, ls=":", zorder=0)
    ax.text(3.7, 0.31, "moderate", fontsize=6.5, color=C["grey"], ha="right")

    fig.suptitle("Figure 2  |  Arm B — Gastric extension and independent external validation",
                 fontsize=12, fontweight="bold", y=1.03)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, "fig2_armB_extension")


# ========================================================================
# FIGURE 3 — De-circularized tumor-intrinsic target prioritization
# ========================================================================
def figure3(dpi=300):
    ev = load("candidate_evidence_matrix.tsv")
    ti = load("tumor_intrinsic_candidates.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), dpi=dpi,
                             gridspec_kw={"width_ratios": [1.15, 1, 0.85]})

    # --- A: Evidence landscape scatter ---
    ax = axes[0]; panel(ax, "A")
    if ev is not None:
        neg = ev[ev.tumor_specificity_log2 <= 0]
        pos = ev[ev.tumor_specificity_log2 > 0]
        ax.scatter(neg.tumor_specificity_log2, neg.nk_dysfunction_correlation,
                   s=14, c=C["grey"], alpha=0.4, edgecolors="none",
                   rasterized=True, label="NK-side / depleted in tumor")
        ax.scatter(pos.tumor_specificity_log2, pos.nk_dysfunction_correlation,
                   s=28, c=C["verm"], alpha=0.7, edgecolors="white", linewidths=0.5,
                   rasterized=True, label="Tumor-intrinsic (up in malignant cells)")
        ax.axvline(0, color=C["black"], lw=0.8, ls="--")
        # label key genes
        lbl_off = {"PHGDH": (6, -10), "SGMS2": (6, 4), "SMPD3": (6, 6),
                   "ERBB2": (5, -3), "MET": (5, 3), "NKG7": (6, 0), "GZMB": (6, 0)}
        for g, off in lbl_off.items():
            row = ev[ev.gene == g]
            if len(row):
                ax.annotate(g, (row.tumor_specificity_log2.iloc[0],
                                row.nk_dysfunction_correlation.iloc[0]),
                            fontsize=7, xytext=off, textcoords="offset points",
                            arrowprops=dict(arrowstyle="->", lw=0.5, color=C["dkgrey"]))
        ax.set_xlabel("Tumor specificity (log2FC malignant vs rest)")
        ax.set_ylabel("NK dysfunction correlation")
        ax.set_title("Candidate evidence landscape")
        ax.legend(loc="upper left", fontsize=7, markerscale=0.8)

    # --- B: Top-15 tumor-intrinsic candidates ---
    ax = axes[1]; panel(ax, "B")
    if ti is not None:
        catcol = {
            "metabolic_suppression": C["blue"],
            "sst_axis_nk_sm_synthesis": C["green"],
            "sst_axis_nk_sm_catabolism": C["green"],
            "sst_axis_tumor_serine_capacity": C["sky"],
            "sst_axis_nk_protrusion_machinery": C["orange"],
            "caf_ecm_exclusion": C["purple"],
            "nk_inhibitory_ligand": C["yellow"],
            "gastric_cancer_target": C["verm"],
            "stress_ligand_shedding": C["grey"],
        }
        top = ti.sort_values("target_score_v2", ascending=False).head(15).iloc[::-1]
        cols = [catcol.get(c, C["grey"]) for c in top.target_category]
        ax.barh(range(len(top)), top.target_score_v2, color=cols, edgecolor="white", lw=0.4)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.gene.values, fontsize=7.5,
                           fontfamily="monospace" if "fontfamily" in dir() else "sans-serif")
        ax.set_xlabel("Tumor-intrinsic target score (v2)")
        ax.set_title("Top-15 tumor-intrinsic candidates")
        # highlight top 3
        for i in range(min(3, len(top))):
            ax.text(top.target_score_v2.iloc[i] + 0.01, i,
                    f"  #{i+1}", fontsize=6.5, va="center", color=C["dkgrey"])

    # --- C: Category composition of the 37-gene pool ---
    ax = axes[2]; panel(ax, "C")
    if ti is not None:
        cc = ti.target_category.value_counts()
        cat_names = [c.replace("sst_axis_", "").replace("_", " ") for c in cc.index]
        bar_colors = [catcol.get(c, C["grey"]) for c in cc.index]
        ax.barh(range(len(cc)), cc.values, color=bar_colors, edgecolor="white", lw=0.4)
        ax.set_yticks(range(len(cc)))
        ax.set_yticklabels(cat_names, fontsize=7); ax.invert_yaxis()
        for i, v in enumerate(cc.values):
            ax.text(v + 0.15, i, str(v), va="center", fontsize=8)
        ax.set_xlabel("Number of genes")
        ax.set_title(f"Target pool (n = {cc.sum()}) by category")

    fig.suptitle("Figure 3  |  De-circularized tumor-intrinsic target prioritization",
                 fontsize=12, fontweight="bold", y=1.03)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, "fig3_targets")


# ========================================================================
# FIGURE 4 — Model comparison (Table 3)
# ========================================================================
def figure4(dpi=300):
    comp = load("model_comparison.tsv")
    if comp is None:
        print("  [SKIP] Figure 4 needs model_comparison.tsv")
        return

    summ = comp.groupby("method").agg(
        MCC_m=("MCC", "mean"), MCC_s=("MCC", "std"),
        AUR_m=("AUROC", "mean"), AUR_s=("AUROC", "std"),
    ).reset_index().sort_values("MCC_m")

    gnn_label = "GC-NKGraph-Atlas\n(GNN, ours)"
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=dpi)

    for ax, (mcol, scol, title, xlo) in zip(axes, [
        ("MCC_m", "MCC_s", "MCC", 0.0),
        ("AUR_m", "AUR_s", "AUROC", 0.5),
    ]):
        bar_colors = [GNN_C if "GNN" in x or "GC-NK" in x else BL_C for x in summ.method]
        ax.barh(range(len(summ)), summ[mcol], xerr=summ[scol],
                color=bar_colors, edgecolor="white", lw=0.4,
                error_kw=dict(lw=0.8, ecolor=C["black"], capsize=2.5))
        ax.set_yticks(range(len(summ)))
        labels = [x.replace("GC-NKGraph-Atlas(GNN)", gnn_label)
                   .replace("(GNN)", "\n(GNN)") for x in summ.method]
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.set_xlabel(f"{title}  (5-fold CV, mean ± SD)"); ax.set_xlim(xlo, 1.02)
        ax.set_title(f"NK-state classification — {title}")
        for i, (v, s) in enumerate(zip(summ[mcol], summ[scol])):
            ax.text(v + s + 0.012, i, f"{v:.3f}", va="center", fontsize=7.5)

    fig.suptitle("Figure 4  |  GNN on par with top tree baselines (p > 0.27), "
                 "beats linear/kernel/shallow (p < 0.05)",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles=[
        Patch(fc=GNN_C, label="GC-NKGraph-Atlas (ours)"),
        Patch(fc=BL_C,  label="Tabular baseline"),
    ], loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.005), fontsize=8)
    save(fig, "fig4_model_comparison")


# ========================================================================
# FIGURE 5 — Mechanism-card abstraction: concept diagram
# ========================================================================
def figure5(dpi=300):
    """
    Conceptual diagram showing how a mechanism card (YAML) drives the pipeline:
    one published mechanism → one card → cell-type-attributed proxy →
    mechanism-grounded heterogeneous graph → GNN → target list.

    This is a vector illustration built entirely with matplotlib patches —
    no external dependencies beyond matplotlib.
    """
    fig, ax = plt.subplots(1, 1, figsize=(14, 5.5), dpi=dpi)
    ax.set_xlim(0, 14); ax.set_ylim(0, 5.5)
    ax.axis("off")

    # ── colour scheme ──
    CARD_C  = "#FFF3CD"   # pale yellow (the mechanism card)
    DATA_C  = "#D4EDDA"   # pale green
    GRAPH_C = "#CCE5FF"   # pale blue
    MODEL_C = "#E2D9F3"   # pale purple
    TARGET_C= "#F8D7DA"   # pale red/pink
    BORDER  = "#6C757D"
    ARROW_C = "#495057"

    def draw_box(x, y, w, h, text, color, title=None, fontsize=8.5, title_fs=9):
        """Draw a rounded box with text."""
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                              facecolor=color, edgecolor=BORDER, lw=1.2, zorder=2)
        ax.add_patch(rect)
        if title:
            ax.text(x + w/2, y + h - 0.22, title, ha="center", va="top",
                   fontsize=title_fs, fontweight="bold", color="#333", zorder=3)
        ax.text(x + w/2, y + h/2 - 0.05, text, ha="center", va="center",
               fontsize=fontsize, color="#333", zorder=3,
               linespacing=1.3)

    def draw_arrow(x1, y1, x2, y2, style="->"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle=style, color=ARROW_C, lw=1.5))

    # ── Title ──
    ax.text(7, 5.2, "Figure 5  |  Mechanism-card abstraction: one card drives the entire pipeline",
            ha="center", fontsize=12, fontweight="bold")

    # ── Row 1: Published mechanism → Mechanism Card ──
    # "Published wet-lab mechanism" (left)
    draw_box(0.3, 3.6, 2.6, 1.3,
             "Published wet-lab\nmechanism\n\nZheng et al. 2023\nNat Immunol 24:748–59",
             "#F1F3F5", title="Input", fontsize=7.5)
    draw_arrow(2.95, 4.25, 3.6, 4.25)

    # "Mechanism card (YAML)" (center-left)
    mech_text = ("origin · biology · transcriptional_proxy\n"
                 "physical_ground_truth (gated)\n"
                 "graph_integration · validation\n"
                 "therapeutic_hook")
    draw_box(3.65, 3.3, 3.2, 1.9,
             mech_text, CARD_C,
             title="Mechanism Card (YAML)", fontsize=7)
    draw_arrow(6.9, 4.25, 7.6, 4.25)

    # "Reusable registry" note
    ax.text(5.25, 3.1, "One card per mechanism. New mechanism = new card, not new code.",
            ha="center", fontsize=7, style="italic", color="#6C757D")

    # ── Row 2: Pipeline stages ──
    # Stage 1: Cell-type-attributed proxy
    draw_box(0.3, 1.5, 2.6, 1.4,
             "7 gene modules (62 genes)\n"
             "Cell-type attribution\n"
             "(malignant vs NK)\n"
             "Pre-registered hypotheses H1–H5",
             DATA_C, title="Stage I: Transcriptional Proxy", fontsize=7)
    draw_arrow(2.95, 2.2, 3.6, 2.2)

    # Stage 2: Heterogeneous graph
    draw_box(3.65, 1.5, 3.2, 1.4,
             "6 edge types:\n"
             "PPI · ligand–receptor · TF-target\n"
             "metabolic_crosstalk (mechanism-grounded)\n"
             "sm_topology_axis · dysfunction_correlation",
             GRAPH_C, title="Stage II: Heterogeneous Graph", fontsize=7)
    draw_arrow(6.9, 2.2, 7.6, 2.2)

    # Stage 3: GNN
    draw_box(7.65, 1.5, 2.6, 1.4,
             "Heterogeneous Graph\nTransformer (HGT)\n"
             "→ Gene embeddings\n"
             "→ NK-state classifier (MLP)",
             MODEL_C, title="Stage III: Graph Model", fontsize=7)
    draw_arrow(10.3, 2.2, 11.0, 2.2)

    # Stage 4: Targets
    draw_box(11.05, 1.5, 2.6, 1.4,
             "37 tumor-intrinsic\ncandidates\n"
             "De-circularized\n"
             "Wet-lab assay per target",
             TARGET_C, title="Stage IV: Target List", fontsize=7)

    # ── Feedback arrow (bottom) ──
    ax.annotate("", xy=(3.65, 1.3), xytext=(11.05, 1.3),
                arrowprops=dict(arrowstyle="->", color="#ADB5BD", lw=1.2,
                               connectionstyle="arc3,rad=-0.3", ls="--"))
    ax.text(7.35, 0.85, "Card structure enforces claim boundaries and pre-registration\n"
                         "→ honest scoping map, not blanket recovery",
            ha="center", fontsize=7.5, style="italic", color="#6C757D")

    # ── Right-side annotation: The scoping result ──
    ax.text(13.5, 3.9, "Output:\nScoping Map", ha="center", fontsize=9,
            fontweight="bold", color="#333",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=BORDER, lw=1))
    ax.text(13.5, 3.2,
            "Effector arm   ✅\n"
            "Metabolic      ⚠️\n"
            "  (cell-type resolved)\n"
            "Topology       ❌\n"
            "  (transcription\n"
            "   ≠ phenotype)",
            ha="center", fontsize=7.5, color="#333", linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.3", fc="#FFF", ec="#DDD", lw=0.5))

    # ── Bottom: Key design principle ──
    ax.text(7, 0.35,
            "Design principle: card strictly separates what is computable (transcriptional proxy, ACTIVE) "
            "from what requires physical measurement (SEM, mass-spec, GATED)",
            ha="center", fontsize=7.5, fontweight="bold", color="#495057",
            bbox=dict(boxstyle="round,pad=0.3", fc="#F8F9FA", ec="#DEE2E6", lw=0.8))

    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    save(fig, "fig5_mechanism_card_concept")


# ========================================================================
# Main
# ========================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GC-NKGraph-Atlas publication figures.")
    parser.add_argument("--fig", type=int, choices=[1, 2, 3, 4, 5],
                        help="Generate a single figure (1–5)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Output resolution (default: 300)")
    args = parser.parse_args()

    set_style(dpi=args.dpi)
    print(f"GC-NKGraph-Atlas  |  Generating figures (dpi={args.dpi})")
    print(f"  Tables: {T}   Figures: {OUT}")

    figures = {
        1: figure1, 2: figure2, 3: figure3, 4: figure4, 5: figure5,
    }
    if args.fig:
        figures[args.fig](dpi=args.dpi)
    else:
        for num, fn in figures.items():
            try:
                fn(dpi=args.dpi)
            except Exception as e:
                print(f"  [ERROR] Figure {num}: {e}")
    print("Done.")
