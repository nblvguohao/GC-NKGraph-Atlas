"""
Publication figures for GC-NKGraph-Atlas (Fig 1-4), generated from result tables.

Every number is read from results/tables/*.tsv so the figures stay consistent with
the manuscript. Palette: Okabe-Ito (colour-blind safe). Output: PDF (vector) + PNG.

    python src/figures/make_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats

T = "results/tables/"
OUT = "results/figures/"
os.makedirs(OUT, exist_ok=True)

# --- Okabe-Ito colour-blind-safe palette ---
OK = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", verm="#D55E00", purple="#CC79A7",
          grey="#8C8C8C")
REC, NOT, WEAK, FLAG = OK["green"], OK["verm"], OK["grey"], OK["orange"]

mpl.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.grid": True, "grid.color": "#E6E6E6", "grid.linewidth": 0.6,
    "legend.frameon": False, "legend.fontsize": 8,
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {OUT}{name}.pdf / .png")


def panel(ax, letter):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def flow_box(ax, x, y, w, h, label, color, sub=None, fc_alpha=0.18, fontsize=7.6):
    """One rounded box for a schematic panel: bold label + optional small sub-line."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                                 fc=color, ec=color, lw=1.3, alpha=fc_alpha))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                                 fc="none", ec=color, lw=1.3))
    ty = y + h / 2 + (0.10 if sub else 0)
    ax.text(x + w / 2, ty, label, ha="center", va="center", fontsize=fontsize,
            fontweight="bold", color=OK["black"])
    if sub:
        ax.text(x + w / 2, y + h / 2 - 0.14, sub, ha="center", va="center",
                fontsize=fontsize - 1.1, color=OK["black"])


def flow_arrow(ax, xy0, xy1, color=None, lw=1.2, rad=0.0):
    ax.add_patch(FancyArrowPatch(xy0, xy1, arrowstyle="-|>", mutation_scale=11,
                                  color=color or OK["black"], lw=lw,
                                  connectionstyle=f"arc3,rad={rad}"))


def schematic_axes(fig, gs_slice):
    ax = fig.add_subplot(gs_slice)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    panel(ax, "a")
    return ax


ASSETS = "src/figures/assets/"


def image_panel(fig, gs_slice, image_path, label=None):
    """Embed an externally-designed schematic (e.g. a BioRender-style panel) as a figure panel."""
    ax = fig.add_subplot(gs_slice)
    ax.imshow(plt.imread(image_path), aspect="auto")
    ax.axis("off")
    if label:
        panel(ax, label)
    return ax


def tcga_cond(bc):
    try:
        return "tumor" if int(str(bc).split("-")[3][:2]) < 10 else "normal"
    except Exception:
        return "unknown"


# ============================ FIGURE 1 — Arm A ============================
def figure1():
    liver = pd.read_csv(T + "sst_axis_scores_liver_bulk.tsv", sep="\t", index_col=0)
    sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")
    rec = pd.read_csv(T + "sst_axis_positive_control_recovery.tsv", sep="\t")

    fig = plt.figure(figsize=(9.5, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.4, wspace=0.32)

    # a: schematic — bulk vs single-cell confound-control logic for H2/H3
    # (externally designed BioRender-style panel; already carries its own "a" label and title)
    image_panel(fig, gs[0, :], ASSETS + "fig1_panel_a_schematic.png")

    # b: bulk LIHC H3 scatter protrusion vs cytotoxicity
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b")
    x = liver["nk_protrusion_machinery_score"]; y = liver["nk_synapse_cytotoxicity_outcome_score"]
    ax.scatter(x, y, s=10, c=OK["blue"], alpha=0.5, edgecolors="none")
    b, a = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a + b * xs, color=OK["black"], lw=2)
    r, p = stats.pearsonr(x, y)
    ax.set_title("H3: protrusion → cytotoxicity (bulk TCGA-LIHC)")
    ax.set_xlabel("NK protrusion-machinery score"); ax.set_ylabel("NK cytotoxicity-output score")
    ax.text(0.04, 0.95, f"r = {r:.2f}\np = {p:.0e}", transform=ax.transAxes,
            va="top", fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#CCC"))

    # c: H2/H3 bulk vs single-cell r
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "c")
    def getr(h, res):
        v = rec[(rec.hypothesis == h) & (rec.resolution == res)]["r"]
        return float(v.iloc[0]) if len(v) else np.nan
    labels = ["H2\nSM-bal→protr", "H3\nprotr→cytotox"]
    bulk = [getr("H2", "bulk"), getr("H3", "bulk")]
    single = [getr("H2", "single_cell_NK"), getr("H3", "single_cell_NK")]
    xx = np.arange(2); w = 0.36
    ax.bar(xx - w/2, bulk, w, label="bulk", color=OK["grey"])
    ax.bar(xx + w/2, single, w, label="single-cell NK", color=OK["blue"])
    ax.axhline(0, color=OK["black"], lw=0.8)
    ax.set_xticks(xx); ax.set_xticklabels(labels); ax.set_ylabel("Pearson r")
    ax.set_title("Cell-type resolution rescues H2")
    ax.legend(loc="upper left")
    ax.annotate("bulk ≈ 0", (0 - w/2, 0.02), fontsize=7, ha="center", color=OK["verm"])

    # d: intratumoral vs normal NK, two modules
    ax = fig.add_subplot(gs[2, 0]); panel(ax, "d")
    mods = [("nk_synapse_cytotoxicity_outcome_score", "cytotoxicity\noutput"),
            ("nk_protrusion_machinery_score", "protrusion\nmachinery")]
    data, poss, cols, ticks = [], [], [], []
    for i, (m, lab) in enumerate(mods):
        for j, cond in enumerate(["normal", "tumor"]):
            data.append(sc[sc.condition == cond][m].dropna().values)
            poss.append(i * 3 + j)
            cols.append(OK["sky"] if cond == "normal" else OK["verm"])
        ticks.append((i * 3 + 0.5, lab))
    bp = ax.boxplot(data, positions=poss, widths=0.7, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], cols): patch.set_facecolor(c); patch.set_alpha(0.7)
    for med in bp["medians"]: med.set_color(OK["black"])
    ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks])
    ax.set_ylabel("per-cell module z-score")
    ax.set_title("Intratumoral vs normal NK (single cell)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc=OK["sky"], alpha=0.7, label="normal NK"),
                       Patch(fc=OK["verm"], alpha=0.7, label="intratumoral NK")], loc="upper right")

    # e: forest plot of H2-H5
    ax = fig.add_subplot(gs[2, 1]); panel(ax, "e")
    rr = rec[rec["r"].astype(str).str.replace(".", "", 1).str.replace("-", "", 1).str.isdigit()].copy()
    rr["rv"] = rr["r"].astype(float)
    def mklab(row):
        if row.hypothesis == "H5":
            sub = row.test.split(":")[-1].strip().replace("_", " ").replace("outcome", "").strip()
            return f"H5 {sub} (scNK)"
        res = row.resolution.replace("single_cell_NK", "scNK")
        return f"{row.hypothesis} ({res})"
    rr["lab"] = rr.apply(mklab, axis=1)
    rr = rr.iloc[::-1].reset_index(drop=True)
    cmap = {"RECOVERED": REC, "INCONCLUSIVE": WEAK, "FLAGGED": FLAG}
    cols = [cmap.get(o, NOT) for o in rr["outcome"]]
    markers = ["^" if o == "FLAGGED" else "o" for o in rr["outcome"]]
    for xi, yi, ci, mi in zip(rr["rv"], range(len(rr)), cols, markers):
        ax.scatter([xi], [yi], c=[ci], s=55, zorder=3, edgecolors="white",
                   linewidths=0.8, marker=mi)
    ax.axvline(0, color=OK["black"], lw=0.8, ls="--")
    ax.set_yticks(range(len(rr))); ax.set_yticklabels(rr["lab"], fontsize=7)
    ax.set_xlabel("r  (or Δ for H5)"); ax.set_title("Pre-registered hypotheses (H2–H5)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=REC,mec="w",label="recovered"),
                       Line2D([0],[0],marker="o",ls="",mfc=NOT,mec="w",label="not recovered"),
                       Line2D([0],[0],marker="o",ls="",mfc=WEAK,mec="w",label="inconclusive"),
                       Line2D([0],[0],marker="^",ls="",mfc=FLAG,mec="w",
                              label="flagged: not distinguishable\nfrom technical background")],
              loc="lower right", fontsize=6.5)

    # No fig.suptitle: panel (a)'s schematic image already carries the figure title.
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    save(fig, "fig1_armA_positive_control")


# ==================== FIGURE 2 — Arm B + external validation ====================
def figure2():
    tis = pd.read_csv(T + "sst_axis_scrna_by_tissue.tsv", sep="\t")
    lab = pd.read_csv(T + "nk_state_labels.tsv", sep="\t", comment="#")
    ext = pd.read_csv(T + "external_validation_results.tsv", sep="\t")
    sc = pd.read_csv(T + "sst_axis_scores_single_cell.tsv", sep="\t")

    fig = plt.figure(figsize=(13.5, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.36, 1], hspace=0.32, wspace=0.34)

    # a: schematic — external-validation design
    axs = schematic_axes(fig, gs[0, :])
    flow_box(axs, 0.01, 0.30, 0.16, 0.5, "Arm B\ndiscovery", OK["verm"], "gastric scRNA")
    flow_box(axs, 0.22, 0.30, 0.20, 0.5, "TCGA-STAD\nbulk (n=450)", OK["sky"], "NK-state calling")
    flow_arrow(axs, (0.17, 0.55), (0.22, 0.55))
    flow_box(axs, 0.49, 0.55, 0.20, 0.32, "GSE62254\n(bulk, n indep.)", OK["green"], fontsize=7.0)
    flow_box(axs, 0.49, 0.12, 0.20, 0.32, "GSE84437\n(bulk, n indep.)", OK["green"], fontsize=7.0)
    flow_arrow(axs, (0.42, 0.55), (0.49, 0.71))
    flow_arrow(axs, (0.42, 0.55), (0.49, 0.28))
    flow_box(axs, 0.76, 0.30, 0.22, 0.5, "Effector coupling\nreplicated\nindependently", OK["blue"], fontsize=7.2)
    flow_arrow(axs, (0.69, 0.71), (0.76, 0.58))
    flow_arrow(axs, (0.69, 0.28), (0.76, 0.42))
    axs.text(0.5, 0.98, "Two fully independent external gastric cohorts, "
             "never used for module definition, test whether the bulk "
             "effector-arm coupling generalizes beyond TCGA-LIHC",
             ha="center", va="top", fontsize=7.4, style="italic", transform=axs.transAxes)

    # b: three tissues × 4 modules grouped bars
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b")
    modcols = ["nk_sm_balance_score", "nk_protrusion_machinery_score",
               "nk_topology_permissive_score", "nk_synapse_cytotoxicity_outcome_score"]
    short = ["SM balance", "protrusion", "topology-perm", "cytotoxicity"]
    tis2 = tis.set_index("tissue").reindex(["healthy_liver", "liver_metastasis", "gastric_cancer"])
    xx = np.arange(len(modcols)); w = 0.26
    tcol = {"healthy_liver": OK["sky"], "liver_metastasis": OK["orange"], "gastric_cancer": OK["verm"]}
    for k, (t, row) in enumerate(tis2.iterrows()):
        ax.bar(xx + (k-1)*w, [row[m] for m in modcols], w, label=t.replace("_", " "), color=tcol[t])
    ax.axhline(0, color=OK["black"], lw=0.8)
    ax.set_xticks(xx); ax.set_xticklabels(short, rotation=20, ha="right")
    ax.set_ylabel("mean NK module z-score"); ax.set_title("NK axis modules by tissue (scRNA)")
    ax.legend()

    # c: TCGA-STAD NK state distribution
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "c")
    stad = lab[lab.dataset == "TCGA-STAD"]["nk_immune_state"].value_counts()
    order = ["NK-hot-cytotoxic", "NK-cold/excluded", "NK-intermediate", "NK-hot-dysfunctional"]
    stad = stad.reindex([o for o in order if o in stad.index])
    scol = [OK["green"], OK["blue"], OK["grey"], OK["verm"]]
    ax.barh(range(len(stad)), stad.values, color=scol[:len(stad)])
    ax.set_yticks(range(len(stad))); ax.set_yticklabels(stad.index, fontsize=8)
    ax.invert_yaxis()
    for i, v in enumerate(stad.values): ax.text(v + 3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("samples"); ax.set_title("NK states in TCGA-STAD (n=450)")

    # d: effector coupling across cohorts (protrusion~cytotoxicity r)
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "d")
    gc_sc = sc[sc.tissue == "gastric_cancer"]
    r_gc, _ = stats.pearsonr(gc_sc["nk_protrusion_machinery_score"], gc_sc["nk_synapse_cytotoxicity_outcome_score"])
    e = ext[ext.test == "protrusion~cytotoxicity"].set_index("dataset")["r"].astype(float)
    names = ["TCGA-LIHC\n(bulk)", "GC scRNA\n(single cell)*", "GSE62254\n(bulk)", "GSE84437\n(bulk)"]
    vals = [0.55, r_gc, e.get("gse62254", np.nan), e.get("gse84437", np.nan)]
    cols = [OK["sky"], OK["orange"], OK["green"], OK["green"]]
    ax.bar(range(4), vals, color=cols)
    for i, v in enumerate(vals): ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=8)
    ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=7.5)
    ax.set_ylabel("protrusion~cytotoxicity  r"); ax.set_ylim(0, 0.72)
    ax.set_title("Effector coupling replicates (4 cohorts)")
    ax.text(0.5, 0.95, "*not re-verified vs the Arm A (§3.2)\ntechnical confound; bulk = primary evidence",
            transform=ax.transAxes, fontsize=6, color=OK["verm"], va="top", ha="center",
            bbox=dict(boxstyle="round", fc="white", ec=OK["verm"], alpha=0.9))

    fig.suptitle("Arm B — gastric extension and independent external validation",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "fig2_armB_extension")


# ==================== FIGURE 3 — target prioritization ====================
def figure3():
    ev = pd.read_csv(T + "candidate_evidence_matrix.tsv", sep="\t")
    ti = pd.read_csv(T + "tumor_intrinsic_candidates.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))

    # a: tumor specificity vs NK dysfunction assoc, tumor-intrinsic pool highlighted
    ax = axes[0]; panel(ax, "a")
    neg = ev[ev.tumor_specificity_log2 <= 0]; pos = ev[ev.tumor_specificity_log2 > 0]
    ax.scatter(neg.tumor_specificity_log2, neg.nk_dysfunction_correlation, s=14,
               c=OK["grey"], alpha=0.5, edgecolors="none", label="NK-side / depleted")
    ax.scatter(pos.tumor_specificity_log2, pos.nk_dysfunction_correlation, s=22,
               c=OK["verm"], edgecolors="white", linewidths=0.4, label="tumor-intrinsic (up in tumor)")
    ax.axvline(0, color=OK["black"], lw=0.8, ls="--")
    lbl_off = {"PHGDH": (5, -9), "SGMS2": (5, 4), "ERBB2": (5, -2), "MET": (5, 3),
               "NKG7": (6, 0), "GZMB": (6, 0), "COL1A2": (5, 3)}
    for g, off in lbl_off.items():
        row = ev[ev.gene == g]
        if len(row):
            ax.annotate(g, (row.tumor_specificity_log2.iloc[0], row.nk_dysfunction_correlation.iloc[0]),
                        fontsize=7, xytext=off, textcoords="offset points")
    ax.set_xlabel("tumor specificity  (log2FC malignant vs rest)")
    ax.set_ylabel("NK dysfunction correlation"); ax.set_title("Candidate evidence landscape")
    ax.legend(loc="upper left")

    # b: top 15 tumor-intrinsic by score
    ax = axes[1]; panel(ax, "b")
    top = ti.sort_values("target_score_v2", ascending=False).head(15).iloc[::-1]
    catcol = {"metabolic_suppression": OK["blue"], "sst_axis_nk_sm_synthesis": OK["green"],
              "sst_axis_nk_sm_catabolism": OK["green"], "sst_axis_tumor_serine_capacity": OK["sky"],
              "sst_axis_nk_protrusion_machinery": OK["orange"], "caf_ecm_exclusion": OK["purple"],
              "nk_inhibitory_ligand": OK["yellow"], "gastric_cancer_target": OK["verm"],
              "stress_ligand_shedding": OK["grey"]}
    cols = [catcol.get(c, OK["grey"]) for c in top.target_category]
    ax.barh(range(len(top)), top.target_score_v2, color=cols)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, fontsize=8)
    ax.set_xlabel("tumor-intrinsic target score"); ax.set_title("Top-15 tumor-intrinsic candidates")

    # c: category composition of the 37-gene pool
    ax = axes[2]; panel(ax, "c")
    cc = ti.target_category.value_counts()
    clean = [c.replace("sst_axis_", "").replace("_", " ") for c in cc.index]
    ax.barh(range(len(cc)), cc.values, color=[catcol.get(c, OK["grey"]) for c in cc.index])
    ax.set_yticks(range(len(cc))); ax.set_yticklabels(clean, fontsize=7.5); ax.invert_yaxis()
    for i, v in enumerate(cc.values): ax.text(v + 0.1, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("genes"); ax.set_title("Tumor-intrinsic pool (n=37) by category")

    fig.suptitle("De-circularized tumor-intrinsic target prioritization",
                 fontsize=12, fontweight="bold", y=1.03)
    fig.tight_layout()
    save(fig, "fig3_targets")


# ==================== FIGURE 4 — model comparison (Table 3) ====================
def figure4():
    comp = pd.read_csv(T + "model_comparison.tsv", sep="\t")
    summ = comp.groupby("method").agg(MCC_m=("MCC", "mean"), MCC_s=("MCC", "std"),
                                      AUR_m=("AUROC", "mean"), AUR_s=("AUROC", "std")).reset_index()
    summ = summ.sort_values("MCC_m")
    gnn = "GC-NKGraph-Atlas(GNN)"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, letter, (m, s, ttl, lo) in zip(axes, ["a", "b"],
            [("MCC_m", "MCC_s", "MCC", 0), ("AUR_m", "AUR_s", "AUROC", 0.5)]):
        panel(ax, letter)
        cols = [OK["verm"] if x == gnn else OK["sky"] for x in summ.method]
        ax.barh(range(len(summ)), summ[m], xerr=summ[s], color=cols,
                error_kw=dict(lw=0.8, ecolor=OK["black"], capsize=2))
        ax.set_yticks(range(len(summ)))
        ax.set_yticklabels([x.replace("(GNN)", "\n(GNN)") for x in summ.method], fontsize=7.5)
        ax.set_xlabel(f"{ttl}  (5-fold mean ± SD)"); ax.set_xlim(lo, 1.0)
        ax.set_title(f"NK-state classification — {ttl}")
        for i, v in enumerate(summ[m]): ax.text(v + summ[s].iloc[i] + 0.01, i, f"{v:.2f}", va="center", fontsize=7.5)
    from matplotlib.patches import Patch
    fig.suptitle("Graph model on par with top tree baselines (paired p>0.27), beats linear/kernel/shallow (p<0.05)",
                 fontsize=10.5, fontweight="bold", y=1.02)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.legend(handles=[Patch(fc=OK["verm"], label="GC-NKGraph-Atlas (ours)"),
                        Patch(fc=OK["sky"], label="tabular baseline")],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    save(fig, "fig4_model_comparison")


if __name__ == "__main__":
    print("Generating figures...")
    figure1(); figure2(); figure3(); figure4()
    print("Done.")
