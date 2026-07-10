"""
Figure 0 (workflow overview) for GC-NKGraph-Atlas.

Standalone from make_figures.py because this is a pipeline schematic, not a
data figure -- it does not read results/tables/*.tsv. Stage names, phase
numbers, and content strictly match Methods Table 1 (Study design overview,
main_manuscript.md / main.tex Section 2.1): five stages, DATA / scRNA / GRAPH
/ MODEL / TARGETS. The metabolic_crosstalk edge is labelled "Zheng 2023" to
match the Heterogeneous graph edge-types table (Section 2.5), not "Reactome".

    python src/figures/make_workflow_figure.py
"""
import os
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "results/figures/"
os.makedirs(OUT, exist_ok=True)

# --- Okabe-Ito colour-blind-safe palette (matches make_figures.py) ---
OK = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", verm="#D55E00", purple="#CC79A7",
          grey="#8C8C8C")

mpl.rcParams.update({
    "font.size": 9, "font.family": "sans-serif",
})

STAGES = [
    ("STAGE I\nDATA", "Phases 1–2", OK["purple"],
     "TCGA-STAD, TCGA-LIHC,\nGEO gastric cohorts\n+ mechanism-card prior"),
    ("STAGE II\nscRNA", "Phases 3–7", OK["orange"],
     "scRNA integration, NK atlas\nannotation, SST-axis proxy\n(7 modules), trajectory"),
    ("STAGE III\nGRAPH", "Phase 8", OK["blue"],
     "Heterogeneous gene graph:\nPPI / LR / TF / metabolic_\ncrosstalk (Zheng 2023) / SST"),
    ("STAGE IV\nMODEL", "Phases 9–10", OK["green"],
     "Baseline comparison +\nGNN gene embedding →\nNK-state classifier"),
    ("STAGE V\nTARGETS", "Phases 11–14R", OK["verm"],
     "SST-axis scoring, candidate\nprioritization, assay\nrecommendation"),
]

fig, ax = plt.subplots(figsize=(14, 6.2))
ax.set_xlim(0, 14)
ax.set_ylim(0, 6.2)
ax.axis("off")

fig.suptitle("GC-NKGraph-Atlas workflow overview", fontsize=15,
             fontweight="bold", y=0.99)

# --- mechanism-card prior feeding into Stage I ---
card = FancyBboxPatch((0.3, 5.15), 2.0, 0.75, boxstyle="round,pad=0.06",
                       fc="white", ec=OK["black"], lw=1.2)
ax.add_patch(card)
ax.text(1.3, 5.52, "mechanism-card\n(machine-readable prior)", ha="center",
        va="center", fontsize=8.2, fontweight="bold")
ax.add_patch(FancyArrowPatch((1.3, 5.15), (1.3, 4.55), arrowstyle="-|>",
                              mutation_scale=14, color=OK["black"], lw=1.3))

# --- five stage boxes ---
box_w, box_h, gap = 2.15, 1.85, 0.32
x0, y0 = 0.3, 2.5
centers = []
for i, (title, phases, color, content) in enumerate(STAGES):
    x = x0 + i * (box_w + gap)
    centers.append(x + box_w / 2)
    box = FancyBboxPatch((x, y0), box_w, box_h, boxstyle="round,pad=0.05",
                          fc=color, ec="none", alpha=0.16)
    ax.add_patch(box)
    header = FancyBboxPatch((x, y0 + box_h - 0.55), box_w, 0.55,
                             boxstyle="round,pad=0.02", fc=color, ec="none")
    ax.add_patch(header)
    ax.text(x + box_w / 2, y0 + box_h - 0.275, title, ha="center", va="center",
            fontsize=9.3, fontweight="bold", color="white")
    ax.text(x + box_w / 2, y0 + box_h - 0.68, phases, ha="center", va="top",
            fontsize=7.3, style="italic", color=OK["black"])
    ax.text(x + box_w / 2, y0 + 0.62, content, ha="center", va="center",
            fontsize=6.9, color=OK["black"])
    if i > 0:
        xa = x0 + (i - 1) * (box_w + gap) + box_w
        ax.add_patch(FancyArrowPatch((xa, y0 + box_h / 2), (x, y0 + box_h / 2),
                                      arrowstyle="-|>", mutation_scale=14,
                                      color=OK["black"], lw=1.3))

# --- two-arm validation band below Stage V ---
arm_y = 0.35
for label, sub, color, x in [
    ("Arm A — Liver", "positive control", OK["sky"], centers[-1] - 1.15),
    ("Arm B — Gastric", "extension", OK["green"], centers[-1] + 1.15),
]:
    box = FancyBboxPatch((x - 1.0, arm_y), 2.0, 0.85, boxstyle="round,pad=0.05",
                          fc="white", ec=color, lw=1.6)
    ax.add_patch(box)
    ax.text(x, arm_y + 0.56, label, ha="center", va="center", fontsize=8.6,
            fontweight="bold", color=color)
    ax.text(x, arm_y + 0.22, sub, ha="center", va="center", fontsize=7.4,
            color=OK["black"])
    ax.add_patch(FancyArrowPatch((centers[-1], y0), (x, arm_y + 0.85),
                                  arrowstyle="-|>", mutation_scale=13,
                                  color=color, lw=1.2,
                                  connectionstyle="arc3,rad=0.0"))

plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}fig0_workflow.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT}fig0_workflow.pdf / .png")
