#!/usr/bin/env python3
"""
Sync Markdown source (main_manuscript.md) with LaTeX fixes applied to main.tex.

T1–T4 changes need to be mirrored in the Markdown source so the two stay
consistent. Run idempotently; safe to re-run.
"""
from __future__ import annotations
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MD = ROOT / "manuscript" / "main_manuscript.md"

txt = MD.read_text(encoding="utf-8")

# --- T2: remove hard-coded "Table N." prefix from Results tables ---
# Table 1 (datasets)
txt = txt.replace(
    "**Table 1. Dataset characteristics.**",
    "**Dataset characteristics.**"
)
# Table 2 (hypothesis outcomes)
txt = txt.replace(
    "**Table 2. Pre-registered hypothesis outcomes (multi-resolution).**",
    "**Pre-registered hypothesis outcomes (multi-resolution).**"
)
# Table 3 (classification)
txt = txt.replace(
    "**Table 3. NK-state classification (TCGA-STAD, 5-fold CV; mean over folds).**",
    "**NK-state classification (TCGA-STAD, 5-fold CV; mean over folds).**"
)
# Table 4 (targets)
txt = txt.replace(
    "**Table 4. Top putative tumor-intrinsic candidate targets in gastric cancer (excerpt).**",
    "**Top putative tumor-intrinsic candidate targets in gastric cancer (excerpt).**"
)
# Table 5 (external validation)
txt = txt.replace(
    "**Table 5. External validation of the axis (independent gastric cohorts).**",
    "**External validation of the axis (independent gastric cohorts).**"
)

# --- T3: "Figure 0" -> "Figure 1" (workflow figure) ---
# The Markdown body references "Figure 0" for the workflow; rename to just
# describe it without a baked-in number (LaTeX provides the number).
txt = txt.replace(
    "**Figure 0.** Workflow overview of the GC-NKGraph-Atlas framework showing the five stages from data acquisition through target prioritization, the mechanism-card abstraction, and the two-arm study design (Arm A: liver positive control; Arm B: gastric cancer extension).",
    "**Workflow overview.** Workflow of the GC-NKGraph-Atlas framework, from data acquisition through target prioritization, showing the mechanism-card abstraction and the two-arm study design (Arm A: liver positive control; Arm B: gastric cancer extension)."
)

# --- T3: add model-comparison figure after the GNN classification text ---
# Insert after the `*GNN not significantly different...` line
model_fig_block = """*GNN not significantly different from LightGBM/XGBoost (paired p>0.27); significantly above ElasticNet/SVM/MLP (paired p<0.05).*

**Model comparison.** Bar chart of NK-state classification performance (TCGA-STAD, 5-fold stratified CV): MCC and AUROC for GC-NKGraph-Atlas (GNN) versus six tabular baselines on identical seed-42 folds; error bars show cross-fold standard deviation. The GNN is statistically on par with LightGBM/XGBoost and significantly above ElasticNet/SVM/MLP."""

old_block = """*GNN not significantly different from LightGBM/XGBoost (paired p>0.27); significantly above ElasticNet/SVM/MLP (paired p<0.05).*

### 3.5 Candidate target prioritization (de-circularized)"""

txt = txt.replace(old_block, model_fig_block + """

### 3.5 Candidate target prioritization (de-circularized)""")

# --- T4: remove "five stages" from figure description ---
txt = txt.replace(
    "five stages from data acquisition through target prioritization",
    "pipeline stages from data acquisition through target prioritization"
)

# --- T1 proxy fix: the Markdown source is clean but the LaTeX converter had bugs ---
# The Markdown source already writes cell-cell correctly; no fix needed here.
# The `metabolic_crosstalk` edge name is correct in the Markdown already.

with MD.open("w", encoding="utf-8", newline="\n") as f:
    f.write(txt)
print(f"Updated {MD}")
