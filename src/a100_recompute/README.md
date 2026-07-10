# A100 Recompute Scripts — P1 Science Tasks

Five self-contained Python scripts for T5/T6/T7/T14/T15. Each writes its output
to `results/tables/` and prints a PASS/FAIL at the end with specific verdicts for
manuscript writing.

## Quick start (A100)

```bash
# SSH to the A100 (via tailscale)
ssh user@100.112.165.109

# Navigate to the project
cd /data/lgh/GC-NKGraph-Atlas
conda activate gc-nkgraph

# Copy these scripts from local first:
#   scp -r src/a100_recompute user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/src/

# Run all five in order
bash src/a100_recompute/run_all.sh
```

## Task summary

| Script | Task | Time | Output |
|--------|------|------|--------|
| `run_h1_result.py` | T5 — H1 (serine capacity ~ sm balance) | <1 min | Appends to `sst_axis_positive_control_recovery.tsv` |
| `run_effect_size_reframe.py` | T15 — Effect sizes, CIs, BH correction | <1 min | `sst_axis_positive_control_recovery_v2.tsv` |
| `run_h3_activation_control.py` | T14 — Partial-correlation control for NK activation | <1 min | `h3_activation_control.tsv` |
| `run_decirc_audit.py` | T6 — De-circ verification audit | <1 min | `tumor_intrinsic_pool_audit.tsv` |
| `run_ablation.py` | T7 — Edge ablation (FULL vs −MC vs −SST) | 10-30 min | `ablation_results.tsv` |

## After A100 run completes

Copy the results back to local:

```powershell
scp -r user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/results/tables/sst_axis_positive_control_recovery.tsv results/tables/
scp -r user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/results/tables/sst_axis_positive_control_recovery_v2.tsv results/tables/
scp -r user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/results/tables/h3_activation_control.tsv results/tables/
scp -r user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/results/tables/tumor_intrinsic_pool_audit.tsv results/tables/
scp -r user@100.112.165.109:/data/lgh/GC-NKGraph-Atlas/results/tables/ablation_results.tsv results/tables/
```

Then regenerate figures and update manuscript:

```powershell
cd G:/cc/GC-NKGraph-Atlas

# Regenerate figures with updated numbers
python src/figures/make_figures.py

# Sync regenerated figures to manuscript + submission dirs
foreach ($d in @("manuscript/figures","submission_package/latex_official/figures","submission_package/03_figures")) {
  Copy-Item -Force results/figures/fig1_armA_positive_control.pdf $d/
  Copy-Item -Force results/figures/fig2_armB_extension.pdf $d/
  Copy-Item -Force results/figures/fig3_targets.pdf $d/
  Copy-Item -Force results/figures/fig4_model_comparison.pdf $d/
}

# Rebuild LaTeX from updated Markdown
cd submission_package/latex_official
python build_latex_submission.py
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## What each script's verdict means for the manuscript

### T5 (H1)
- **If H1 is null (r ≈ 0):** write "tumor_serine_capacity and nk_sm_balance are uncorrelated at the bulk level (H1: r=..., P=...), consistent with the metabolic coupling operating at a post-transcriptional level."
- **If H1 is significant:** write "[direction] correlation observed; however this may reflect [confounding]."

### T14 (H3 activation control)
- **If r_partial > 0.15:** the "effector arm recovers" claim is validated — protrusion~cytotoxicity is not just co-activation.
- **If r_partial 0.05–0.15:** the claim partially survives; write "partially independent of generic NK activation (partial r=...)."
- **If r_partial < 0.05:** DOWNGRADE to "protrusion-machinery and cytotoxicity transcripts co-vary with NK activation state."

### T7 (Ablation)
- **If ΔMCC(FULL − −MC) < 0.02:** the graph's value is in the embedding, not accuracy → reinforces the "comparable accuracy, added interpretability" framing already in the manuscript.
- **If ΔMCC > 0.02:** the metabolic_crosstalk edge genuinely contributes → strengthens the graph novelty argument.
- **If ΔMCC is negative:** the edge hurts performance — need to investigate.

### T15 (Effect sizes)
- **H2 single-cell: r ≈ 0.03, r² < 0.001** → the manuscript MUST frame this as "statistically detectable but of negligible magnitude."
- **BH correction:** any test that loses significance at FDR<0.05 must be explicitly noted.

### T6 (De-circ audit)
- **If any NK-effector genes (RAC1, NKG7, etc.) are in the n=37 tumor pool:** the de-circularization script needs to be patched to exclude them, then re-run split_target_lists.py.
- **If NK-side modules dominate the category count:** the category labeling in Fig 3C needs a clarifying note in the figure legend.
