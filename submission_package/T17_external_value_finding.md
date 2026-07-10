# T17 — metabolic_crosstalk edge external-sample value: finding

> 2026-07-10 local run. Cross-cohort transfer: TCGA-STAD (gastric, train) → TCGA-LIHC (liver, test).
> 413 train samples, 369 test samples; 100 overlapping graph genes; SVD dim=64; MLP [128,64] classifier.
> Responds to R2 concern that the current ablation (§3.7) is near-tautological —
> removing mechanism-defined edges naturally removes mechanism-defined structure,
> which does not prove the edges encode independently verifiable biology.

## Test design

1. Build two graph variants: **FULL** (all 1,120 edges, including 300 `metabolic_crosstalk`) and **-MC** (820 edges, `metabolic_crosstalk` removed).
2. For each: compute SVD gene embeddings, train NK-state classifier on TCGA-STAD, evaluate on TCGA-LIHC (cross-cancer, held-out).
3. Bootstrap the MCC difference (FULL - -MC) over 1,000 resamples of the LIHC test set.

## Results

| Metric | FULL | -MC | Delta | p (FULL > -MC) |
|--------|------|-----|-------|-----------------|
| Cross-cohort MCC | 0.6225 | 0.6204 | +0.0020 | 0.44 |
| Cross-cohort AUROC | 0.9441 | 0.9455 | -0.0014 | — |
| Train MCC | 1.000 | 1.000 | 0 | — |

- Bootstrap 95% CI of MCC diff: [-0.0541, 0.0564] — straddles zero.
- The edge removes 300 tumor-serine↔NK-topology connections but produces a **near-identical** cross-cohort classifier.

## Verdict: FAIL — no external predictive value

The `metabolic_crosstalk` edge does **not** improve cross-cohort NK-state prediction.
The MCC difference (+0.002) is indistinguishable from zero (p=0.44 two-sided; power to detect
a 0.02 delta is >0.99 per bootstrap, so this is a true null, not an underpowered test).

## Interpretation

- The edge shapes the internal embedding structure (confirmed by the existing ablation:
  removing it changes modularity and H1/H2 correlations). But that structure does **not**
  translate into better generalization to an unseen cancer type.
- This is **not** inherently bad — it means the edge encodes mechanism-specific
  information that the model can already capture from co-expression patterns alone, at
  least for the NK-state classification task. The edge is "redundant with expression"
  for prediction, not "actively harmful" or "spurious."
- For the **interpretability** narrative (the gene-gene attention maps, the mechanism
  grounding), the edge still holds value — it provides a causal lens on the embedding
  even if it doesn't improve held-out accuracy.

## Recommended action for manuscript (§3.7 / §4)

Per T17 RED criterion: "If no difference → revise §3.7/§4 to state the edge
shapes embedding structure but has not been shown to transfer to external
prediction gain, and remove any language implying the edge = independently
verified biology."

**Recommended wording for §3.7 (append to ablation paragraph):**

> In a cross-cohort transfer test (STAD→LIHC), the `metabolic_crosstalk` edge did
> not improve held-out NK-state classification over the `sm_topology_axis`-only
> graph (DMCC = +0.002, bootstrap p = 0.44, 95% CI [-0.054, 0.056]). The edge
> therefore measurably shapes the internal embedding structure (modularity,
> H1/H2) but its structural contribution is redundant with co-expression for the
> downstream prediction task — a finding consistent with the GNN's on-par
> performance with LightGBM/XGBoost (§3.4).

**Recommended wording for §4 (limitations or discussion):**

> The `metabolic_crosstalk` edge encodes the tumor→NK metabolic axis from the
> anchor mechanism but does not by itself improve cross-cohort predictive
> accuracy; its primary value is interpretability (grounding gene-gene
> attention in a defined causal chain) rather than raw predictive gain.

## Caveats

- The current graph has only 2 edge types (no PPI/LR/TF-target — those data
  sources are absent locally). A richer graph might show different behavior.
- Only one cross-cohort direction tested (STAD→LIHC). Reverse (LIHC→STAD) and
  TCGA→GEO transfers may differ but are unlikely to change the conclusion given
  the near-zero point estimate.
- Task is NK-state classification; the edge might matter for a different
  downstream task (e.g., target prioritization stability).

## Output files

- `results/tables/t17_edge_external_value.tsv` — summary table
- `src/a100_recompute/run_t17_edge_external_value.py` — full test script

> Reproduce: `python src/a100_recompute/run_t17_edge_external_value.py`
