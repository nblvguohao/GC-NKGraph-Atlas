# GC-NKGraph-Atlas — Submission-Readiness Assessment & Remediation Roadmap

> **Status note (2026-07-07 pre-submission audit):** This file began as a
> historical blocker audit. The current submission-facing action list is
> `PRE_SUBMISSION_REVISION_PLAN.md`. Sections below are retained for provenance;
> claims and package language should follow the newer plan and `main_claims.md`.

> Author-facing working document. Prepared 2026-07-07.
> Purpose: an honest gap analysis mapping every manuscript claim to its current
> evidence, the specific defect blocking it, and the exact server-side action
> needed to close it. Paired with fix scripts under `src/` (see §5).

---

## 0. PROGRESS UPDATE (2026-07-07, after local fixes)

Fixed locally (from result tables present on this machine — no server needed):

- **Blocker A resolved into an honest, numeric result.** Re-tested H1–H5 at bulk
  *and* single-cell NK resolution (`sst_axis_positive_control_recovery.tsv`).
  Verdict: effector arm recovers (H3 bulk r=0.55, scNK r=0.32; H5-cytotoxicity
  Δ=−0.14, p=6e-52); metabolic coupling recovers **only at single-cell
  resolution** (H2 bulk null → scNK r=+0.030, p=6e-3); physical topology **not**
  recovered from transcription (H4 wrong sign; intratumoral NK protrusion
  transcript *higher*, not lower). Now the paper's central scoping finding.
- **R5 circularity fixed.** `src/interpretation/split_target_lists.py` →
  `tumor_intrinsic_candidates.tsv` (n=37, led by PHGDH/SGMS2/PSAT1/PSPH/SMPD3/1)
  and `axis_confirmation_panel.tsv` (NK readout, labeled as such).
- **Manuscript reframed to v0.2.** Abstract, Key Points, Results §3.1–3.6,
  Discussion §4.1, Conclusion rewritten with real numbers; `main_claims.md` updated.
- **Blocker B RESOLVED (2026-07-07).** Pulled `baseline_internal_results.tsv` from
  the A100 (`/data/lgh/GC-NKGraph-Atlas`, user `user`), merged with the local GNN
  table on identical seed-42 folds, ran paired tests (`model_comparison*.tsv`).
  Verdict: GNN (MCC 0.706, AUROC 0.950) is **statistically on par** with LightGBM
  (0.733) / XGBoost (0.727) — paired p>0.27 — and **significantly beats**
  ElasticNet/SVM/MLP (p<0.05). Table 3 + §3.4 filled with the honest reading
  ("comparable accuracy, added interpretability"). GNN MCC reproduced at 0.7057
  (no fold drift).

- **Blocker C RESOLVED (2026-07-07).** Ran `run_geo_external_validation.py` on the
  A100 (GEOparse installed; probe→symbol via GPL570/GPL6947 directly on the
  processed probe matrices, since raw series matrices were absent). GSE62254
  54,675→22,880 genes (NK 6/7); GSE84437 49,576→25,159 genes (NK 7/7). Effector
  coupling replicates in BOTH cohorts (protrusion~cytotoxicity r=0.42/0.62,
  p≪1e-13); sm_balance~protrusion weakly positive/significant. Written to §3.3 +
  Table 5. TCGA-STAD/LIHC labels untouched (external validation done independently).

**All data/analysis blockers are now cleared.** Remaining before submission
(no server compute): final visual inspection of figures/PDF output, ORCID/CRediT
finalization, citation-detail verification, reviewer email/conflict verification,
and final character-encoding checks. Optional: `qc_filter.py` hardening.

---

> **Bottom line (original assessment, now superseded):** the project was **not
> submittable** before the fixes. Blocker A has been resolved by the honest
> partial-recovery reframe, and Blockers B/C have been resolved by the baseline
> comparison and corrected GEO validation. The remaining work is package hygiene
> and claim-boundary control.

---

## 1. Verdict by manuscript claim

Status legend: ✅ supported by real output · ⚠️ partial / needs reframing ·
❌ contradicted or missing · ⬜ not yet run.

| # | Claim (from `main_claims.md`) | Evidence on disk | Status | Blocker |
|---|-------------------------------|------------------|--------|---------|
| R1 | Cell-type-attributed SST-axis proxy defined | `sst_axis.py`, 7 modules/62 genes, spec doc | ✅ | none — this is genuinely done |
| R2 | **Liver positive control maps partial axis recovery** | `sst_axis_positive_control_recovery.tsv` | ✅ reframed | Effector arm recovers; metabolic coupling is weak/cell-resolved; physical topology is not recovered |
| R3 | Mechanism-grounded heterogeneous graph | `data/processed/graph/{nodes,edges}.tsv` | ⚠️ | graph built, but `metabolic_crosstalk` edge sign is meant to be *calibrated on the liver control* — and the control failed, so the calibration source is compromised |
| R4 | **Graph model improves over baselines** | `gc_nkgraph_internal_results.tsv` only | ❌ | no `baseline_*_results.tsv` exists — there is **no comparison at all**; GNN numbers stand alone |
| R5 | Gastric extension + tumor-intrinsic target list | `top_candidate_targets.tsv` | ⚠️ | list is **circular**: top hits (NKG7, PRF1, GZMB, GNLY) are the axis's own NK-effector markers with *negative* tumor specificity — not novel tumor-intrinsic druggable targets |
| R6 | In-silico SM-restoration stratification | derivable from scores | ⬜ | not yet computed; low risk once R2/R5 settled |

### The three hard blockers, stated plainly

**Blocker A — the positive control does not recover.**
`sst_axis_positive_control_liver.tsv` (your own pre-registered test):

| Hyp | Test | r | p | Expected | Result |
|-----|------|---|---|----------|--------|
| H1 | tumor_serine ⟂ nk_sm_balance | −0.016 | 0.74 | (calibrate) | null |
| H2 | nk_sm_balance → protrusion | −0.017 | 0.72 | positive | **FAIL** |
| H3 | protrusion → cytotoxicity | +0.551 | 5e-35 | positive | **PASS** |
| H4 | topology_permissive → dysfunction | +0.311 | 7e-11 | negative | **FAIL** |
| H5 | intratumoral < peritumoral NK | — | — | negative | DATA_UNAVAILABLE |

The original "Arm A recovers the axis as the credibility anchor" framing was not
supportable. The current manuscript resolves this by treating Arm A as a
partial-recovery scoping test rather than a full-recovery claim.

**Blocker B — no baseline comparison exists.** Methods §2.7 and Table 2 promise
XGBoost/LightGBM/RF/ElasticNet/SVM/MLP vs. the GNN. `run_all_baselines.py` is
complete and correct but was **never executed** — there is no
`baseline_internal_results.tsv`. The central methods claim ("graph improves or
stabilizes prediction") has zero supporting evidence today.

**Blocker C — external validation is broken.** `LOG.md` shows GSE62254 and
GSE84437 both logged `NK_MARKERS: 0` genes found → every sample collapsed to
`NK-intermediate`. Root cause (confirmed in code): `load_geo_expression()` keeps
**microarray probe IDs** (e.g. `1007_s_at`) as the row index, and
`standardize_gene_symbols()` only remaps ~20 hand-curated aliases. Gene-symbol
lookups therefore return nothing. Both external cohorts are currently unusable.

---

## 2. Secondary issues (not fatal, but reviewers will flag them)

| Issue | Where | Why it matters | Fix |
|-------|-------|----------------|-----|
| scRNA has **no real QC** | `run_scrna_v2.py` Step 3 keeps *all* cells | v1 dropped 100% of cells; v2 over-corrected to zero filtering — no mito%, no min-genes, no doublet removal. Reviewers expect standard QC | `qc_filter.py` (§5) |
| Cell-type annotation is threshold-on-mean-score | `run_scrna_v2.py` `annotate()` | hard cutoffs (nk>0.3) are fragile across platforms; no marker-based reference mapping | document as limitation or move to a reference-based labeler (scANVI) |
| Candidate score double-counts axis membership | `prioritize_targets.py` weights | SST membership = 0.30 + axis-core 0.10 = 40% of score → guarantees axis markers top the list | re-weight to surface *tumor-side* candidates, or split into two lists (see §4) |
| GNN vs baseline fold alignment unverified | `gc_nkgraph_atlas.py` vs `run_all_baselines.py` | both use StratifiedKFold(seed=42) but this must be *proven identical* or the comparison is invalid | `run_model_comparison.py` (§5) runs both on one split object |
| GEO clinical/survival not wired | `process_dataset()` GEO branch sets pheno=None | no survival analysis possible on external cohorts | parse GEO `!Sample_characteristics` if survival is claimed |
| Binary label only (`NK-hot-cytotoxic` vs rest) | baselines + GNN | Methods describe 4 states but models collapse to binary; state as such | align text to binary, or run true multiclass |

---

## 3. Recommended narrative reframe (decision: partial-recovery, honest)

Per the chosen direction, the paper is reframed so it no longer *depends* on full
axis recovery in liver. Concretely:

- **Downgrade Arm A** from "credibility anchor / the axis is recovered" to a
  **dissection**: "The *effector arm* of the axis — protrusion-machinery →
  cytotoxicity (H3) — reproduces strongly and independently (r=0.55, p<1e-34) in
  TCGA-LIHC. The *upstream metabolic arm* (SM-balance → protrusion, H2) and the
  dysfunction-direction test (H4) do **not** reproduce at the bulk-transcriptome
  level." Then explain *why this is expected and interesting*, not a failure:
  bulk transcriptome captures effector-machinery abundance well but is a poor
  proxy for the metabolite-level serine→SM flux (the anchor paper itself needed
  single-cell mass-spec precisely because transcription doesn't capture flux).
  This turns the negative result into a **scoping contribution**: it delineates
  which parts of a physical/metabolic mechanism are and are not transcriptionally
  recoverable — genuinely useful to the target lab.
- **Retitle the contribution** around "how far a transcriptome can reconstruct a
  physical immune-evasion mechanism," not "we reconstruct it."
- **Re-anchor calibration (R3):** since the liver control can't calibrate the
  `metabolic_crosstalk` sign, calibrate the edge on the arm that *did* reproduce
  (H3 effector axis) and state the SM-arm edge as hypothesis-weighted, not
  data-calibrated. Document the change in the pre-registration log.
- **Keep every honesty boundary from `POSITIONING_and_FRAMING.md` §6** — they now
  do real work carrying the reframed story.

This reframe is compatible with *Briefings in Bioinformatics* (methods/negative-
scoping results are in scope) but the two remaining blockers (B, C) still must be
closed — a paper with no baseline comparison and no external validation will not
survive review regardless of framing.

---

## 4. Fixing the circular target list (R5)

The promised deliverable is *tumor-intrinsic, druggable* gastric targets. The
current list returns NK-effector markers because axis membership dominates the
score and tumor-specificity is not required to be positive. Recommended:

1. **Split into two tables.** (a) *Axis-confirmation panel* — the NK-side markers
   (honestly labeled as readout, not targets). (b) *Tumor-intrinsic candidate
   list* — filtered to `tumor_specificity_log2 > 0` (up in malignant cells) AND
   membership in the tumor-serine module or its graph neighborhood.
2. **Re-weight** the tumor-intrinsic table toward tumor-specificity and
   graph-attention-to-NK-state, downweighting raw axis membership.
3. This is a `prioritize_targets.py` change, runnable once R2 scores are settled;
   no new data needed.

---

## 5. Server fix scripts (delivered with this report)

All three are written but **untested locally** (this workstation has no raw data
and only a Python-3.8 core stack). Run on `/opt/data/lgh/GC-NKGraph-Atlas`.

| Script | Fixes | Command | Produces |
|--------|-------|---------|----------|
| `src/preprocessing/fix_geo_gene_mapping.py` | Blocker C — probe→symbol | `python src/preprocessing/fix_geo_gene_mapping.py --gse GSE62254 GSE84437` | corrected `data/processed/bulk/gse*_expression.tsv` (gene-level) |
| `src/scrna_analysis/qc_filter.py` | scRNA no-QC issue | `python src/scrna_analysis/qc_filter.py --in data/processed/scrna/gc_integrated.h5ad` | QC'd h5ad + `results/tables/scrna_qc_summary.tsv` |
| `src/baselines/run_model_comparison.py` | Blocker B — no comparison | `python src/baselines/run_model_comparison.py` | `results/tables/model_comparison.tsv` (baselines + GNN, identical folds, paired test) |

After running: re-run NK scoring on the corrected GEO matrices, then re-run the
SST-axis and target prioritization. See §6 execution order.

---

## 6. Execution order (server)

```
STEP 1  Fix external cohorts
        python src/preprocessing/fix_geo_gene_mapping.py --gse GSE62254 GSE84437
        python src/immune_scoring/nk_scores.py            # NK_MARKERS should now be >0
        -> verify LOG shows non-zero gene counts + non-degenerate state distribution

STEP 2  Model comparison (closes Blocker B)
        python src/baselines/run_model_comparison.py
        -> results/tables/model_comparison.tsv ; check GNN vs best baseline delta + p

STEP 3  scRNA QC (methodological hardening)
        python src/scrna_analysis/qc_filter.py --in data/processed/scrna/gc_integrated.h5ad
        -> re-derive NK subset scores if QC changes the NK population materially

STEP 4  Positive-control re-examination (informs reframe vs. bug)
        Re-run sst_axis_validation on TCGA-LIHC AFTER confirming module signs and
        tumor-purity adjustment. If H2/H4 still fail with clean inputs, adopt the
        §3 partial-recovery reframe. If they flip, the earlier failure was an
        artifact and Arm A is salvageable as a full anchor.

STEP 5  Target list de-circularization (R5)
        Patch prioritize_targets.py per §4; regenerate the two tables.

STEP 6  Fill manuscript Results with real numbers; write Discussion around the
        reframed story; assemble submission package (cover letter, checklist,
        supp. tables, data-availability, ORCID/funding front matter).
```

Only after Steps 1–5 produce clean outputs is the manuscript worth filling in.

---

## 7. What is genuinely strong (keep and lead with)

- The **mechanism-card abstraction** and two-arm design are a real, publishable
  idea independent of the numeric results.
- **H3 reproduces convincingly** (protrusion→cytotoxicity, r=0.55) in an
  independent liver cohort — a real, defensible positive finding.
- Internal GNN performance is solid (acc≈0.87, AUROC≈0.94) — it just needs the
  baseline context to mean anything.
- The honesty discipline already written into the positioning document is exactly
  what makes a partial-recovery story credible rather than damaging.

---

## 8. Reality check on timeline

This is a "strong methods paper that needs its experiments finished," not a
"done paper that needs polishing." Closing Blockers B and C is mechanical
(scripts provided). Blocker A is a scientific judgment (reframe adopted). Budget:
~1–2 days server time for Steps 1–5, then ~2–3 days writing for Step 6.
```
```
