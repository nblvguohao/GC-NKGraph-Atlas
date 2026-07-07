# Phase 14R — Serine–Sphingomyelin–Topology transcriptional axis (SST-Axis)

> This phase **replaces** the old `Phase 14 — CROWN topology interface, GATED`.
> It converts the topology idea from a gated future-promise into a **real,
> computable, mechanism-grounded module**, while keeping the *physical* topology
> layer (imaging / single-cell mass-spectrometry) honestly gated.
>
> **Scientific anchor.** Zheng et al., *Tumors evade immune cytotoxicity by
> altering the surface topology of NK cells*, Nat Immunol 2023
> (DOI 10.1038/s41590-023-01462-9, PMID 36959292); and the follow-up framing of
> sphingomyelin as a metabolic immune checkpoint (Clin Transl Med 2023,
> DOI 10.1002/ctm2.1395). Reported mechanistic chain:
> tumor serine-metabolism dysregulation → reduced sphingomyelin (SM) in NK
> membranes → loss of membrane protrusions/microvilli → failure to form lytic
> immune synapses → loss of cytotoxicity; rescued by inhibiting SM catabolism
> (ASM/SMPD1, NSMASE1-3/SMPD2-4) and synergistic with Tim3 (HAVCR2) blockade.
>
> **Execution-order dependencies.** Runs AFTER Phase 3–6 (scRNA integration,
> NK atlas, NK-state scoring, trajectory). Feeds Phase 8 (graph construction),
> Phase 12 (ablation), and Phase 13 (candidate prioritization).

---

## 14R.0 Design split (read first)

This phase has two strictly separated layers. **Do not let them contaminate each other.**

```text
Layer 14R-A  — TRANSCRIPTIONAL AXIS  (status: ACTIVE / computable)
    What is computable from bulk + scRNA transcriptomes:
    the molecular MACHINERY / CAPACITY for the serine→SM→protrusion→cytotoxicity axis.
    This is the real scientific contribution of this phase.

Layer 14R-B  — PHYSICAL GROUND TRUTH  (status: GATED)
    Membrane protrusion density, microvilli density, membrane roughness,
    immune-synapse contact area (SEM / super-resolution imaging), and
    single-cell membrane SM content (single-immunocyte mass spectrometry).
    No real values currently available. Loader must raise a clear error if
    real targets are absent. Mock only for code-path tests, named *_MOCK_*,
    never entering scientific claims.
```

**Non-negotiable honesty rule for this whole phase:**
transcriptional proxy ≠ physical topology ≠ metabolite flux.
Never state that the model "predicts membrane topology" or "predicts SM content"
from transcriptome. The only permitted claim is that a transcriptional program
is **permissive-of / associated-with** the topology phenotype. Every output
table and figure caption produced here must carry this qualifier.

---

## 14R.1 Config — `configs/sst_axis_config.yaml`

```yaml
sst_axis:
  status: ACTIVE
  anchor_papers:
    - title: "Tumors evade immune cytotoxicity by altering the surface topology of NK cells"
      journal: "Nat Immunol"
      year: 2023
      doi: "10.1038/s41590-023-01462-9"
      pmid: "36959292"
    - title: "Sphingomyelin is a prospective metabolic immune checkpoint for natural killer cells"
      journal: "Clin Transl Med"
      year: 2023
      doi: "10.1002/ctm2.1395"

  # ---- Mechanistic gene modules ----
  # NOTE: gene sets below are STARTING sets and MUST pass marker_validation
  # (see scrna_config.annotation.marker_validation). Log any curation as NEEDS_REVIEW.
  modules:

    tumor_serine_capacity:            # TUMOR-CELL side (attribute to malignant cells)
      role: tumor_side
      cell_type_attribution: malignant
      genes: [PHGDH, PSAT1, PSPH, SHMT1, SHMT2, MTHFD1, MTHFD2, MTHFD1L, SLC1A4, SLC1A5]
      # expected_direction w.r.t. NK SM availability is a CROSSTALK term and is
      # NOT hard-coded. It must be CALIBRATED on the positive-control cohort
      # (see 14R.5). Default sign = NEEDS_REVIEW.
      expected_direction: NEEDS_REVIEW

    nk_sm_synthesis:                  # NK-CELL side
      role: immune_side
      cell_type_attribution: nk
      genes: [SGMS1, SGMS2]
      expected_direction: higher_is_more_topology_permissive

    nk_sm_catabolism:                 # NK-CELL side (drug targets; high = SM loss)
      role: immune_side
      cell_type_attribution: nk
      genes: [SMPD1, SMPD2, SMPD3, SMPD4]
      expected_direction: higher_is_less_topology_permissive

    nk_denovo_sphingolipid:           # serine-consuming de novo branch, NK side
      role: immune_side
      cell_type_attribution: nk
      genes: [SPTLC1, SPTLC2, SPTLC3, SPTSSA, CERS2, CERS4, CERS5, CERS6, DEGS1]
      expected_direction: context_dependent

    nk_protrusion_machinery:          # cortical actin / microvilli machinery, NK side
      role: immune_side
      cell_type_attribution: nk
      genes: [EZR, MSN, RDX,                      # ERM linkers
              ACTR2, ACTR3, ARPC1B, ARPC2, ARPC3, ARPC4, ARPC5,   # Arp2/3
              WAS, WASL, WASF1, WASF2, WASF3, WIPF1,              # WASP/WAVE NPFs
              CDC42, RAC1, RHOA,                  # Rho GTPases
              DIAPH1, DIAPH3, FMNL1,              # formins
              BAIAP2, PACSIN2]                    # I-BAR / membrane curvature
      expected_direction: higher_is_more_topology_permissive

    nk_synapse_cytotoxicity_outcome:  # OUTCOME anchors (NK side)
      role: outcome_anchor
      cell_type_attribution: nk
      genes: [NKG7, GNLY, GZMB, PRF1, IFNG,       # cytotoxicity
              LCP2, LAT, VAV1, TLN1, ITGAL, ITGB2] # synapse / adhesion
      expected_direction: axis_positive_correlate

    checkpoint_link:                  # therapeutic hook (NK side)
      role: therapeutic_hook
      cell_type_attribution: nk
      genes: [HAVCR2]                             # Tim3
      expected_direction: higher_is_less_topology_permissive

  # ---- Derived scores (deterministic) ----
  derived_scores:
    tumor_serine_capacity_score: mean_zscore(tumor_serine_capacity, cells=malignant)
    nk_sm_balance_score: "mean_zscore(nk_sm_synthesis) - mean_zscore(nk_sm_catabolism)"
    nk_protrusion_machinery_score: mean_zscore(nk_protrusion_machinery, cells=nk)
    nk_topology_permissive_score: >
      composite(nk_sm_balance_score, nk_protrusion_machinery_score)   # NK-intrinsic permissiveness
    sst_axis_score: >
      integrated(tumor_serine_capacity_score, nk_topology_permissive_score,
                 nk_synapse_cytotoxicity_outcome)   # sign of tumor term set by 14R.5 calibration

  # ---- Cell-type attribution (MANDATORY) ----
  attribution:
    require_celltype_resolution_before_axis_claim: true
    method: scrna_anchored          # attribute tumor vs NK modules using Phase 3-4 outputs
    bulk_deconvolution_fallback: [CIBERSORTx_or_local, quanTIseq_or_local]
    # If attribution is unavailable, axis scores may ONLY be reported as
    # whole-tissue mixtures and clearly labeled MIXED_UNRESOLVED. No tumor-vs-NK
    # crosstalk claim is permitted on unresolved data.

  # ---- Graph integration (feeds Phase 8) ----
  graph:
    new_node_types:
      - tumor_serine_program          # malignant-cell serine-capacity node
      - nk_topology_state             # NK SM/protrusion capacity node
    new_edge_types:
      - metabolic_crosstalk           # tumor_serine_program -> nk_topology_state (signed, weighted)
      - sm_topology_axis              # within-axis gene-gene edges
    rules:
      - metabolic_crosstalk edges require cell-type-resolved endpoints
      - crosstalk edge sign is CALIBRATED, never assumed; log as NEEDS_REVIEW until 14R.5 passes
      - no crosstalk edge derived from external-validation labels

  positive_control:
    cohort: TCGA-LIHC                 # liver: where the mechanism was proven
    scrna_intratumoral_vs_peritumoral_nk_dataset: SEARCH_REQUIRED
    # Task: check Data Availability of DOI 10.1038/s41590-023-01462-9 for deposited
    # scRNA/bulk of intratumoral vs peritumoral NK. If a real accession exists, use it.
    # DO NOT fabricate an accession. If none found, log DATA_UNAVAILABLE_lihc_nk.md.

  extension:
    cohort_group: gastric_cancer      # novel application (see Phase 1-2 GC cohorts)
```

---

## 14R.2 Tasks — Layer A (transcriptional axis, ACTIVE)

- [ ] Create `src/topology/sst_axis.py`.
- [ ] Load module gene sets from `configs/sst_axis_config.yaml`; run marker_validation; log curation to `results/logs/LOG.md` as `NEEDS_REVIEW`.
- [ ] **Attribute modules to cell types** using Phase 3–4 scRNA outputs:
  - tumor_serine_capacity → malignant cells only;
  - all NK-side modules → NK subset only.
- [ ] Compute per-cell (scRNA) and per-sample (bulk, cell-type-adjusted) derived scores.
- [ ] In scRNA: compare **intratumoral vs peritumoral/peripheral NK** for `nk_sm_balance_score` and `nk_protrusion_machinery_score` (this is the direct scRNA test of the Zheng phenotype at the molecular level).
- [ ] Fit all thresholds/standardization **on training folds only**; save; reuse for external cohorts (obey global leakage rules).
- [ ] Emit `results/tables/label_definition_sst_axis.md` documenting every scoring assumption.

Expected outputs (Layer A):

```text
data/processed/scrna/gc_nk_sst_scores_single_cell.tsv
results/tables/sst_axis_scores_bulk.tsv
results/tables/sst_axis_module_scores.tsv
results/tables/nk_intratumoral_vs_peritumoral_sst.tsv
results/tables/label_definition_sst_axis.md
results/figures/fig9a_sst_axis_scrna.pdf
results/figures/fig9b_sst_axis_bulk.pdf
```

---

## 14R.3 Tasks — Layer B (physical ground truth, GATED)

- [ ] Keep `configs/topology_schema.yaml` (physical targets) as the GATED schema.
- [ ] `src/topology/topology_infer.py`: real loader MUST raise a clear error if physical/MS targets are absent.
- [ ] Provide a `*_MOCK_*` dataset **only** for code-path testing.
- [ ] Document in README exactly how real paired data (SEM protrusion metrics, single-immunocyte MS SM values) would be formatted, and how they would attach to `nk_topology_state` graph nodes.

Expected outputs (Layer B):

```text
configs/topology_schema.yaml
src/topology/topology_infer.py
results/logs/TOPOLOGY_GATED_STATUS.md
tests/test_topology_mock_pipeline.py
```

Hard rules (Layer B):

```text
No physical topology or MS-SM value may enter any scientific table or figure.
Real loader raises error if real targets are missing.
Mock files are named *_MOCK_* and are excluded from all claims.
```

---

## 14R.4 Graph edges written for Phase 8

Add to `data/processed/graph/edges.tsv` (only after cell-type attribution passes):

```text
edge_type = metabolic_crosstalk
    source = tumor_serine_program node
    target = nk_topology_state node
    weight = |calibrated association| ; direction = calibrated sign (14R.5)
    evidence = "Zheng 2023 NatImmunol serine->SM->topology axis (transcriptional proxy)"

edge_type = sm_topology_axis
    within-axis gene-gene edges (SM synthesis/catabolism/protrusion machinery)
    weight = coexpression within NK subset (train folds only)
```

---

## 14R.5 Positive-control validation protocol (Arm A — LIVER)

**Purpose.** Before making any gastric-cancer claim, demonstrate that the framework
**recovers the published Zheng mechanism** in the system where it was proven (liver).
This is the credibility anchor of the whole project.

Pre-register the following hypotheses **and their expected directions** in
`manuscript/notes/sst_axis_prereg.md` *before* looking at test data:

```text
H1 (crosstalk):  tumor_serine_capacity_score  ⟂  nk_sm_balance_score
                 (direction to be CALIBRATED on Arm A; report the learned sign)
H2 (SM→machinery): nk_sm_balance_score  (+)  nk_protrusion_machinery_score
H3 (machinery→function): nk_protrusion_machinery_score  (+)  cytotoxicity anchors
H4 (dysfunction): nk_topology_permissive_score  (−)  HAVCR2 / dysfunction score
H5 (scRNA phenotype): intratumoral NK show LOWER sm_balance + protrusion machinery
                 than peritumoral/peripheral NK
```

Data:

```text
Bulk liver:  TCGA-LIHC (+ optional GEO HCC cohorts, no fabricated accessions).
scRNA liver: intratumoral-vs-peritumoral NK dataset from the Nat Immunol paper's
             Data Availability if it exists (SEARCH_REQUIRED); else public HCC
             scRNA with NK subset; else DATA_UNAVAILABLE_lihc_nk.md.
```

Recovery metric:

```text
For H2-H5, report effect sizes + CIs across seeds/folds and state PASS/FAIL vs
pre-registered direction. For H1, report the CALIBRATED sign and its stability.
"Recovery" = H2-H5 pass in the pre-registered direction in liver.
```

Outputs:

```text
manuscript/notes/sst_axis_prereg.md
results/tables/sst_axis_positive_control_liver.tsv
results/figures/fig10_positive_control_liver.pdf
```

Acceptance for Arm A:

```text
- H2-H5 recovered in pre-registered direction in liver, OR
- an explicit, logged explanation of which failed and why (NEEDS_REVIEW).
- The calibrated tumor-serine crosstalk sign from H1 is the sign used for all
  downstream metabolic_crosstalk edges and sst_axis_score. It is NOT assumed.
```

---

## 14R.6 Extension (Arm B — GASTRIC)

- [ ] Repeat 14R.2–14R.5 hypotheses in gastric cancer cohorts (Phase 1–2).
- [ ] Report whether the same axis operates in GC (this is the **novel** result).
- [ ] Feed gastric SST-axis genes into the candidate pool (Phase 13) with an
      `sst_axis_membership` evidence column.

Outputs:

```text
results/tables/sst_axis_gastric.tsv
results/figures/fig11_sst_axis_gastric.pdf
```

---

## 14R.7 Candidate-prioritization hooks (feeds Phase 13)

Add these columns to `candidate_evidence_matrix.tsv`:

```text
sst_axis_membership          # {tumor_serine, nk_sm_synthesis, nk_sm_catabolism,
                             #  nk_denovo_sphingolipid, nk_protrusion_machinery, none}
sst_axis_direction_consistent # does the gene's sample-level behavior match calibrated axis?
positive_control_recovered   # did this gene's axis behavior replicate in liver?
```

Add the wet-lab-actionable assays used by the anchor lab to
`recommended_validation_assay` vocabulary:

```text
NK-tumor co-culture cytotoxicity assay
SEM / super-resolution membrane-protrusion imaging
single-immunocyte mass spectrometry (membrane SM)
lytic immune-synapse formation assay
sphingomyelinase-inhibitor rescue (+/- Tim3/HAVCR2 blockade)
```

---

## 14R.8 Optional: in-silico perturbation / stratification head (feeds Phase 10 & 13)

If the graph model supports it, add a perturbation readout that mirrors the
anchor lab's therapeutic logic:

```text
For each sample, estimate predicted change in NK cytotoxicity / topology_permissive
score under an in-silico "SM-restoration" perturbation
(raise nk_sm_synthesis / suppress nk_sm_catabolism).

Output a stratification table ranking samples by predicted benefit from
"SM-restoration + Tim3 blockade" combination logic.
```

```text
results/tables/sst_axis_insilico_perturbation.tsv
results/tables/sst_axis_combination_stratification.tsv
```

**Honesty rule.** This is a model-based hypothesis for patient stratification,
NOT a validated predictor. Label every output accordingly.

---

## 14R.9 Tests (extend Phase 17)

```text
tests/test_sst_axis_celltype_attribution.py
    - tumor_serine module scored only on malignant cells
    - NK-side modules scored only on NK subset
    - MIXED_UNRESOLVED flag raised when attribution unavailable

tests/test_sst_axis_leakage.py
    - axis thresholds fit on train folds only
    - external cohorts use saved thresholds
    - no metabolic_crosstalk edge derived from external-validation labels

tests/test_sst_axis_gating.py
    - no physical topology / MS-SM value appears in scientific outputs
    - real topology loader raises error when real targets absent
    - crosstalk edge sign is present and equals the 14R.5-calibrated sign
```

---

## 14R.10 Acceptance criteria (add to Section 18)

```text
[ ] SST-axis modules validated and scored with cell-type attribution.
[ ] scRNA intratumoral-vs-peritumoral NK comparison generated (Zheng phenotype, molecular level).
[ ] Positive-control (liver) recovery reported with PASS/FAIL vs pre-registration.
[ ] Tumor-serine crosstalk sign CALIBRATED (not assumed) and reused downstream.
[ ] Gastric-cancer extension reported as the novel result.
[ ] Physical topology + MS-SM remain GATED; no mock values in claims.
[ ] Every SST output carries the "transcriptional proxy, not physical topology" qualifier.
```
