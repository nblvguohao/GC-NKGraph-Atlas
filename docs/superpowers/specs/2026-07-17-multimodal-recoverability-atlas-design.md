# Real Multimodal Mechanism-Recoverability Atlas

## Objective

Upgrade the manuscript from a single-mechanism transcriptomic boundary study
to a real-data, cross-mechanism atlas of which immune-evasion layers can be
recovered by transcriptomic, protein, metabolite, and spatial assays. The
deliverable is evidence about measurement boundaries, not a claim that a more
complex graph model improves classification.

## Scientific question and claim boundary

The primary question is whether recoverability is associated with biological
layer across independently authored cancer immune-evasion mechanisms. The four
cards are serine--sphingomyelin--membrane topology, TGF-beta--SMAD--NK
exclusion, adenosine--A2AR--cAMP/PKA suppression, and NKG2D--MICA/B shedding.

Every card-module result has exactly one of four statuses:

- `recovered`: pre-specified direction, FDR < 0.05, and concordant direction
  in at least two independent transcriptomic cohorts; or direct support from a
  real, modality-appropriate assay.
- `partially_recovered`: one cohort supports the direction, or an association
  survives controls but lacks the required cross-cohort support.
- `not_recovered`: tested with adequate gene/feature coverage but does not
  meet the above evidence threshold.
- `not_measured`: the needed real modality or analyte is unavailable. This is
  never converted to a negative biological conclusion.

A manuscript-level cross-mechanism statement is allowed only if at least
three of the four cards show the same ordered layer pattern in two or more
independent transcriptomic cohorts, and at least one real orthogonal modality
(protein, metabolite, or spatial assay) supports a relevant boundary. Otherwise
the manuscript reports a comparative atlas without claiming a universal rule.

## Evidence architecture

Each card is mapped to the following layer vocabulary:

| Layer | Transcriptomic evidence | Orthogonal evidence | Permitted interpretation |
|---|---|---|---|
| Upstream driver | bulk/scRNA module score | metabolite or protein measurement where available | machinery/capacity, not flux |
| Receptor/signaling | scRNA or bulk module score | protein abundance or phosphorylation | association, not receptor occupancy |
| Downstream effector | NK cytotoxicity/receptor module score | protein abundance where available | transcriptional effector state |
| Metabolite/protease execution | no RNA proxy may substitute | LC-MS metabolite or direct protease/readout | direct molecular measurement only |
| Spatial/physical phenotype | spatial module adjacency only | spatial coordinates or imaging | spatial co-localization, never cell-contact distance from Visium spots |

The existing TCGA-STAD, TCGA-LIHC, GSE62254, and GSE84437 matrices remain the
transcriptomic backbone. Existing real scRNA data provide cell-type attribution.
New public data are admitted only after a data-manifest check:

- GSE122401: 80 paired early-onset gastric tumors and adjacent tissues with
  public RNA-seq; linked proteomics is used only if a sample-level processed
  protein matrix is publicly obtainable.
- MTBLS3303: real LC-MS gastric tumor/adjacent-tissue metabolomics; it is used
  only if metabolite identities and sample groups are downloadable.
- GSE251950: 9 primary gastric cancers (10 Visium sections); it is used for
  spatial module adjacency/co-localization, not for inferred single-cell
  distances.
- Recent controlled-access multi-omics studies may be documented in the
  manifest but cannot enter results unless their processed matrices are openly
  downloadable and pass the same manifest checks.

## No-synthetic-data contract

Formal analyses and submission artifacts must never use synthetic, mock, demo,
or fallback data. A `real_data_manifest` records accession, source URL, assay,
species, sample count, file hashes, retrieval timestamp, license/access notes,
and the exact processed input file. The production loader rejects paths or
metadata containing `synthetic`, `mock`, or `demo`, rejects missing manifest
entries, and raises a clear error rather than generating a substitute dataset.

Small local fixtures remain permitted only under `tests/`. Their outputs are
prohibited from `results/` and `submission_bundle_BiB/`.

## Analysis design

1. Validate each source's species, sample-level units, case/control labels,
   feature identifiers, coverage of card modules, and checksum before analysis.
2. Run the same pre-specified module-scoring convention for each card and
   cohort. Do not tune gene sets or thresholds after seeing results.
3. Test card hypotheses with direction, effect size, nominal P value, BH-FDR,
   coverage, and an NK-lineage/purity sensitivity analysis where the endpoint
   is an NK program.
4. Summarize each card-layer across cohorts using the four-status vocabulary.
   Missing modality is `not_measured`, never `not_recovered`.
5. Run modality-specific analyses only for their direct endpoints: RNA-protein
   concordance on matched samples; metabolite tumor-vs-adjacent contrasts for
   named analytes; and spatial spot-level module adjacency in Visium data.
6. Generate a primary recoverability heatmap (card x layer) with evidence-type
   markers, a source/coverage table, and per-card result tables. The heatmap is
   descriptive; it does not pool incompatible effect sizes into a fabricated
   scalar score.

## Quality gates and failure behavior

- A non-human, cell-line-only, unpaired, inaccessible, or insufficiently
  annotated source is excluded with a recorded reason.
- A source lacking the needed named analyte is `not_measured` for that layer.
- Card hypotheses containing `NEEDS_REVIEW` directions remain calibration-only
  and cannot count toward recovery until direction is pre-specified from an
  independent positive-control rule.
- A result cannot say "physical validation" unless it uses a direct physical,
  protein, metabolite, or spatial readout appropriate to that endpoint.
- Model performance is secondary and cannot determine biological recovery.

## Manuscript integration

The abstract and conclusion will be changed only if the cross-mechanism gate
passes. Results will introduce the atlas after the individual Zheng/TGF-beta
analyses, report the gate verdict in one short paragraph, and place all
per-cohort tables, hashes, and coverage diagnostics in the supplement. The
main figure will show the card-by-layer evidence map and distinguish observed,
unmeasured, and gated layers visually.

## Non-goals

- No synthetic data generation or synthetic-data benchmark is part of this
  project.
- No full GCN, HGT, or Transformer implementation is required.
- No claim that public transcriptomic data directly measure metabolite flux,
  soluble MICA/B, phospho-SMAD/PKA, membrane topology, or cell-cell distance.
- No download of controlled-access data without both open processed matrices
  and a documented access basis.
