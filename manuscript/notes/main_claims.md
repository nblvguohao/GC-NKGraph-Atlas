# Main Claims — GC-NKGraph-Atlas

> Each claim maps to a figure/table in `manuscript/main_manuscript.md`.
> Checked claims have been verified by pipeline output or curated submission tables.
> This file records the current submission-facing claim boundary after the
> pre-submission reviewer audit.

---

## Claim 1 — R1. Transcriptional proxy definition
**Statement:** We define a reproducible, cell-type-attributed transcriptional proxy for the serine → sphingomyelin → membrane-protrusion-machinery → cytotoxicity axis of NK immune evasion.

| Item | Status |
|------|--------|
| SST-axis gene modules defined (7 modules, 62 genes) | ✅ Done |
| Module-to-cell-type attribution (malignant vs NK) | ✅ Done |
| Claim boundary documented (transcriptional proxy ≠ physical topology) | ✅ Done |
| Figure/Table | Table in Methods 2.3, Figure 1 |
| Pipeline verification | `src/topology/sst_axis.py` computes all modules |

---

## Claim 2 — R2. Liver positive control: PARTIAL recovery (reframed, honest)
**Statement (revised):** In liver cancer, the framework recovers the **effector
layer** of the axis (protrusion→cytotoxicity) at both bulk and single-cell
resolution, detects a **weak metabolic coupling** (SM-balance→protrusion) only
when cell-type is resolved, and does **not** recover the **physical topology
phenotype** from machinery transcription. Reported as a scoping map, not full
recovery.

| Item | Status |
|------|--------|
| TCGA-LIHC expression processed (n=423) | ✅ Done |
| H1–H5 pre-registered + tested (`sst_axis_positive_control_recovery.tsv`) | ✅ Done |
| H2: sm_balance (+) protrusion | ⚠️ null in bulk; weak but significant in single-cell NK (r=+0.030, p=6e-3) |
| H3: protrusion (+) cytotoxicity | ✅ **recovered** (bulk r=0.55; scNK r=0.32) |
| H4: topology_permissive (−) dysfunction/HAVCR2 | ❌ not recovered (wrong sign both resolutions) |
| H5: intratumoral < normal NK | ⚠️ cytotoxicity **recovered** (Δ=−0.14, p=6e-52); protrusion wrong sign |
| Recovery verdict | ✅ Reframed: effector recovered; cell-resolved metabolic signal weak/partial; topology not recovered |
| Figure/Table | Figure 1, Table 2 |

---

## Claim 3 — R3. Mechanism-grounded heterogeneous graph
**Statement:** We construct a single-cell-informed tumor–NK heterogeneous graph in which a mechanism-grounded `metabolic_crosstalk` edge is justified by the biology, not by generic priors.

| Item | Status |
|------|--------|
| 6 edge types specified | ✅ Done |
| metabolic_crosstalk edge connects tumor_serine ↔ NK topology | ✅ Done |
| Edge sign calibrated on liver control (not hard-coded) | ✅ Designed, ⬜ Executed |
| No edge from external-validation labels | ✅ Enforced in code |
| Cell-type-resolved endpoints required | ✅ Enforced in config rules |
| Figure/Table | Methods 2.5 |

---

## Claim 4 — R4. Graph model: on par with top baselines + interpretability (reframed)
**Statement (revised):** On internal CV the graph model is **statistically on par**
with the strongest gradient-boosting baselines and **significantly better** than
linear/kernel/shallow baselines, while additionally providing a
mechanism-structured gene embedding. We do **not** claim SOTA accuracy.

| Item | Status |
|------|--------|
| 6 baselines implemented + run on A100 (`baseline_internal_results.tsv`) | ✅ Done |
| GNN 5-fold CV on TCGA-STAD | ✅ Done (Acc 0.864, MCC 0.706, AUROC 0.950) |
| GNN vs baselines, identical folds + paired tests (`model_comparison*.tsv`) | ✅ Done |
| Verdict: GNN ≈ LightGBM(0.733)/XGBoost(0.727) MCC (p>0.27); > ElasticNet/SVM/MLP (p<0.05) | ✅ Done |
| External validation on GSE62254/GSE84437 | ✅ Done — effector coupling replicates (r=0.42/0.62, p≪1e-13) after GEO probe fix (`external_validation_results.tsv`) |
| Figure/Table | Table 3, Table 5 |

---

## Claim 5 — R5. Gastric cancer extension + target list
**Statement:** The recoverable effector layer is tested in gastric cancer; where it holds, we prioritize putative tumor-intrinsic candidate targets supported by multi-evidence patterns.

| Item | Status |
|------|--------|
| Gastric cancer axis tested | ✅ Done: effector coupling replicated in GSE62254/GSE84437 (`external_validation_results.tsv`) |
| Candidate pool assembled (SST genes + seed candidates) | ✅ Done |
| Multi-evidence matrix (5 dimensions) | ✅ Done |
| Composite target score with weighting | ✅ Done |
| Circularity fixed: split tumor-intrinsic vs axis-readout (`split_target_lists.py`) | ✅ Done |
| Putative tumor-intrinsic candidates (n=37) with assays; led by PHGDH/SGMS2/PSAT1/PSPH/SMPD3 | ✅ Done (`tumor_intrinsic_candidates.tsv`) |
| Gold-standard druggable enzymes surfaced at top | ✅ Done (PHGDH, SMPD1/3, PSAT1) |
| Figure/Table | Table 4, Figure 3 |

---

## Claim 6 — R6. In-silico SM-restoration stratification (secondary hypothesis)
**Statement:** A model-based in-silico "SM-restoration" readout is retained only as a hypothesis for experimental testing, not as a validated predictor or primary manuscript claim.

| Item | Status |
|------|--------|
| Readout defined (nk_sm_catabolism + HAVCR2) | ✅ Designed |
| Per-sample stratification | ⬜ Not promoted to a primary claim before experimental validation |
| Qualifier language enforced | ✅ |
| Figure/Table | Optional / not required for the current submission package |

---

## Execution checklist

| Phase | Script | Output | Status |
|-------|--------|--------|--------|
| DATA | `download_tcga.py` | TCGA-STAD + TCGA-LIHC in `data/raw/` | ⬜ |
| DATA | `download_geo.py` | GSE62254/84437 in `data/raw/` | ⬜ |
| PREPROCESS | `run_bulk_preprocessing.py` | Standardized expr in `data/processed/` | ⬜ |
| scRNA | `run_scrna_pipeline.py` | NK atlas + states in `results/` | ⬜ |
| SST | `sst_axis.py` | SST scores in `results/tables/` | ⬜ |
| GRAPH | `build_heterograph.py` | `data/processed/graph/{nodes,edges}.tsv` | ⬜ |
| BASELINES | `run_all_baselines.py` | `results/tables/baseline_internal_results.tsv` | ⬜ |
| MODEL | `gc_nkgraph_atlas.py` | `results/tables/gc_nkgraph_gnn_internal_results.tsv` | ⬜ |
| TARGETS | `prioritize_targets.py` | `results/tables/candidate_evidence_matrix.tsv` | ⬜ |
