# Supplementary Information — Index

Manuscript: *Transcriptional Reach of the Serine–Sphingomyelin–
Membrane-Topology Axis in NK-Cell Immune Evasion* (GC-NKGraph-Atlas).

All supplementary tables are tab-separated (`.tsv`) unless noted. Column meanings
follow the manuscript Methods (§2). Values are the final results reported in the
main text and figures.

## Supplementary Methods

These notes expand the main-text Methods (§2) for the parameters most relevant to
reproduction. All settings are fixed in the repository
(`nblvguohao/GC-NKGraph-Atlas`); the pointers below identify the exact files.

**S.M.1 NK immune-state label definition and thresholds** (main text §2.5;
`label_definition.md`, `nk_state_thresholds.json`). Each sample receives four
transcriptional scores computed as the mean z-score of a literature-derived gene
set: an NK-infiltration score (NCAM1, NCR1, KLRD1, KLRK1, KLRF1, NKG7, GNLY, GZMB,
PRF1, FCGR3A, XCL1, XCL2, CCL5, IFNG, TYROBP), an NK-cytotoxicity score (NKG7,
GNLY, GZMB, PRF1, IFNG, XCL1, XCL2, CCL5), a dysfunction score (mean-z of the
dysfunction set minus the cytotoxicity set), and an exclusion score (mean-z of a
CAF/ECM set minus the infiltration score). Samples are assigned to one of four
states — NK-hot-cytotoxic, NK-hot-dysfunctional, NK-cold/excluded, or
NK-intermediate — by thresholding these scores. Thresholds are the **median values
of the primary training cohort (TCGA-STAD)** (`nk_state_thresholds.json`:
infiltration 0.0388, cytotoxicity 0.00032, dysfunction 0.0425, exclusion 0.0214);
external cohorts are labeled with these same saved training thresholds, never
re-fit per cohort. These signatures are transcriptional proxies, not functional
measurements, and the gene sets are starting points (see main-text Limitations).

**S.M.2 SST-axis modules and scoring** (main text §2.3;
`configs/mechanism_cards/zheng_nk_sm_topology.yaml`). The
serine–sphingomyelin–topology axis is operationalized as seven gene modules
(~100 genes total) derived from the anchor paper and its follow-up: tumor serine
capacity (10 genes), NK SM synthesis (2), NK SM catabolism (4), NK de-novo
sphingolipid (9), NK protrusion machinery (24), NK synapse/cytotoxicity outcome
(11), and a checkpoint link (HAVCR2). Per-cell (or per-sample) module scores are
the mean z-score of the module's constituent genes. Three derived composites are
defined: `nk_sm_balance` = mean-z(SM synthesis) − mean-z(SM catabolism);
`nk_topology_permissive`, combining SM balance and protrusion machinery; and
`sst_axis_score`, integrating tumor serine capacity (sign calibrated on the liver
control), topology-permissive score, and cytotoxicity outcome. Per the manuscript's
honesty rule, these quantify a *transcriptional program associated with* the
topology phenotype and are never interpreted as predicting physical membrane
topology or SM metabolite content.

**S.M.3 scRNA-seq QC, integration, and clustering** (main text §2.4;
`src/scrna_analysis/run_scrna_v2.py`). Single-cell data (GSE246662) comprise nine
samples across healthy liver, gastric cancer, and gastric-cancer liver metastasis
(166,829 cells). Per-sample matrices are orientation-auto-detected and concatenated
on the gene intersection (inner join). Per-cell QC metrics (detected genes,
mitochondrial fraction) are computed but **no hard per-cell threshold or
doublet-based exclusion is applied**; all 166,829 cells proceed to normalization.
Technical variance (library size, detected-gene count) is instead handled
downstream by the count-depth residualization diagnostics (§3.2). Counts are
library-size-normalized to 10⁴ and log1p-transformed; 3,000 highly variable genes
are selected (Seurat v3 flavor, variance-based fallback). Batch effects across the
nine samples are corrected with **scVI** (`sample_id` as batch key, 30 latent
dimensions, 2 layers, up to 200 epochs with early stopping). Neighbors, UMAP, and
**Leiden clustering (resolution 1.0)** are computed on the scVI latent space
(SCANPY). NK cells are separated from T cells by a marker-threshold rule (NK score
high, T score low), yielding **8,310 NK cells** for the axis analyses
(per-sample counts in `gc_scrna_dataset_summary.tsv`).

**S.M.4 Model-comparison paired-test protocol** (main text §2.6–§3.4;
`model_comparison.tsv`, `model_comparison_stats.tsv`). The GNN and all six tabular
baselines are trained and evaluated on the **same 5-fold stratified
cross-validation splits**. For each pairwise contrast (GNN vs each baseline) the
per-fold metric vectors (MCC, AUROC) are compared with both a paired *t*-test and a
Wilcoxon signed-rank test across the five folds; both statistics are reported so
that the conclusion does not rest on the parametric assumption. The NK-state
classifier hyperparameters (learning rate 1.7×10⁻³, weight decay 5.6×10⁻⁶, dropout
0.6) were selected by a 100-trial Bayesian (TPE) search maximizing MCC
(`gc_nkgraph_bayesian_trials.tsv`, `gc_nkgraph_best_hyperparams.tsv`).

## Supplementary Tables

| File | Maps to | Content |
|------|---------|---------|
| `sst_axis_positive_control_recovery.tsv` | Table 2 / Fig 1 | Arm A per-hypothesis (H1–H5) outcomes at bulk and single-cell resolution (r, p, expected sign, recovery verdict). |
| `axis_confirmation_panel.tsv` | §3.2 | Consolidated axis-confirmation statistics across resolutions. |
| `sst_axis_condition_comparison.tsv` | Table / Fig 2 | SST-module scores: healthy-liver vs gastric-cancer NK, per module, with p-values. |
| `sst_axis_gc_vs_hl.tsv` | §3.3 | Gastric-cancer vs healthy-liver NK contrasts (Arm B). |
| `sst_axis_scrna_by_tissue.tsv` | §3.2–3.3 | Single-cell SST-module scores by tissue (intratumoral / normal / metastatic). |
| `external_validation_results.tsv` | Table 5 / Fig 2 | Effector- and metabolic-arm replication in GSE62254 and GSE84437 (r, p, NK-marker coverage). |
| `gc_scrna_dataset_summary.tsv` | Table 1 | scRNA dataset composition (cells, NK cells, clusters, tissues). |
| `baseline_internal_results.tsv` | Table 3 | Six tabular baselines, 5-fold CV metrics. |
| `gc_nkgraph_internal_results.tsv` | Table 3 | GC-NKGraph-Atlas (GNN) 5-fold CV metrics. |
| `model_comparison.tsv` | Table 3 / Fig 4 | Per-fold metrics for all methods (GNN + baselines). |
| `model_comparison_summary.tsv` | Table 3 | Mean ± SD summary across folds. |
| `model_comparison_stats.tsv` | §3.4 | Paired significance tests (GNN vs each baseline; t-test + Wilcoxon). |
| `ablation_results.tsv` | Table (§3.7) | Graph ablation on the enriched real graph (FULL / −MC / −SST): edge counts, embedding-coupling H1/H2, and modularity. |
| `tumor_intrinsic_candidates.tsv` | Table 4 / Fig 3 | De-circularized tumor-intrinsic candidate list (n=37): rank, gene, category, score, tumor-specificity, druggability stage, recommended assay. |
| `candidate_evidence_matrix.tsv` | Fig 3 | Multi-dimensional evidence matrix underlying candidate scoring (DepMap / DrugBank / Open Targets). |
| `label_definition.md` | §2.5 | NK-state label definitions. |
| `nk_state_thresholds.json` | §2.5 | Numeric thresholds used to assign NK states. |
| `h3_scoring_method_diagnostic.tsv` / `h3_scoring_method_diagnostic_summary.md` | §3.2 | H3 count-depth/scVI-latent residualization and expression-matched permutation-null diagnostics (the "does the single-cell H3 number survive confound control" analysis). |
| `h3_activation_control.tsv` / `h3_activation_matched_subset.tsv` | §3.2 | Generic NK-activation partialling and activation-matched-quintile control for the single-cell H3 coupling. |
| `h3_module_permutation_test.tsv` | §3.2 | 10,000-draw random-gene-module permutation null for H3 (naive mean-zscore version). |
| `h3_leave_one_sample_out.tsv` | §4.3 (Limitations item 5) | Leave-one-sample-out sensitivity analysis for the pooled H3 meta-analytic estimate. |
| `sst_axis_count_depth_control.tsv` | §3.2 | Per-cell library-size / detected-gene-count covariates underlying the count-depth residualization. |
| `sst_axis_pseudoreplication_corrected.tsv` | §2.9 / Table 2 | Naive per-cell values retained for transparency alongside the corrected meta-analytic values reported in-text. |
| `mc_edge_sign_calibration_audit.tsv` | §2.5 | Audit distinguishing the fixed `metabolic_crosstalk` edge weight (uncalibrated) from the separately-calibrated `sst_axis_score` sign term. |
| `mechanism_card_comparison.tsv` / `mechanism_card_gene_overlap.tsv` | §4.2 | Cross-card comparison and pairwise gene-set Jaccard overlap for the four registered mechanism cards. |
| `bulk_h3_purity_control.tsv` | §3.2 / §4.3 (Limitations item 17) | Bulk effector-arm (protrusion~cytotoxicity) correlation before vs after partialling out a clean NK-lineage fraction proxy, across TCGA-LIHC / GSE62254 / GSE84437 (zero-order r, partial r, 95% CI, attenuation %). Shows the bulk coupling is ~50% NK-abundance-driven and does not survive in GSE84437. |
| `mechanism_card_tgfb_recovery.tsv` | §4.1 / §4.2 | Second mechanism card (TGFβ→SMAD→NK exclusion) run end-to-end on the gastric cohorts: pre-registered hypotheses H2–H5, zero-order and NK-fraction-controlled H3, per cohort. Demonstrates the transcriptional-reach boundary generalizes across two mechanisms. |
| `domain_baselines_per_fold.tsv` / `domain_baselines_summary.tsv` / `domain_baselines_tests.tsv` | §3.4 | NK-marker-signature and SST-module-signature baselines (per-fold results, summary, and paired significance tests against the GNN). |
| `target_validation_depmap.tsv` / `target_validation_nk_state_de.tsv` / `target_validation_v2_merged.tsv` | §3.5 / Table 4 | DepMap 26Q1 CERES essentiality, NK-state DE results, and the merged evidence table underlying the evidence-tiered candidate list. |
| `nk_state_de_external_replication.tsv` / `nk_state_de_external_concordance.tsv` / `nk_state_de_external_replication_summary.md` | §3.5 | External replication of the NK-state DE test in GSE62254/GSE84437 and directional concordance with TCGA-STAD. |
| `trivial_baseline_comparison.tsv` / `trivial_baseline_overlap.tsv` / `trivial_baseline_summary.md` | §3.5 | Comparison of the five-dimension target-scoring against a mechanism-card-membership-only trivial baseline. |
| `geneset_separation_audit.tsv` / `geneset_separation_audit_summary.md` | §4.3 (Limitations item 15) | Gene-set overlap audit between the NK-state classification label, the SST-axis modules, and the 37-gene candidate list (label-leakage / circularity check). |
| `gc_nkgraph_bayesian_trials.tsv` / `gc_nkgraph_best_hyperparams.tsv` | §2.6 | 100-trial Bayesian (TPE) hyperparameter search log and selected configuration (learning rate, weight decay, dropout, edge weights) for the NK-state classifier. |

*Note:* the interactive web playground (`web/index.html`, also hosted via GitHub
Pages) lets reviewers browse the mechanism cards and the candidate list with the
external-database evidence linked above.
