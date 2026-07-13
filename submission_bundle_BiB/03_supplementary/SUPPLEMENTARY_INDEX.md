# Supplementary Information — Index

Manuscript: *Transcriptional Reach of the Serine–Sphingomyelin–
Membrane-Topology Axis in NK-Cell Immune Evasion* (GC-NKGraph-Atlas).

All supplementary tables are tab-separated (`.tsv`) unless noted. Column meanings
follow the manuscript Methods (§2). Values are the final results reported in the
main text and figures.

## Supplementary Methods (to add as prose before submission)
A short Supplementary Methods note should accompany these tables, covering: (i)
NK-state label definition and thresholds (see `label_definition.md`,
`nk_state_thresholds.json`); (ii) SST-axis module gene lists and scoring
(configs/mechanism_cards/zheng_nk_sm_topology.yaml); (iii) scRNA QC thresholds and
scVI/Leiden settings (§2.4); (iv) the paired-test protocol for model comparison
(§2.7). These are fully specified in the repository; the note simply points to them.

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
| `ablation_edge_types.tsv` | Table / §3.7 | Edge-type ablation: performance with each heterograph edge type removed. |
| `ablation_results.tsv` | §3.7 | Full ablation metrics per configuration. |
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
| `domain_baselines_per_fold.tsv` / `domain_baselines_summary.tsv` / `domain_baselines_tests.tsv` | §3.4 | NK-marker-signature and SST-module-signature baselines (per-fold results, summary, and paired significance tests against the GNN). |
| `target_validation_depmap.tsv` / `target_validation_nk_state_de.tsv` / `target_validation_v2_merged.tsv` | §3.5 / Table 4 | DepMap 26Q1 CERES essentiality, NK-state DE results, and the merged evidence table underlying the evidence-tiered candidate list. |
| `nk_state_de_external_replication.tsv` / `nk_state_de_external_concordance.tsv` / `nk_state_de_external_replication_summary.md` | §3.5 | External replication of the NK-state DE test in GSE62254/GSE84437 and directional concordance with TCGA-STAD. |
| `trivial_baseline_comparison.tsv` / `trivial_baseline_overlap.tsv` / `trivial_baseline_summary.md` | §3.5 | Comparison of the five-dimension target-scoring against a mechanism-card-membership-only trivial baseline. |
| `geneset_separation_audit.tsv` / `geneset_separation_audit_summary.md` | §4.3 (Limitations item 15) | Gene-set overlap audit between the NK-state classification label, the SST-axis modules, and the 37-gene candidate list (label-leakage / circularity check). |

*Note:* the interactive web playground (`web/index.html`, also hosted via GitHub
Pages) lets reviewers browse the mechanism cards and the candidate list with the
external-database evidence linked above.
