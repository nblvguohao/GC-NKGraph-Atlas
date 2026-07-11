# Supplementary Information — Index

Manuscript: *Mapping the Transcriptional Reach of the Serine–Sphingomyelin–
Membrane-Topology Axis of NK-Cell Immune Evasion* (GC-NKGraph-Atlas).

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

*Note:* the interactive web playground (`web/index.html`, also hosted via GitHub
Pages) lets reviewers browse the mechanism cards and the candidate list with the
external-database evidence linked above.
