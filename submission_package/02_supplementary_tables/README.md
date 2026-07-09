# Supplementary Tables — Index

> This document maps each supplementary table to its manuscript section.
> Last updated: 2026-07-09.

---

## File inventory

### S1. Dataset summaries

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `gc_scrna_dataset_summary.tsv` | Per-sample scRNA-seq QC summary: cells, NK count, tissue, condition | Methods §2.4, Table 1 |
| `nk_state_thresholds.json` | NK immune state classification z-score thresholds | Methods §2.4 |

### S2. SST-axis results (Arm A — liver positive control)

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `sst_axis_positive_control_recovery.tsv` | H1–H5 outcomes: r, p, direction, verdict per hypothesis/resolution | §3.2, Table 2, Figure 1 |
| `sst_axis_condition_comparison.tsv` | Cross-condition SST module score comparison (tumor vs normal) | §3.2 |
| `sst_axis_gc_vs_hl.tsv` | Gastric cancer vs healthy liver NK SST module scores | §3.3 |
| `sst_axis_scrna_by_tissue.tsv` | Per-tissue SST module score means in single NK cells | §3.3, Figure 2 |

### S3. Model comparison

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `model_comparison.tsv` | 7 methods × 5 folds: all metrics, identical seed-42 splits | §3.4, Table 3, Figure 4 |
| `model_comparison_stats.tsv` | Paired tests: GNN vs each baseline, ΔMCC, Wilcoxon/t-test p | §3.4, Table 3 footnote |
| `model_comparison_summary.tsv` | Mean ± SD over folds per method × metric | §3.4 |
| `gc_nkgraph_internal_results.tsv` | GNN per-fold results | §3.4 |
| `baseline_internal_results.tsv` | Six baselines per-fold results | §3.4 |

### S4. Candidate targets

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `tumor_intrinsic_candidates.tsv` | 37 tumor-intrinsic candidates: rank, gene, scores, druggability, assay | §3.5, Table 4, Figure 3 |
| `axis_confirmation_panel.tsv` | 36 NK-side genes (readout, not targets) | §3.5 |
| `candidate_evidence_matrix.tsv` | Full multi-evidence matrix with raw and normalized scores | §2.8, §3.5 |

### S5. External validation

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `external_validation_results.tsv` | GSE62254 + GSE84437: axis correlations, NK marker coverage | §3.3, Table 5, Figure 2 |

### S6. Label definitions

| File | Description | Manuscript reference |
|------|-------------|---------------------|
| `label_definition.md` | NK state definitions: gene sets, scoring, threshold calibration | Methods §2.4 |

---

## Column legend

| Column | Meaning |
|--------|---------|
| `r` | Pearson correlation coefficient |
| `p` | Two-sided p-value |
| `MCC` | Matthews correlation coefficient |
| `AUROC` | Area under ROC curve |
| `AUPRC` | Area under precision–recall curve |
| `tumor_specificity_log2FC` | log2 FC malignant vs non-malignant (scRNA) |
| `druggability` | Drug development stage |
| `composite_score` | Weighted 5-dimension sum (§2.8) |
| `MIXED_UNRESOLVED` | Cell-type attribution unavailable |

---

## Supplementary Methods

### Probe-to-gene remapping for external validation

GSE62254 (GPL570, 54,675 probes) and GSE84437 (GPL6947, 49,576 probes) were
remapped to gene symbols using platform GPL SOFT annotation files from GEO.
Multi-probe mapping: the probe with highest mean expression was retained
(max-mean). NK markers recovered: 6/7 (GSE62254) and 7/7 (GSE84437). Original
probe-level matrices backed up as `*_probe_level.bak.tsv`. Automated by
`src/preprocessing/fix_geo_gene_mapping.py` (executed on A100, 2026-07-07).

### NK immune state classification

Four states defined from scRNA signatures; binary formulation
(NK-hot-cytotoxic vs rest) justified in Methods §2.4. Thresholds in
`nk_state_thresholds.json`.

### Composite target score weights

Five dimensions: tumor specificity 0.30, NK dysfunction 0.20, SST-axis
membership 0.30, axis core 0.10, literature 0.10. Rank stability under
±0.10 perturbation confirmed (see target rank stability analysis).

---

*Generated 2026-07-09. Verify against actual file listing before submission.*
