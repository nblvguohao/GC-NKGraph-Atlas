# GC-NKGraph-Atlas Codex Execution Document (v2)

> **v2 integration note.** This version reframes the project around the
> operationalization of a specific published immune-evasion mechanism
> (Zheng et al., Nat Immunol 2023, DOI 10.1038/s41590-023-01462-9:
> tumor serine metabolism → NK-cell sphingomyelin loss → collapse of membrane
> protrusions/microvilli → failure of lytic immune synapse → loss of
> cytotoxicity; rescued by inhibiting SM catabolism, synergistic with Tim3/HAVCR2
> blockade). The old gated "CROWN topology" module is replaced by an ACTIVE,
> computable transcriptional axis (Phase 14R) plus a still-GATED physical layer.
> A two-arm design (liver = positive control that recovers the mechanism where it
> was proven; gastric = novel extension) is now the spine of the paper. The
> pipeline is driven by a reusable **mechanism-card** abstraction so the same
> engine can operationalize other published mechanisms.
>
> **Data-availability reality check.** The anchor paper deposited NO transcriptomic
> data (data "available upon reasonable request"; source data are per-figure
> statistical XLSX; its single-cell data are single-immunocyte mass-spectrometry
> lipidomics, not scRNA-seq). Therefore the liver positive-control arm uses
> INDEPENDENT public liver transcriptomes (TCGA-LIHC + public HCC scRNA with an
> NK subset), never fabricated accessions. The transcriptional footprint of this
> mechanism is genuinely unexplored — that is the additive contribution.

## 0. Project Identity

Project name:

**GC-NKGraph-Atlas**

Full title:

**GC-NKGraph-Atlas: a single-cell atlas-guided heterogeneous graph learning framework for identifying NK-cell dysfunction and tumor immune-evasion programs in gastric cancer**

Manuscript target:

**Briefings in Bioinformatics / Cell Reports Methods / Genome Medicine-style computational method, benchmark, or problem-solving protocol.**

Core positioning:

This project merges the strengths of two previous directions:

1. **NK-GCGraph** contributes the main BIB-style method framework:
   - gastric cancer bulk cohorts;
   - single-cell and spatial validation;
   - NK immune-state labels;
   - tumor-NK heterogeneous graph;
   - benchmark, ablation, target prioritization, reproducibility tests.

2. **HCC-NK** contributes the biological evidence layer:
   - NK/ILC annotation discipline;
   - single-cell NK reference atlas;
   - scANVI/scArches reference mapping;
   - dysfunction / tissue-resident trajectory;
   - reversible-state proxy labels;
   - future transcriptome-to-membrane-topology interface.

The merged project is **not** a generic gastric cancer prognosis-signature project and is **not** a public single-cell atlas only. The central deliverable is a reproducible NK-aware multimodal graph-learning protocol for gastric cancer immune-evasion target discovery.

**v2 central positioning (mechanism operationalization).** The project's spine is
the computational reconstruction of the serine–sphingomyelin–membrane-topology
axis of NK immune evasion (Zheng et al., Nat Immunol 2023) from tumor
transcriptomes, with a two-arm design:

```text
Arm A — POSITIVE CONTROL (liver / HCC)
    Recover the published serine->SM->protrusion-machinery->cytotoxicity axis in
    the system where it was proven. Cohorts: TCGA-LIHC + public HCC scRNA with an
    NK subset (NO fabricated accessions; the anchor paper has no deposited omics).
    This is the credibility anchor.

Arm B — NOVEL EXTENSION (gastric cancer)
    Test whether the SAME axis operates in gastric cancer — a digestive-tract
    cancer NOT yet on the mechanism's published extension list (which already
    covers lung, colon, ovarian). Output gastric-specific tumor-intrinsic targets.
```

The engine is driven by a reusable **mechanism-card** (see 6.7 and
`configs/mechanism_cards/`): the framework is a machine that turns a published
wet-lab mechanism into a scalable computational target-discovery run, not a
one-off pipeline.

---

## 1. One-Sentence Scientific Objective

Use gastric cancer bulk transcriptomics, gastric cancer single-cell RNA-seq, optional spatial transcriptomics, ligand-receptor knowledge, PPI/pathway/TF-target prior networks, and a curated NK-cell reference atlas to define reproducible NK immune states, construct a tumor-NK heterogeneous graph, train a graph-learning model to predict NK dysfunction/exclusion programs, and prioritize tumor-intrinsic candidate targets that may drive NK-cell immune evasion.

---

## 2. Main Biological Question

The project should answer:

> Which tumor-intrinsic genes, pathways, and microenvironmental programs are associated with NK-cell dysfunction, NK exclusion, and impaired cytotoxicity in gastric cancer, and can a knowledge-guided graph model rank candidate targets for downstream experimental validation?

---

## 3. Final Project Shape

The merged project has four evidence layers.

```text
Layer 1: Bulk gastric cancer cohorts
    TCGA-STAD + GSE62254/ACRG + GSE84437 + optional cohorts
    → NK infiltration / cytotoxicity / dysfunction / exclusion scores
    → sample-level NK immune-state labels

Layer 2: Gastric cancer single-cell atlas
    GC scRNA-seq datasets
    → NK/ILC/T/NKT separation
    → NK subtypes and dysfunction trajectory
    → scRNA-derived marker/state evidence

Layer 3: Spatial validation
    GC spatial transcriptomics if available
    → NK niche / tumor adjacency / CAF-ECM exclusion / ligand-receptor context
    → if unavailable, create DATA_UNAVAILABLE.md and keep pipeline ready

Layer 4: Tumor-NK heterogeneous graph model
    bulk + scRNA + spatial + ligand-receptor + PPI + pathway + TF-target
    → NK-state prediction
    → dysfunction regression
    → candidate gene target scoring
    → edge/module importance
```

---

## 4. Hard Rules

1. Do not build a simple TCGA + LASSO Cox + nomogram project.
2. Do not use overall survival as the only label.
3. Do not use only TCGA-STAD.
4. Do not mix training and external validation data during normalization.
5. Do not compute thresholds, marker selection, feature selection, GWAS-like statistics, or signature thresholds on test or external validation data.
6. Do not claim biological discovery unless supported by at least:
   - bulk + single-cell evidence, or
   - bulk + spatial evidence, or
   - bulk + single-cell + literature/druggability evidence.
7. Do not silently drop failed datasets. Every failed dataset must be logged.
8. Do not fabricate GEO, ArrayExpress, CNGB, TISCH2, or other accession identifiers.
9. Do not fabricate topology values.
10. Any mock topology data must be explicitly named `*_MOCK_*` and must be excluded from scientific claims.
11. All scripts must be executable from command line.
12. All generated tables must include source dataset, processing date, script name, and config file.
13. Every model result must include random seed, split ID, hyperparameters, and git commit if available.
14. Every biologically ambiguous decision must be marked `NEEDS_REVIEW` in `results/logs/LOG.md`.
15. (v2) Never claim the pipeline predicts physical membrane topology, microvilli density, membrane roughness, immune-synapse area, or sphingomyelin metabolite content from transcriptome. The only permitted wording is "transcriptional program permissive-of / associated-with the topology phenotype." Every SST-axis output table and figure caption must carry this qualifier.
16. (v2) Transcription captures machinery/capacity only, not metabolite flux. Any serine/SM crosstalk statement must acknowledge that metabolite-level validation requires metabolomics / the origin lab's single-cell mass spectrometry.
17. (v2) Cell-type attribution is mandatory before any tumor-vs-NK axis claim. Results computed on cell-type-unresolved data are labeled `MIXED_UNRESOLVED` and may not support crosstalk claims.
18. (v2) The direction (sign) of the tumor-serine → NK-SM crosstalk term is CALIBRATED on the liver positive-control cohort, never hard-coded. Until calibrated, it is `NEEDS_REVIEW`, and every downstream `metabolic_crosstalk` edge and `sst_axis_score` inherits that status.
19. (v2) The liver positive-control arm uses independent public transcriptomes only. The anchor paper has no deposited omics; do not invent an accession for it.

---

## 5. Repository Structure

Create this structure if it does not already exist.

```text
GC-NKGraph-Atlas/
├── README.md
├── environment.yml
├── configs/
│   ├── data_config.yaml
│   ├── scrna_config.yaml
│   ├── graph_config.yaml
│   ├── model_config.yaml
│   ├── experiment_config.yaml
│   ├── topology_schema.yaml
│   ├── sst_axis_config.yaml
│   └── mechanism_cards/
│       ├── registry.yaml
│       └── zheng_nk_sm_topology.yaml
├── data/
│   ├── raw/
│   │   ├── bulk/
│   │   ├── scrna/
│   │   ├── spatial/
│   │   ├── prior_networks/
│   │   └── topology/
│   ├── interim/
│   │   ├── bulk/
│   │   ├── scrna/
│   │   └── spatial/
│   ├── processed/
│   │   ├── bulk/
│   │   ├── scrna/
│   │   ├── spatial/
│   │   ├── graph/
│   │   └── topology/
│   └── metadata/
├── src/
│   ├── common/
│   ├── data_download/
│   ├── preprocessing/
│   ├── immune_scoring/
│   ├── scrna_analysis/
│   ├── nk_atlas/
│   ├── trajectory/
│   ├── spatial_analysis/
│   ├── graph_construction/
│   ├── models/
│   ├── baselines/
│   ├── training/
│   ├── evaluation/
│   ├── interpretation/
│   └── topology/
├── notebooks/
├── results/
│   ├── tables/
│   ├── figures/
│   ├── checkpoints/
│   ├── models/
│   └── logs/
├── manuscript/
│   ├── figures/
│   ├── supplementary_tables/
│   └── notes/
└── tests/
```

Initialization command:

```bash
mkdir -p configs \
data/raw/{bulk,scrna,spatial,prior_networks,topology} \
data/interim/{bulk,scrna,spatial} \
data/processed/{bulk,scrna,spatial,graph,topology} \
data/metadata \
src/{common,data_download,preprocessing,immune_scoring,scrna_analysis,nk_atlas,trajectory,spatial_analysis,graph_construction,models,baselines,training,evaluation,interpretation,topology} \
notebooks results/{tables,figures,checkpoints,models,logs} \
manuscript/{figures,supplementary_tables,notes} tests
```

---

## 6. Config Design

### 6.1 `configs/data_config.yaml`

```yaml
project_name: GC-NKGraph-Atlas
primary_cancer: gastric_cancer
primary_bulk_dataset: TCGA-STAD
seed: 0

bulk_datasets:
  - name: TCGA-STAD
    role: train_primary
    expression_path: data/processed/bulk/tcga_stad_expression.tsv
    clinical_path: data/processed/bulk/tcga_stad_clinical.tsv
  - name: GSE62254
    role: external_validation
    expression_path: data/processed/bulk/gse62254_expression.tsv
    clinical_path: data/processed/bulk/gse62254_clinical.tsv
  - name: GSE84437
    role: external_validation
    expression_path: data/processed/bulk/gse84437_expression.tsv
    clinical_path: data/processed/bulk/gse84437_clinical.tsv

optional_bulk_datasets:
  - GSE15459
  - GSE26942

# (v2) Liver positive-control arm (Arm A). Independent public transcriptomes only.
positive_control_bulk_datasets:
  - name: TCGA-LIHC
    role: positive_control_liver
    cancer: liver_hcc
    purpose: recover_serine_SM_topology_axis_where_proven
    expression_path: data/processed/bulk/tcga_lihc_expression.tsv
    clinical_path: data/processed/bulk/tcga_lihc_clinical.tsv
  # optional additional public HCC cohorts may be added; never fabricate accessions.

# (v2) scRNA datasets. Verify contents/QC before use; log failures.
scrna_datasets:
  # Bridges gastric + liver + NK dysfunction via gastric-cancer liver metastasis.
  # Oncogene 2024. VERIFY accession, NK-cell count, and metadata before relying on it.
  - name: GSE246662
    role: candidate_gastric_and_liver_bridge
    status: VERIFY_REQUIRED
    note: "gastric cancer liver metastasis scRNA with impaired NK function"
  # Add additional GC and HCC scRNA datasets with an NK subset as found in Phase 1.

normalization:
  bulk_expression: log2_tpm_plus_1_or_platform_specific
  external_validation: normalize_independently
  batch_correction: false_by_default

gene_id:
  preferred: HGNC_symbol
  duplicate_strategy: keep_highest_mean_expression

logging:
  log_file: results/logs/LOG.md
  failed_dataset_table: results/tables/dataset_status.tsv
```

### 6.2 `configs/scrna_config.yaml`

```yaml
scrna_project:
  primary_cancer: gastric_cancer
  required_min_datasets: 1
  preferred_platforms:
    - 10x
    - snRNA-seq
  required_metadata_fields:
    - patient_id
    - tissue
    - dataset
    - platform
    - condition

qc:
  use_mad_filtering: true
  remove_doublets: true
  doublet_method: scrublet
  min_genes: null
  max_mito_percent: null

integration:
  default_method: scvi
  fallback_method: harmony
  batch_keys:
    - dataset
    - patient_id

annotation:
  celltypist: true
  marker_validation: true
  needs_review_log: true

nk_ilc_markers:
  keep:
    - NCAM1
    - KLRD1
    - NKG7
    - GNLY
    - KLRF1
    - EOMES
  exclude_t:
    - CD3D
    - CD3E
    - CD3G
    - TRAC
  tissue_resident_ilc1_like:
    - CXCR6
    - ITGA1
    - ITGAE
    - ZNF683
  circulating_nk:
    - S1PR5
    - FCGR3A
```

### 6.3 `configs/graph_config.yaml`

```yaml
node_types:
  - gene
  - nk_receptor
  - ligand
  - pathway
  - cell_state
  - sample
  - spatial_niche
  - malignant_program
  - topology_state

edge_types:
  - ppi
  - ligand_receptor
  - pathway_membership
  - tf_target
  - coexpression
  - cell_state_marker
  - spatial_adjacency
  - tumor_specific_expression
  - scrna_state_gene
  - dysfunction_correlation
  - exclusion_correlation
  - topology_association

prior_networks:
  required:
    - PPI
    - ligand_receptor
    - pathway_membership
    - TF_target
    - NK_marker_receptor_list
    - gastric_cancer_target_list

candidate_gene_pool:
  include:
    - nk_receptors_ligands
    - hypoxia_genes
    - glycolysis_lactate_genes
    - adenosine_pathway_genes
    - tgfb_pathway_genes
    - caf_ecm_genes
    - immune_checkpoint_genes
    - gastric_cancer_target_genes
    - malignant_cell_de_genes
    - nk_cytotoxicity_negative_correlated_genes
    - nk_dysfunction_positive_correlated_genes
    - spatial_colocalized_genes
    - trajectory_driver_genes
  max_size: 2000

seed_genes:
  metabolic_suppression:
    - LDHA
    - SLC16A3
    - SLC16A1
    - CA9
  adenosine_pathway:
    - NT5E
    - ENTPD1
    - ADORA2A
  nk_inhibitory_ligand:
    - PVR
    - NECTIN2
    - HLA-E
  stress_ligand_shedding:
    - MICA
    - MICB
    - ULBP1
    - ULBP2
    - ULBP3
    - ADAM10
    - ADAM17
  tgfb_ecm_exclusion:
    - TGFB1
    - TGFBR1
    - TGFBR2
    - COL1A1
    - COL1A2
    - FN1
    - FAP
  gastric_targets:
    - CLDN18
    - ERBB2
    - FGFR2
    - MET
```

### 6.4 `configs/model_config.yaml`

```yaml
model_name: GC-NKGraph-Atlas
task_type: multi_task

inputs:
  sample_features: true
  gene_features: true
  graph_edges: true
  edge_types: true
  scrna_state_features: true
  spatial_features_if_available: true
  topology_features_gated: true

outputs:
  nk_state_classification: true
  nk_dysfunction_regression: true
  nk_exclusion_regression: true
  candidate_target_score: true
  edge_importance: true
  state_gene_importance: true

architecture:
  preferred:
    - heterograph_transformer
    - edge_type_attention
  fallback:
    - HGT
    - HAN
    - GAT
    - GraphSAGE

training:
  seeds: [1234, 2345, 3456, 4567, 5678]
  epochs: 200
  patience: 20
  batch_size: 64
  learning_rate: 0.001
  weight_decay: 0.0001

loss:
  classification: focal_loss
  regression: huber
  ranking: pairwise_margin
  graph_regularization: edge_type_l2
```

### 6.5 `configs/experiment_config.yaml`

```yaml
splits:
  internal:
    dataset: TCGA-STAD
    strategy: stratified_kfold
    n_folds: 5
  external:
    datasets:
      - GSE62254
      - GSE84437

thresholding:
  train_only_quantiles: true
  save_thresholds: true
  external_uses_saved_thresholds: true

metrics:
  classification:
    - AUROC
    - AUPRC
    - MacroF1
    - BalancedAccuracy
    - MCC
  regression:
    - RMSE
    - MAE
    - Spearman
    - Pearson
  survival_if_available:
    - C_index
    - logrank_p
    - hazard_ratio
  target_ranking:
    - Precision_at_10
    - Precision_at_20
    - enrichment_odds_ratio
    - literature_support_rate
    - druggability_support_rate

ablation:
  - full_model
  - without_ligand_receptor_edges
  - without_pathway_edges
  - without_ppi_edges
  - without_spatial_edges
  - without_scrna_derived_labels
  - without_trajectory_features
  - without_candidate_prior_pool
  - bulk_only
  - graph_only
  - random_edges
```

### 6.6 `configs/topology_schema.yaml`

```yaml
status: GATED
description: >
  This schema is for future paired transcriptome-to-membrane-topology data.
  No real topology values are currently available. Do not use mock values for
  scientific claims.

sample_fields:
  - sample_id
  - cell_id
  - patient_id
  - cancer_type
  - tissue
  - assay_type

transcriptome_fields:
  - expression_matrix_path
  - gene_id_type
  - normalization_method

topology_targets:
  - membrane_protrusion_density
  - microvilli_density
  - membrane_roughness
  - immune_synapse_contact_area
  - topology_state_label

rules:
  allow_mock_for_code_test_only: true
  mock_file_pattern: "*_MOCK_*"
  real_training_requires_topology_targets: true
  missing_real_data_behavior: raise_clear_error
```

---

### 6.7 (v2) `configs/sst_axis_config.yaml` and `configs/mechanism_cards/`

The serine–sphingomyelin–topology axis (Phase 14R) is driven by
`configs/sst_axis_config.yaml`, and the whole pipeline is parameterized by a
reusable **mechanism card** under `configs/mechanism_cards/`. See Phase 14R for
the full `sst_axis_config.yaml`. The mechanism card encodes one published
mechanism as: origin/provenance, the computable transcriptional-proxy modules
(with cell-type attribution and expected directions), the GATED physical
ground-truth layer, graph-integration node/edge types, the positive-control +
extension validation plan, gold-standard genes, and the therapeutic hook /
wet-lab assays. The reference instance is
`configs/mechanism_cards/zheng_nk_sm_topology.yaml`.

```yaml
# configs/mechanism_cards/registry.yaml
active_mechanism_card: zheng_nk_sm_topology
cards_dir: configs/mechanism_cards
rules:
  - transcriptional_proxy != physical_ground_truth != metabolite_flux
  - crosstalk directions are calibrated on the positive-control cohort
  - physical ground-truth layer stays GATED unless real data are provided
```

---

## 7. Phase Plan

## Phase 0 — Scaffold, environment, logging

Tasks:

- [ ] Create repository structure.
- [ ] Create `.gitignore`.
- [ ] Create `environment.yml`.
- [ ] Add core packages:
  - scanpy
  - anndata
  - scvi-tools
  - scarches
  - squidpy
  - leidenalg
  - harmonypy
  - scrublet
  - celltypist
  - scikit-learn
  - shap
  - cellrank
  - torch
  - torch-geometric
  - xgboost
  - lightgbm
  - pandas
  - pyarrow
  - networkx
  - pyyaml
  - matplotlib
  - seaborn
- [ ] Add `src/common/logging.py`.
- [ ] Add `src/common/seed.py`.
- [ ] Add `src/common/io_utils.py`.
- [ ] Add import smoke test.

Acceptance:

```bash
python -c "import scanpy, scvi, squidpy, shap, torch, sklearn"
python -c "import yaml; yaml.safe_load(open('configs/data_config.yaml'))"
```

---

## Phase 1 — Dataset discovery, verification, and provenance

Tasks:

- [ ] Verify and download bulk gastric cancer cohorts:
  - TCGA-STAD;
  - GSE62254 / ACRG;
  - GSE84437;
  - optional GSE15459 / GSE26942.
- [ ] Verify and download at least one gastric cancer scRNA-seq dataset.
- [ ] Search for gastric cancer spatial transcriptomics datasets.
- [ ] If spatial data unavailable, create `results/logs/DATA_UNAVAILABLE_spatial.md`.
- [ ] Download prior networks:
  - PPI;
  - ligand-receptor;
  - pathway membership;
  - TF-target;
  - NK receptor / marker list;
  - gastric cancer target list.
- [ ] Create `data/metadata/dataset_manifest.tsv`.
- [ ] Create `results/tables/data_provenance.tsv`.
- [ ] Create `results/tables/dataset_status.tsv`.

Acceptance:

- every data source has accession / URL / date / script name;
- every failed source has `SKIPPED` reason;
- no invented accession.

---

## Phase 2 — Bulk preprocessing

Tasks:

- [ ] Standardize gene IDs to HGNC symbols.
- [ ] Resolve duplicates by highest mean expression.
- [ ] Normalize TCGA and external cohorts independently.
- [ ] Do not batch-correct external validation into training.
- [ ] Align clinical and expression samples.
- [ ] Export processed matrices.

Expected outputs:

```text
data/processed/bulk/tcga_stad_expression.tsv
data/processed/bulk/tcga_stad_clinical.tsv
data/processed/bulk/gse62254_expression.tsv
data/processed/bulk/gse62254_clinical.tsv
data/processed/bulk/gse84437_expression.tsv
data/processed/bulk/gse84437_clinical.tsv
results/tables/bulk_dataset_summary.tsv
```

---

## Phase 3 — Gastric cancer single-cell atlas and NK/ILC annotation

This phase imports the HCC-NK atlas discipline into the gastric cancer project.

Tasks:

- [ ] Read scRNA datasets into AnnData.
- [ ] Preserve raw counts in `.layers["counts"]`.
- [ ] Standardize metadata:
  - patient_id;
  - tissue;
  - dataset;
  - platform;
  - condition;
  - treatment if available.
- [ ] QC using MAD-based adaptive thresholds.
- [ ] Remove doublets with Scrublet.
- [ ] Integrate using scVI by default; Harmony fallback.
- [ ] Annotate major lineages with CellTypist + marker validation.
- [ ] Separate NK/ILC/T/NKT carefully.

NK/ILC/T rules:

```text
Keep NK/ILC candidates:
NCAM1, KLRD1, NKG7, GNLY, KLRF1, EOMES

Exclude or separately label T/NKT:
CD3D, CD3E, CD3G, TRAC

Tissue-resident / ILC1-like features:
CXCR6, ITGA1, ITGAE, ZNF683

Circulating NK features:
S1PR5, FCGR3A
```

Required output:

```text
data/processed/scrna/gc_integrated.h5ad
data/processed/scrna/gc_nk_subset.h5ad
results/tables/gc_scrna_dataset_summary.tsv
results/tables/nk_annotation_evidence.tsv
results/figures/scrna_umap_major_lineages.pdf
results/figures/scrna_umap_nk_subtypes.pdf
```

Acceptance:

- NK and T cells are not merged silently.
- NK/ILC1/NKT boundary decisions are logged as `NEEDS_REVIEW`.
- Annotation evidence table exists.

---

## Phase 4 — NK reference atlas and query mapping

Tasks:

- [ ] Reintegrate NK/ILC subset.
- [ ] Cluster NK/ILC cells.
- [ ] Annotate NK states:
  - CD56bright-like;
  - CD56dim-like;
  - cytotoxic NK;
  - tissue-resident NK;
  - ILC1-like;
  - proliferating NK;
  - dysfunctional NK.
- [ ] Train scANVI or scArches reference model.
- [ ] Save query-mapping model to `results/models/nk_reference/`.
- [ ] Test query mapping using held-out cells or held-out dataset.
- [ ] Score cells for:
  - cytotoxicity;
  - dysfunction/exhaustion;
  - tissue residency;
  - IL-15-response / reversible-like state if configured.

Expected outputs:

```text
results/models/nk_reference/
data/processed/scrna/gc_nk_reference.h5ad
results/tables/nk_subtype_markers.tsv
results/tables/nk_state_scores_single_cell.tsv
results/figures/fig2_nk_state_atlas.pdf
```

Acceptance:

- reference model can be saved and reloaded;
- query mapping works on held-out cells;
- state scores exist in AnnData.

---

## Phase 5 — NK immune-state definitions in bulk

Create `src/immune_scoring/nk_scores.py`.

Core signatures:

```python
NK_MARKERS = [
    "NCAM1", "NCR1", "KLRD1", "KLRK1", "KLRF1",
    "NKG7", "GNLY", "GZMB", "PRF1", "FCGR3A",
    "XCL1", "XCL2", "CCL5", "IFNG", "TYROBP"
]

NK_CYTOTOXICITY_GENES = [
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "XCL1", "XCL2", "CCL5"
]

NK_DYSFUNCTION_GENES = [
    "KLRC1", "TIGIT", "CD96", "HAVCR2", "TOX",
    "ENTPD1", "TGFB1", "HIF1A", "NT5E", "ADORA2A"
]

CAF_ECM_TGFB_GENES = [
    "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN",
    "ACTA2", "FAP", "TAGLN", "TGFBI", "TGFB1", "CXCL12"
]
```

Definitions:

```text
NK_infiltration_score = mean_zscore(NK_MARKERS)

NK_cytotoxicity_score = mean_zscore(NK_CYTOTOXICITY_GENES)

NK_dysfunction_score =
mean_zscore(NK_DYSFUNCTION_GENES) - mean_zscore(NK_CYTOTOXICITY_GENES)

NK_exclusion_score =
mean_zscore(CAF_ECM_TGFB_GENES) - NK_infiltration_score
```

Four deterministic immune states:

```text
NK-hot-cytotoxic:
  NK infiltration high
  NK cytotoxicity high
  NK dysfunction low/intermediate

NK-hot-dysfunctional:
  NK infiltration high
  NK cytotoxicity low/intermediate
  NK dysfunction high

NK-cold/excluded:
  NK infiltration low
  NK exclusion high

NK-intermediate:
  remaining samples
```

Rules:

- thresholds must be fitted only on training folds;
- thresholds must be saved;
- external validation must use saved thresholds;
- all label assumptions must be written to `results/tables/label_definition.md`.

Expected outputs:

```text
results/tables/nk_scores_bulk.tsv
results/tables/nk_state_labels.tsv
results/tables/nk_state_thresholds.json
results/tables/label_definition.md
```

---

## Phase 6 — Dysfunction trajectory and reversible-state proxy

Tasks:

- [ ] Run diffusion pseudotime / CellRank on GC NK atlas.
- [ ] Infer trajectory:
  - circulating/cytotoxic NK → tissue-resident NK → dysfunctional/excluded-like NK.
- [ ] Identify branch points.
- [ ] Identify trajectory-associated driver genes.
- [ ] Include cytoskeleton / membrane remodeling genes:
  - ACTB;
  - MSN / EZR / RDX;
  - WAS / WASF2;
  - Arp2/3 complex genes;
  - BAR-domain proteins;
  - mTOR-related genes.
- [ ] Define reversible-like versus terminal-like proxy label:
  - default: IL-15 response / rescue-associated signature;
  - alternative: branch-before versus branch-after;
  - alternative: dysfunction score regression target.
- [ ] Export label documentation and limitations.

Expected outputs:

```text
results/tables/trajectory_drivers.tsv
results/tables/reversible_proxy_label_definition.md
results/figures/nk_dysfunction_trajectory.pdf
```

Acceptance:

- no claim of true reversibility without experimental data;
- proxy label assumptions are explicit;
- trajectory genes feed into candidate gene pool.

---

## Phase 7 — Spatial validation arm

Tasks:

- [ ] Search and download gastric cancer spatial transcriptomics if available.
- [ ] If unavailable, create `results/logs/DATA_UNAVAILABLE_spatial.md`.
- [ ] Deconvolve or score NK/cytotoxic/dysfunctional/exclusion programs.
- [ ] Run spatial neighborhood analysis with Squidpy.
- [ ] Quantify tumor-NK adjacency, CAF-ECM exclusion, ligand-receptor niches.
- [ ] Generate spatial niche nodes for graph construction if spatial data exists.

Expected outputs if available:

```text
data/processed/spatial/spatial_metadata.tsv
data/processed/spatial/spatial_expression.tsv
data/processed/spatial/spatial_coordinates.tsv
data/processed/spatial/spatial_niche_scores.tsv
results/figures/fig7_spatial_validation.pdf
```

---

## Phase 8 — Heterogeneous graph construction

Create `src/graph_construction/build_heterograph.py`.

Node schema:

```text
node_id
node_name
node_type
source
description
```

Allowed node types:

```text
gene
nk_receptor
ligand
pathway
cell_state
sample
spatial_niche
malignant_program
topology_state
```

Edge schema:

```text
source_node_id
target_node_id
edge_type
weight
confidence
database_source
direction
evidence
```

Allowed edge types:

```text
ppi
ligand_receptor
pathway_membership
tf_target
coexpression
cell_state_marker
spatial_adjacency
tumor_specific_expression
scrna_state_gene
dysfunction_correlation
exclusion_correlation
topology_association
```

Expected outputs:

```text
data/processed/graph/nodes.tsv
data/processed/graph/edges.tsv
data/processed/graph/node_features.parquet
data/processed/graph/edge_features.parquet
data/processed/graph/graph_statistics.tsv
```

Important:

- no test-label-derived edges;
- no edges calculated from external validation labels;
- graph construction must be reproducible from config.

---

## Phase 9 — Baselines

Required baselines:

```text
ssGSEA NK score
xCell or local approximation
MCP-counter or local approximation
ElasticNet
RandomForest
XGBoost
LightGBM
SVM
MLP
GCN
GAT
GraphSAGE
HGT or HAN
```

Each baseline must write:

```text
results/tables/baseline_<method>_internal.tsv
results/tables/baseline_<method>_external.tsv
```

Required columns:

```text
method
dataset
split_id
seed
label_type
AUROC
AUPRC
MacroF1
BalancedAccuracy
MCC
hyperparameters
timestamp
script_name
config_file
```

If a package cannot be installed, implement a local approximation and document it in:

```text
results/logs/baseline_implementation_notes.md
```

---

## Phase 10 — GC-NKGraph model

Create:

```text
src/models/gc_nkgraph_atlas.py
src/training/train_gc_nkgraph.py
src/evaluation/evaluate_gc_nkgraph.py
src/interpretation/explain_gc_nkgraph.py
```

Model must support:

```text
heterogeneous node types
heterogeneous edge types
edge-type attention or edge-type embeddings
sample-level NK state prediction
NK dysfunction score regression
NK exclusion score regression
gene-level target scoring
edge importance extraction
state-gene importance extraction
```

Minimum acceptable fallback:

1. HGT/HAN backbone.
2. Edge-type embedding.
3. Sample-level multi-task heads.
4. Candidate-gene scoring MLP.
5. Edge importance extraction.
6. Ablations preserved.

Training outputs:

```text
results/checkpoints/gc_nkgraph_seed<SEED>_fold<FOLD>.pt
results/tables/gc_nkgraph_internal_results.tsv
results/tables/gc_nkgraph_external_results.tsv
results/tables/gc_nkgraph_predictions.tsv
results/tables/gc_nkgraph_gene_importance.tsv
results/tables/gc_nkgraph_edge_importance.tsv
```

---

## Phase 11 — Evaluation and leakage tests

Metrics:

```text
AUROC
AUPRC
MacroF1
BalancedAccuracy
MCC
RMSE
MAE
C-index if survival data available
log-rank P if survival data available
Precision@10 for candidate targets
Precision@20 for candidate targets
pathway enrichment odds ratio
```

Rules:

1. AUPRC and Macro-F1 are more important than AUROC for imbalanced state labels.
2. Report mean ± standard deviation across seeds/folds.
3. External validation must be reported separately.
4. Do not tune external validation thresholds.
5. Do not hide failed seeds or failed datasets.
6. Save failed runs with reason.

---

## Phase 12 — Ablation and robustness

Create `src/training/run_ablation.py`.

Run:

```text
full_model
without_ligand_receptor_edges
without_pathway_edges
without_ppi_edges
without_spatial_edges
without_scrna_derived_labels
without_trajectory_features
without_candidate_prior_pool
bulk_only
graph_only
random_edges
```

Outputs:

```text
results/tables/ablation_<variant>.tsv
results/tables/ablation_summary.tsv
results/figures/fig5_ablation.pdf
```

For each ablation, write one sentence explaining which biological or computational problem the removed module addresses.

---

## Phase 13 — Candidate target prioritization

Create `src/interpretation/prioritize_targets.py`.

Outputs:

```text
results/tables/top_candidate_targets.tsv
results/tables/candidate_evidence_matrix.tsv
results/figures/fig6_candidate_targets.pdf
```

Candidate evidence matrix columns:

```text
gene
target_score
rank
tumor_cell_specificity
nk_cytotoxicity_correlation
nk_dysfunction_correlation
nk_exclusion_correlation
trajectory_driver_score
spatial_colocalization_score
ligand_receptor_axis
pathway_membership
ppi_degree
edge_attention_score
scrna_state_support
literature_category
druggability_category
recommended_validation_assay
```

Candidate categories:

```text
metabolic_suppression
adenosine_pathway
nk_inhibitory_ligand
stress_ligand_shedding
caf_ecm_exclusion
gastric_cancer_target
trajectory_driver
unknown_candidate
```

Recommended validation assay examples:

```text
qPCR in malignant cells
IHC / multiplex IF
co-culture with NK cells
NK cytotoxicity assay
ligand blocking assay
CRISPR knockdown / overexpression
flow cytometry for ligand expression
```

---

## Phase 14R — Serine–Sphingomyelin–Topology transcriptional axis (SST-Axis)

> (v2) Replaces the old gated "CROWN topology interface." Converts the topology
> idea from a gated promise into a real, computable, mechanism-grounded module,
> while keeping the physical-imaging/MS layer honestly gated. Anchor: Zheng et al.
> Nat Immunol 2023 (DOI 10.1038/s41590-023-01462-9) + Clin Transl Med 2023
> (DOI 10.1002/ctm2.1395). Full config: `configs/sst_axis_config.yaml`; card:
> `configs/mechanism_cards/zheng_nk_sm_topology.yaml`.
>
> Execution order: runs AFTER Phase 3–6; feeds Phase 8 (graph), 12 (ablation),
> 13 (targets).

### 14R.0 Design split (read first)

```text
Layer 14R-A  — TRANSCRIPTIONAL AXIS  (status: ACTIVE / computable)
    The molecular MACHINERY / CAPACITY for serine->SM->protrusion->cytotoxicity,
    computable from bulk + scRNA. This is the real scientific contribution.

Layer 14R-B  — PHYSICAL GROUND TRUTH  (status: GATED)
    Protrusion/microvilli density, membrane roughness, synapse contact area (SEM/
    super-resolution) and single-cell MS sphingomyelin. No real values available.
    Loader raises a clear error if absent. Mock only for code tests, named
    *_MOCK_*, never in claims.
```

Honesty rule (whole phase): transcriptional proxy ≠ physical topology ≠ metabolite
flux. Only permitted claim: a transcriptional program **permissive-of /
associated-with** the phenotype.

### 14R.1 Gene modules (starting sets; must pass marker_validation)

```text
tumor_serine_capacity  (attribute to MALIGNANT cells):
  PHGDH, PSAT1, PSPH, SHMT1, SHMT2, MTHFD1, MTHFD2, MTHFD1L, SLC1A4, SLC1A5
  direction w.r.t. NK SM availability = NEEDS_REVIEW (calibrate in 14R.5)

nk_sm_synthesis  (NK):  SGMS1, SGMS2                      (higher = more permissive)
nk_sm_catabolism (NK):  SMPD1, SMPD2, SMPD3, SMPD4        (higher = less permissive)
nk_denovo_sphingolipid (NK): SPTLC1, SPTLC2, SPTLC3, SPTSSA, CERS2/4/5/6, DEGS1
nk_protrusion_machinery (NK):
  EZR, MSN, RDX, ACTR2, ACTR3, ARPC1B, ARPC2, ARPC3, ARPC4, ARPC5,
  WAS, WASL, WASF1, WASF2, WASF3, WIPF1, CDC42, RAC1, RHOA,
  DIAPH1, DIAPH3, FMNL1, BAIAP2, PACSIN2                  (higher = more permissive)
nk_synapse_cytotoxicity_outcome (NK):
  NKG7, GNLY, GZMB, PRF1, IFNG, LCP2, LAT, VAV1, TLN1, ITGAL, ITGB2
checkpoint_link (NK): HAVCR2                              (higher = less permissive)
```

Derived scores:

```text
tumor_serine_capacity_score   = mean_zscore(tumor_serine_capacity | malignant)
nk_sm_balance_score           = mean_zscore(nk_sm_synthesis) - mean_zscore(nk_sm_catabolism)
nk_protrusion_machinery_score = mean_zscore(nk_protrusion_machinery | nk)
nk_topology_permissive_score  = composite(nk_sm_balance_score, nk_protrusion_machinery_score)
sst_axis_score                = integrated(tumor_serine_capacity_score[calibrated sign],
                                           nk_topology_permissive_score, outcome)
```

### 14R.2 Tasks — Layer A (ACTIVE)

- [ ] Create `src/topology/sst_axis.py`; load modules from config; run marker_validation; log curation `NEEDS_REVIEW`.
- [ ] Attribute modules to cell types (tumor_serine → malignant; NK modules → NK subset) using Phase 3–4 outputs.
- [ ] Compute per-cell (scRNA) and per-sample (bulk, cell-type-adjusted) scores.
- [ ] scRNA test of the Zheng phenotype at molecular level: compare intratumoral vs peritumoral/peripheral NK for `nk_sm_balance_score` and `nk_protrusion_machinery_score`.
- [ ] Fit thresholds/standardization on training folds only; save; reuse for external cohorts.

Outputs:

```text
data/processed/scrna/gc_nk_sst_scores_single_cell.tsv
results/tables/sst_axis_scores_bulk.tsv
results/tables/sst_axis_module_scores.tsv
results/tables/nk_intratumoral_vs_peritumoral_sst.tsv
results/tables/label_definition_sst_axis.md
results/figures/fig9_sst_axis.pdf
```

### 14R.3 Tasks — Layer B (GATED)

- [ ] Keep `configs/topology_schema.yaml`; `src/topology/topology_infer.py` real loader raises clear error if physical/MS targets absent.
- [ ] `*_MOCK_*` dataset only for code-path tests. No physical/MS value in any scientific output.

Outputs: `configs/topology_schema.yaml`, `src/topology/topology_infer.py`, `results/logs/TOPOLOGY_GATED_STATUS.md`, `tests/test_topology_mock_pipeline.py`.

### 14R.4 Graph edges (feeds Phase 8)

```text
metabolic_crosstalk : tumor_serine_program node -> nk_topology_state node
    weight = |calibrated association| ; sign = calibrated (14R.5)
    requires cell-type-resolved endpoints; never from external-validation labels
sm_topology_axis    : within-axis gene-gene edges (coexpression in NK subset, train folds only)
```

### 14R.5 Positive-control protocol (Arm A — LIVER)

Pre-register in `manuscript/notes/sst_axis_prereg.md` BEFORE touching test data:

```text
H1 tumor_serine_capacity  ⟂  nk_sm_balance      (sign CALIBRATED and reported)
H2 nk_sm_balance          (+) nk_protrusion_machinery
H3 nk_protrusion_machinery(+) cytotoxicity anchors
H4 nk_topology_permissive (−) HAVCR2 / dysfunction
H5 intratumoral NK < peritumoral/peripheral NK in sm_balance & protrusion machinery (scRNA)
```

Data: TCGA-LIHC (+ optional public HCC cohorts, no fabricated accessions);
public HCC scRNA with NK subset; else `results/logs/DATA_UNAVAILABLE_lihc_nk.md`.

Recovery = H2–H5 pass in the pre-registered direction in liver. The calibrated H1
sign is the sign used for all downstream crosstalk edges and `sst_axis_score`.

Outputs:

```text
manuscript/notes/sst_axis_prereg.md
results/tables/sst_axis_positive_control_liver.tsv
results/figures/fig10_positive_control_liver.pdf
```

### 14R.6 Extension (Arm B — GASTRIC)

- [ ] Repeat 14R.5 hypotheses in gastric cohorts; report whether the axis operates in GC (the novel result).
- [ ] Feed gastric SST-axis genes into the candidate pool (Phase 13) with `sst_axis_membership`.

Outputs: `results/tables/sst_axis_gastric.tsv`, `results/figures/fig11_sst_axis_gastric.pdf`.

### 14R.7 Candidate hooks (feeds Phase 13)

Add columns to `candidate_evidence_matrix.tsv`: `sst_axis_membership`,
`sst_axis_direction_consistent`, `positive_control_recovered`. Add wet-lab assays
to the assay vocabulary: NK–tumor co-culture cytotoxicity; SEM/super-resolution
protrusion imaging; single-immunocyte MS (membrane SM); lytic synapse formation;
sphingomyelinase-inhibitor rescue ± Tim3/HAVCR2 blockade.

### 14R.8 In-silico perturbation / stratification (optional; feeds Phase 10 & 13)

```text
Per sample, estimate predicted change in cytotoxicity / topology_permissive under
an in-silico "SM-restoration" perturbation (raise nk_sm_synthesis / suppress
nk_sm_catabolism). Rank samples by predicted benefit from
"SM-restoration + Tim3 blockade" logic.
Outputs: results/tables/sst_axis_insilico_perturbation.tsv,
         results/tables/sst_axis_combination_stratification.tsv
Honesty: model-based hypothesis, NOT a validated predictor. Label accordingly.
```

### 14R.9 Tests (extend Phase 17)

```text
tests/test_sst_axis_celltype_attribution.py
tests/test_sst_axis_leakage.py
tests/test_sst_axis_gating.py
```

### 14R.10 Acceptance (add to Section 18)

```text
[ ] SST modules validated and scored with cell-type attribution.
[ ] scRNA intratumoral-vs-peritumoral NK comparison generated.
[ ] Liver positive-control recovery reported PASS/FAIL vs pre-registration.
[ ] Tumor-serine crosstalk sign CALIBRATED (not assumed) and reused downstream.
[ ] Gastric extension reported as the novel result.
[ ] Physical topology + MS-SM remain GATED; no mock values in claims.
[ ] Every SST output carries the "transcriptional proxy, not physical topology" qualifier.
```

---

## Phase 15 — Figure generation

Create one script per figure under `src/interpretation/figures/`.

Required figures:

```text
Figure 1: integrated workflow and model overview
Figure 2: gastric cancer NK single-cell atlas
Figure 3: tumor-NK heterogeneous graph and model architecture
Figure 4: benchmark comparison
Figure 5: ablation and robustness analysis
Figure 6: candidate target prioritization
Figure 7: single-cell and spatial validation
Figure 8: (v2) SST-axis schematic + mechanism card (serine->SM->topology->cytotoxicity)
Figure 9: (v2) SST-axis scores in scRNA and bulk (intratumoral vs peritumoral NK)
Figure 10: (v2) liver positive-control recovery of the axis
Figure 11: (v2) gastric-cancer extension of the axis
```

Expected outputs:

```text
results/figures/fig1_workflow.pdf
results/figures/fig2_nk_state_atlas.pdf
results/figures/fig3_model_architecture.pdf
results/figures/fig4_benchmark.pdf
results/figures/fig5_ablation.pdf
results/figures/fig6_candidate_targets.pdf
results/figures/fig7_validation.pdf
results/figures/fig8_sst_axis_schematic.pdf
results/figures/fig9_sst_axis.pdf
results/figures/fig10_positive_control_liver.pdf
results/figures/fig11_sst_axis_gastric.pdf
```

---

## Phase 16 — Manuscript-oriented deliverables

Create:

```text
manuscript/notes/main_claims.md
manuscript/notes/method_summary.md
manuscript/notes/result_summary.md
manuscript/notes/limitations.md
manuscript/supplementary_tables/table_s1_dataset_manifest.tsv
manuscript/supplementary_tables/table_s2_marker_genes.tsv
manuscript/supplementary_tables/table_s3_graph_statistics.tsv
manuscript/supplementary_tables/table_s4_baseline_results.tsv
manuscript/supplementary_tables/table_s5_ablation_results.tsv
manuscript/supplementary_tables/table_s6_candidate_targets.tsv
manuscript/supplementary_tables/table_s7_topology_schema.tsv
```

`main_claims.md` must use only generated evidence.

Format:

```text
Claim:
Evidence:
Related figure/table:
Limitations:
```

---

## Phase 17 — Tests

Create:

```text
tests/test_data_integrity.py
tests/test_no_leakage.py
tests/test_graph_edges.py
tests/test_reproducibility.py
tests/test_scrna_annotation_outputs.py
tests/test_topology_mock_pipeline.py
tests/test_sst_axis_celltype_attribution.py
tests/test_sst_axis_leakage.py
tests/test_sst_axis_gating.py
```

Required checks:

`test_data_integrity.py`

```text
No duplicated sample IDs.
No missing labels in training data.
Gene symbols are unique after preprocessing.
Clinical rows match expression columns.
```

`test_no_leakage.py`

```text
External validation samples are never used for training.
External validation samples are not used to compute thresholds.
Feature selection is performed only on training folds.
Normalization parameters are not transferred from external validation to training.
scRNA-derived markers used in bulk model are generated without external validation labels.
```

`test_graph_edges.py`

```text
All edge endpoints exist in nodes.tsv.
All edge types are recognized.
All edge weights are finite.
No test-label-derived edges exist.
```

`test_reproducibility.py`

```text
Same seed reproduces same split.
Same seed reproduces same model metrics within tolerance.
All config files are saved with results.
```

`test_scrna_annotation_outputs.py`

```text
NK subset exists.
NK/T separation evidence table exists.
Required marker columns exist.
AnnData can be read back.
```

`test_topology_mock_pipeline.py`

```text
Mock topology files are explicitly named *_MOCK_*.
Mock topology pipeline runs end-to-end.
Real topology loader raises error if topology targets are absent.
Mock topology outputs are not included in candidate target claims.
```

---

## 18. Final Acceptance Criteria

The project is acceptable only if all items below are complete.

```text
[ ] At least 3 bulk gastric cancer cohorts included.
[ ] At least 1 gastric cancer single-cell dataset included.
[ ] At least 1 spatial validation dataset included or documented as unavailable.
[ ] NK/ILC/T annotation evidence table generated.
[ ] NK reference atlas generated.
[ ] Query mapping model saved and reloadable.
[ ] NK immune-state labels generated reproducibly.
[ ] Training-only thresholds saved and reused for external cohorts.
[ ] Dysfunction trajectory and trajectory-driver table generated.
[ ] Heterogeneous graph created with typed nodes and typed edges.
[ ] At least 8 baselines implemented.
[ ] GC-NKGraph trained across at least 5 seeds or 5 folds.
[ ] External validation reported separately.
[ ] Ablation experiments completed.
[ ] Top candidate target table generated.
[ ] Candidate evidence matrix generated.
[ ] At least 6 manuscript figures generated.
[ ] No data leakage tests pass.
[ ] Topology interface is clearly marked GATED and does not fabricate true values.
[ ] (v2) SST-axis modules scored with mandatory cell-type attribution.
[ ] (v2) Liver positive-control recovery of the axis reported PASS/FAIL vs pre-registration.
[ ] (v2) Tumor-serine crosstalk sign CALIBRATED on liver, not assumed, and reused downstream.
[ ] (v2) Gastric-cancer extension of the axis reported as the novel result.
[ ] (v2) Every SST output carries the "transcriptional proxy, not physical topology" qualifier.
[ ] (v2) Active mechanism card present and consumed from configs/mechanism_cards/.
[ ] README can reproduce the pipeline.
```

---

## 19. MVP Priority

If time is limited, execute in this order:

1. Repository structure, configs, and active mechanism card.
2. TCGA-STAD + GSE62254 + GSE84437 preprocessing; TCGA-LIHC (liver control) preprocessing.
3. Bulk NK scores and four immune-state labels.
4. One gastric cancer scRNA dataset annotation (verify GSE246662 or alternative).
5. NK/ILC/T separation and NK subtype validation.
6. scRNA-derived NK marker/state table.
7. (v2) SST-axis Layer A: cell-type-attributed serine/SM/protrusion module scores (scRNA + bulk).
8. (v2) Liver positive control (14R.5): pre-register + test H1–H5; calibrate crosstalk sign.
9. Ligand-receptor + PPI + pathway graph, incl. metabolic_crosstalk + sm_topology_axis edges.
10. XGBoost, MLP, GAT, HGT/HAN baselines.
11. GC-NKGraph full model.
12. Ablations:
    - without LR;
    - without pathway;
    - without scRNA labels;
    - without SST-axis edges;
    - bulk-only;
    - random-edges.
13. (v2) Gastric extension of the axis (14R.6) + top 20 candidate targets with sst_axis_membership.
14. Figures 1–11.
15. README and tests.
16. Spatial validation if data is available.
17. Physical topology (Layer B) remains a GATED extension.

Do not make any physical-topology or wet-lab claim, and do not treat the crosstalk
sign as known, before the liver positive control (step 8) is complete.

---

## 20. README Command Skeleton

```bash
conda env create -f environment.yml
conda activate gc-nkgraph-atlas

python src/preprocessing/run_bulk_preprocessing.py --config configs/data_config.yaml
python src/scrna_analysis/run_scrna_pipeline.py --config configs/scrna_config.yaml
python src/nk_atlas/build_nk_reference.py --config configs/scrna_config.yaml
python src/immune_scoring/nk_scores.py --config configs/data_config.yaml
python src/trajectory/run_nk_trajectory.py --config configs/scrna_config.yaml
python src/topology/sst_axis.py --config configs/sst_axis_config.yaml --card configs/mechanism_cards/zheng_nk_sm_topology.yaml
python src/spatial_analysis/run_spatial_validation.py --config configs/data_config.yaml
python src/graph_construction/build_heterograph.py --config configs/graph_config.yaml
python src/baselines/run_all_baselines.py --config configs/experiment_config.yaml
python src/training/train_gc_nkgraph.py --config configs/model_config.yaml
python src/training/run_ablation.py --config configs/experiment_config.yaml
python src/interpretation/prioritize_targets.py --config configs/model_config.yaml
python src/interpretation/generate_all_figures.py
pytest tests/
```

---

## 21. Manuscript Framing (v2)

Recommended title:

**GC-NKGraph-Atlas: reconstructing the serine–sphingomyelin–membrane-topology axis of NK-cell immune evasion from tumor transcriptomes, from liver to gastric cancer**

Alternatives:

- **A single-cell-informed tumor–NK heterogeneous graph reconstructs a physical immune-evasion mechanism and prioritizes gastric-cancer targets**
- **Operationalizing NK surface-topology evasion: a knowledge-guided graph framework recovers a metabolic-topology axis in liver and extends it to gastric cancer**

Core claims (conservative, mechanism-grounded):

1. We define a reproducible, cell-type-attributed transcriptional proxy for the serine → sphingomyelin → membrane-protrusion-machinery → cytotoxicity axis of NK immune evasion.
2. In liver cancer (the mechanism's origin system), the framework **recovers** the axis: pre-registered relationships replicate in the expected direction from independent public cohorts (the anchor paper deposited no omics, so this is the first transcriptional test of the mechanism).
3. We construct a single-cell-informed tumor–NK heterogeneous graph whose `metabolic_crosstalk` edge is justified by the biology, giving the heterogeneous-graph architecture a principled reason to exist.
4. The graph model improves or stabilizes NK-state prediction versus strong tabular and graph baselines under internal and external validation.
5. The same axis is tested in gastric cancer; where it holds, we prioritize tumor-intrinsic candidate targets with multi-evidence support and a recommended wet-lab assay.
6. A model-based in-silico "SM-restoration" readout stratifies samples by predicted benefit from an SM-restoration + Tim3-blockade logic — a hypothesis for testing, not a validated predictor.

What is a control vs what is novel (state explicitly):

```text
Positive control (essential, not novel): recovery of the axis in liver.
Novel: scalable transcriptional reconstruction; the mechanism-grounded graph
       (metabolic_crosstalk edge); the gastric extension + target list; the
       transferable mechanism-card engine.
```

Avoid overclaiming:

- Do not claim prediction of physical membrane topology, microvilli density, or SM metabolite content from transcriptome; say "transcriptional program permissive-of / associated-with the topology phenotype."
- Do not present transcription as a substitute for metabolite-level crosstalk (needs metabolomics / single-cell MS).
- Do not hard-code the tumor-serine → NK-SM crosstalk sign; report the sign calibrated on liver.
- Do not say targets are experimentally validated.
- Do not say true NK reversibility is proven.
- Do not say graph attention is causal evidence.
- Do not present survival association as mechanism.
- Do not make tumor-vs-NK axis claims on cell-type-unresolved data.

---

## 22. Codex Execution Summary

Execute this as an engineering pipeline, not as a literature essay.

The final output should let a human run:

```bash
pytest tests/
python src/interpretation/generate_all_figures.py
```

and obtain the main tables and figures for a BIB-level manuscript.

Central deliverable:

```text
A reproducible, single-cell-informed, NK-aware multimodal graph learning protocol
that reconstructs the serine->sphingomyelin->membrane-topology axis of NK immune
evasion from tumor transcriptomes, recovers it in liver (positive control),
extends it to gastric cancer (novel), and prioritizes tumor-intrinsic
immune-evasion targets with wet-lab-ready assays — driven by a reusable
mechanism-card so the same engine can operationalize other published mechanisms.
```
