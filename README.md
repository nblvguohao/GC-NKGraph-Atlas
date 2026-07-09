# GC-NKGraph-Atlas

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-manuscript%20submitted-lightgrey.svg)]()

**Reconstructing the serine–sphingomyelin–membrane-topology axis of NK-cell immune
evasion from tumor transcriptomes — a single-cell-informed heterogeneous graph
framework, from liver to gastric cancer.**

This repository accompanies the manuscript submitted to *Briefings in
Bioinformatics*. It operationalizes a published wet-lab immune-evasion mechanism
(Zheng et al., *Nat Immunol* 2023) as a reusable **mechanism-card**, reconstructs
the axis from public transcriptomes under a two-arm design (liver positive
control + gastric-cancer extension), benchmarks a mechanism-grounded heterogeneous
graph model against tabular baselines, and prioritizes tumor-intrinsic candidate
targets.

> **Honest scope.** The framework produces a *map of the mechanism's
> transcriptional reach*, not a blanket recovery: the effector layer
> (protrusion→cytotoxicity) recovers and generalizes; the metabolic coupling
> recovers only at single-cell resolution; the physical topology phenotype is not
> captured by machinery transcription. See the manuscript for the full framing.

---

## Table of Contents

- [Overview](#overview)
- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Quick start (synthetic data)](#quick-start-synthetic-data)
- [Data](#data)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Mechanism cards](#mechanism-cards)
- [Key result files](#key-result-files)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

### The biological question

A landmark wet-lab study (Zheng et al., *Nature Immunology* 2023) showed that
tumors evade NK-cell cytotoxicity through a physical mechanism: dysregulated tumor
serine metabolism depletes sphingomyelin in NK membranes, collapsing the membrane
protrusions required for lytic immune synapses. This was proven with single-cell
mass spectrometry and super-resolution imaging — techniques that cannot scale to
cohort-level analyses.

**This project asks:** how much of the serine→sphingomyelin→membrane-topology→
cytotoxicity axis can be reconstructed from widely available public transcriptomes?
And where does transcription stop being a valid proxy?

### What we built

GC-NKGraph-Atlas is a computational framework with three innovations:

1. **Mechanism-card abstraction** — a machine-readable YAML specification that
   encodes a published wet-lab mechanism as a recipe the pipeline consumes. New
   mechanisms = new cards, not new code.
2. **Mechanism-grounded heterogeneous graph** — a tumor–NK gene graph where every
   edge type carries an explicit biological justification, including a
   `metabolic_crosstalk` edge derived directly from the Zheng 2023 mechanism.
3. **Honest scoping** — a two-arm design (liver positive control + gastric cancer
   extension) that explicitly measures which layers of the mechanism the
   transcriptome can reach, rather than assuming full recovery.

### Key findings

| Layer | Recovery | Best evidence |
|-------|----------|---------------|
| Effector arm (protrusion → cytotoxicity) | ✅ Robust | r = 0.55 (bulk), r = 0.32 (single-cell), replicates in 3 independent gastric cohorts |
| Metabolic coupling (SM-balance → protrusion) | ⚠️ Cell-type-resolved only | Invisible in bulk (r = −0.02), recovers in single NK cells (r = +0.03, p = 6×10⁻³) |
| Physical topology phenotype | ❌ Not recovered | Intratumoral NK show *higher* protrusion-machinery transcript — opposite to the physical collapse |

---

## Repository layout

```
GC-NKGraph-Atlas/
├── README.md
├── LICENSE                             # MIT
├── environment.yml                     # Conda environment (GPU-capable)
├── requirements.txt                    # Minimal pip requirements
├── pyproject.toml                      # Python project configuration
│
├── configs/                            # YAML configuration files
│   ├── data_config.yaml
│   ├── experiment_config.yaml
│   ├── graph_config.yaml
│   ├── model_config.yaml
│   ├── scrna_config.yaml
│   ├── sst_axis_config.yaml
│   ├── topology_schema.yaml
│   └── mechanism_cards/                # Reusable mechanism specifications
│       ├── registry.yaml               # Active mechanism card selection
│       ├── mechanism_card.template.yaml # Template for new mechanisms
│       └── zheng_nk_sm_topology.yaml    # Zheng 2023 NK SM-topology card
│
├── src/                                # Source modules
│   ├── common/                         # Config, I/O, logging, seeds, synthetic data
│   ├── data_download/                  # TCGA + GEO downloaders
│   ├── preprocessing/                  # Bulk preprocessing, GEO probe→gene fix
│   ├── immune_scoring/                 # NK-state scoring & labels
│   ├── scrna_analysis/                 # scRNA pipeline (scVI, Leiden) + QC
│   ├── topology/                       # SST-axis proxy, validation, gastric extension
│   ├── graph_construction/             # Heterogeneous gene graph builder
│   ├── baselines/                      # 6 tabular baselines + model comparison
│   ├── models/                         # GC-NKGraph-Atlas GNN (HGT + MLP)
│   ├── interpretation/                 # Candidate prioritization + de-circularization
│   ├── figures/                        # Publication figure generation
│   └── pipeline.py                     # Master launcher (checkpointed)
│
├── tests/                              # pytest unit tests
├── manuscript/                         # Manuscript, figures, notes, cover letter
├── submission_package/                 # Self-contained paper snapshot
├── data/                               # Raw/processed (git-ignored) + synthetic
└── results/                            # Generated tables/figures/logs (git-ignored)
```

---

## Installation

### Option A: Full environment (recommended, GPU-capable)

```bash
conda env create -f environment.yml
conda activate gc-nkgraph
```

### Option B: Minimal (CPU, pip)

```bash
pip install -r requirements.txt
# Then install PyTorch and PyG for your platform:
#   https://pytorch.org/get-started/
#   https://pyg.org/whl/
```

**Requirements:** Python ≥ 3.10. Core analyses (bulk scoring, SST axis, baselines,
GEO validation, figures) run CPU-only. scVI integration and the heterogeneous graph
transformer (HGT) benefit from a GPU but fall back to CPU / spectral encoding.

### Verify installation

```bash
pytest -q                                # Unit tests
python src/pipeline.py --synthetic       # End-to-end on toy data
```

---

## Quick start (synthetic data)

No real data needed — the synthetic mode exercises every code path:

```bash
# Full pipeline on generated toy data (checkpointed; ~2 min CPU)
python src/pipeline.py --synthetic

# Unit tests
pytest -q

# Generate all publication figures from synthetic results
python src/figures/make_figures.py
```

This produces a complete set of tables and figures under `results/` that mirrors
the real-data output structure. All file paths are identical; switching to real
data means replacing `configs/data_config.yaml` source paths.

---

## Data

No data are committed to this repository. All inputs are publicly available:

| Dataset | Source | Samples / Cells | Role |
|---------|--------|-----------------|------|
| TCGA-LIHC | [GDC](https://portal.gdc.cancer.gov/) / [UCSC Xena](https://xenabrowser.net/) | 423 | Positive control (Arm A) |
| TCGA-STAD | [GDC](https://portal.gdc.cancer.gov/) / [UCSC Xena](https://xenabrowser.net/) | 450 | Training (Arm B) |
| GSE62254 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62254) | 300 | External validation |
| GSE84437 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437) | 483 | External validation |
| GSE246662 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246662) | 166,829 (8,310 NK) | scRNA NK atlas |

### Downloading the data

```bash
# TCGA (requires GDC data-transfer tool or UCSC Xena download)
python src/data_download/download_tcga.py --cohorts LIHC STAD

# GEO
python src/data_download/download_geo.py --accessions GSE62254 GSE84437 GSE246662
```

Place files under `data/raw/` per the layout in `configs/data_config.yaml`, or
edit the config to point at your existing copies.

---

## Reproducing the analysis

### One command (checkpointed — skips completed stages)

```bash
python src/pipeline.py --config configs/experiment_config.yaml
```

### Stage by stage (mirrors the manuscript)

```bash
# STAGE I — DATA
python src/data_download/download_tcga.py
python src/data_download/download_geo.py
python src/preprocessing/run_bulk_preprocessing.py
python src/immune_scoring/nk_scores.py

# STAGE II — scRNA
python src/scrna_analysis/run_scrna_v2.py
python src/scrna_analysis/qc_filter.py --in data/processed/scrna/gc_integrated.h5ad

# STAGE III — GRAPH
python src/graph_construction/build_heterograph.py

# STAGE IV — MODEL (Table 3, Fig 4)
python src/baselines/run_all_baselines.py
python src/models/gc_nkgraph_atlas.py
python src/baselines/run_model_comparison.py

# STAGE V — TARGETS (Table 2, 4, 5; Fig 1–3)
python src/topology/sst_axis.py
python src/topology/sst_axis_validation.py
python src/interpretation/prioritize_targets.py
python src/interpretation/split_target_lists.py
python src/preprocessing/run_geo_external_validation.py

# FIGURES
python src/figures/make_figures.py
```

Outputs land in `results/tables/` and `results/figures/`. Each script logs to
`results/logs/LOG.md`. A curated snapshot of key outputs is kept under
`submission_package/`.

---

## Mechanism cards

The pipeline is driven by machine-readable mechanism cards under
`configs/mechanism_cards/`. A mechanism card is a YAML document that declares, for
one published wet-lab mechanism:

- Its molecular chain, gene modules, and expected directions
- Cell-type attribution rules for each module
- Physical ground-truth targets that are gated (out of scope for transcriptomes)
- Graph node/edge types it introduces
- Pre-registered validation hypotheses with an explicit recovery definition

**Applying the framework to a new published mechanism requires only authoring a
new card** — not editing the pipeline core.

### Active card

| Card | Mechanism | Reference |
|------|-----------|-----------|
| `zheng_nk_sm_topology` | Serine→sphingomyelin→membrane-topology→cytotoxicity axis of NK immune evasion | Zheng et al., *Nat Immunol* 2023 |

### Adding a new mechanism

```bash
cp configs/mechanism_cards/mechanism_card.template.yaml \
   configs/mechanism_cards/your_mechanism.yaml
# Edit the card with your mechanism's specifics
# Update configs/mechanism_cards/registry.yaml to activate it
```

---

## Key result files

| File | Content |
|------|---------|
| `results/tables/sst_axis_positive_control_recovery.tsv` | Arm A hypothesis outcomes (bulk + single-cell, H1–H5) |
| `results/tables/sst_axis_scrna_by_tissue.tsv` | Cross-tissue NK module means |
| `results/tables/model_comparison.tsv` | GNN vs 6 baselines, 5-fold CV |
| `results/tables/model_comparison_stats.tsv` | Paired statistical tests (Wilcoxon + t-test) |
| `results/tables/external_validation_results.tsv` | Axis replication in GSE62254/GSE84437 |
| `results/tables/tumor_intrinsic_candidates.tsv` | De-circularized target list (n = 37) |
| `results/tables/axis_confirmation_panel.tsv` | NK-side readout panel (n = 36, labeled as such) |
| `results/tables/candidate_evidence_matrix.tsv` | Full evidence matrix for all candidates |
| `results/figures/fig1_armA_positive_control.pdf` | Arm A — partial recovery of SST axis in liver |
| `results/figures/fig2_armB_extension.pdf` | Arm B — gastric extension + external validation |
| `results/figures/fig3_targets.pdf` | De-circularized tumor-intrinsic target prioritization |
| `results/figures/fig4_model_comparison.pdf` | Model comparison — GNN vs 6 baselines |
| `results/figures/fig5_mechanism_card_concept.pdf` | Mechanism-card abstraction — concept diagram |

---

## Citation

If you use this code or the mechanism-card formalism in your work, please cite:

```bibtex
@article{lyu2026gcnkgraph,
  title   = {{GC-NKGraph-Atlas}: reconstructing the serine--sphingomyelin--membrane-topology
             axis of {NK}-cell immune evasion from tumor transcriptomes},
  author  = {Lyu, Guohao and Xia, Yingchun and Liu, Huichao and Zhu, Xiaolei and
             Yang, Shuai and Zhou, Ailian and Gu, Lichuan},
  journal = {Briefings in Bioinformatics},
  year    = {2026},
  note    = {Under review}
}
```

BibTeX will be updated with volume/page/DOI upon acceptance.

---

## License

Released under the [MIT License](LICENSE). © 2026 Guohao Lyu, Lichuan Gu, Ailian
Zhou, and contributors.

---

## Contact

**Corresponding authors:**
- Prof. Lichuan Gu — [glc@ahau.edu.cn](mailto:glc@ahau.edu.cn)
- Prof. Ailian Zhou — [zhouailian@caas.cn](mailto:zhouailian@caas.cn)

**Code issues:** [GitHub Issues](https://github.com/nblvguohao/GC-NKGraph-Atlas/issues)

School of Artificial Intelligence, Anhui Agricultural University, Hefei 230036, China
