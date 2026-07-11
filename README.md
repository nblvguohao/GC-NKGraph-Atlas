# GC-NKGraph-Atlas

**Mapping the transcriptional reach of the serine–sphingomyelin–membrane-topology
axis of NK-cell immune evasion with a single-cell-informed heterogeneous graph
framework, from liver to gastric cancer.**

This repository accompanies a manuscript prepared for *Briefings in
Bioinformatics*. It operationalizes a published wet-lab immune-evasion mechanism
(Zheng et al., *Nat Immunol* 2023) as a reusable **mechanism-card**, maps the
axis's transcriptional reach under a two-arm design (liver positive control +
gastric-cancer extension), benchmarks a mechanism-grounded heterogeneous graph
model against tabular baselines, and prioritizes putative tumor-intrinsic
candidate targets.

> **Honest scope.** The framework produces a *map of the mechanism's
> transcriptional reach*, not a blanket recovery: the effector layer
> (protrusion→cytotoxicity) recovers and generalizes; the metabolic coupling is
> weakly detectable only at single-cell resolution; the physical topology
> phenotype is not captured by machinery transcription. See the manuscript for
> the full framing.

Repo: `https://github.com/nblvguohao/GC-NKGraph-Atlas`

> **🌐 Interactive Mechanism Card Playground:** [Open in browser](https://nblvguohao.github.io/GC-NKGraph-Atlas/web/) —
> browse all mechanism cards, compare gene modules across mechanisms, and explore
> the prioritized target list with external database evidence (DepMap, DrugBank, Open Targets).

---

## Highlights

- **4 mechanism cards** (serine–SM–topology, adenosine–A2AR, TGFβ–SMAD–ECM, NKG2D–MICA/B shedding) demonstrating reusable formalism — one (serine–SM–topology) executed end-to-end, the others authored to demonstrate the format
- **37 tumor-intrinsic candidate targets** with DepMap/DrugBank/Open Targets cross-validation
- **Edge-type ablation study** quantifying each graph edge's contribution to model performance
- **Two-arm design** (liver positive control + gastric cancer extension) with pre-registered hypotheses
- **Interactive web playground** (`web/index.html`) — browse cards and targets without installing anything

---

## Repository layout

```
configs/               YAML configs + mechanism cards (configs/mechanism_cards/)
src/
  data_download/       TCGA + GEO downloaders
  preprocessing/       bulk preprocessing; GEO probe->gene fix; external validation
  immune_scoring/      NK-state scoring & labels
  scrna_analysis/      scRNA pipeline (v2) + QC hardening
  topology/            SST-axis proxy, gastric extension, validation
  graph_construction/  heterogeneous gene graph builder
  baselines/           6 tabular baselines + GNN-vs-baseline comparison
  models/              GC-NKGraph-Atlas GNN
  interpretation/      candidate prioritization + de-circularized target split
  figures/             publication figure generation
  common/              config/IO/logging/seeds + synthetic data
  pipeline.py          master launcher (checkpointed; supports --synthetic)
tests/                 pytest unit tests
manuscript/            main text, cover letter, BiB checklist, figures/, notes/
submission_package/    self-contained paper snapshot (tables, figures, scripts)
results/               generated tables/figures/logs (git-ignored; regenerated)
data/                  raw/interim/processed (git-ignored; see Data below)
```

---

## Installation

**Recommended (full environment, GPU-capable):**

```bash
conda env create -f environment.yml
conda activate gc-nkgraph
```

**Minimal (pip, CPU):**

```bash
pip install -r requirements.txt
# then install PyTorch + PyG for your platform:
#   https://pytorch.org/get-started/   https://pyg.org
```

Python 3.10. Core analyses (bulk scoring, SST axis, baselines, GEO validation,
figures) run CPU-only; scVI integration and the HGT graph model benefit from a
GPU but fall back to CPU / spectral encoding.

---

## Data

No data are committed (raw/interim/processed and `*.h5ad`/`*.tsv` are git-ignored).
All inputs are public:

| Dataset | Source | Role |
|---------|--------|------|
| TCGA-LIHC, TCGA-STAD | GDC (https://portal.gdc.cancer.gov/) / UCSC Xena | positive control / train |
| GSE62254, GSE84437 | GEO (https://www.ncbi.nlm.nih.gov/geo/) | external validation |
| GSE246662 | GEO | scRNA NK atlas |

Download helpers: `python src/data_download/download_tcga.py` and
`python src/data_download/download_geo.py` (or place files under `data/raw/` per
`configs/data_config.yaml`).

**Quick check without any real data** — synthetic mode exercises the full code
path:

```bash
python src/pipeline.py --synthetic      # end-to-end on generated toy data
pytest -q                               # unit tests
```

---

## Reproducing the analysis

Run the master pipeline (checkpointed; skips completed stages):

```bash
python src/pipeline.py --config configs/experiment_config.yaml
```

…or stage by stage (mirrors the manuscript):

```bash
# 1. Bulk preprocessing + NK-state labels
python src/preprocessing/run_bulk_preprocessing.py
python src/immune_scoring/nk_scores.py

# 2. scRNA pipeline + QC hardening
python src/scrna_analysis/run_scrna_v2.py
python src/scrna_analysis/qc_filter.py --in data/processed/scrna/gc_integrated.h5ad

# 3. SST-axis proxy + positive-control / single-cell tests (Arm A, Table 2)
python src/topology/sst_axis.py
python src/topology/sst_axis_validation.py

# 4. Heterogeneous graph + models (Table 3, Fig 4)
python src/graph_construction/build_heterograph.py
python src/baselines/run_all_baselines.py
python src/models/gc_nkgraph_atlas.py
python src/baselines/run_model_comparison.py        # GNN vs baselines + paired tests

# 5. Targets — de-circularized (Table 4, Fig 3)
python src/interpretation/prioritize_targets.py
python src/interpretation/split_target_lists.py

# 6. External validation in gastric cohorts (Arm B, Table 5)
python src/preprocessing/run_geo_external_validation.py

# 7. Publication figures (Fig 1–4)
python src/figures/make_figures.py
```

Outputs land in `results/tables/` and `results/figures/`. A curated snapshot of
the key tables and figures is kept under `submission_package/`.

---

## Key result files

| File | Content |
|------|---------|
| `results/tables/sst_axis_positive_control_recovery.tsv` | Arm A hypothesis outcomes (bulk + single-cell) |
| `results/tables/model_comparison.tsv` / `_stats.tsv` | GNN vs baselines + paired significance |
| `results/tables/external_validation_results.tsv` | axis replication in GSE62254/GSE84437 |
| `results/tables/tumor_intrinsic_candidates.tsv` | de-circularized target list (n=37) |
| `results/figures/fig1–fig4` | publication figures |

---

## Mechanism cards

The pipeline is driven by machine-readable mechanism cards under
`configs/mechanism_cards/` (`zheng_nk_sm_topology.yaml`, with a `registry.yaml`
and a `mechanism_card.template.yaml`). Applying the framework to a new published
mechanism means writing a new card — not editing the pipeline core.

---

## Citation

If you use this code, please cite the accompanying manuscript (details on
acceptance). Funding: NSFC (32472007, 62301006, 62301008); Natural Science
Foundation of Anhui Province (2308085MF217, 2308085QF202); Anhui Province Key
Laboratory of Intelligent Agricultural Technology and Equipment.

**Contact:** glc@ahau.edu.cn (L. Gu), zhouailian@caas.cn (A. Zhou).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
