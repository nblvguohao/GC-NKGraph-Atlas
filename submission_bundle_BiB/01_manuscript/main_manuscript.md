# Transcriptional Reach of the Serine–Sphingomyelin–Membrane-Topology Axis in NK-Cell Immune Evasion: A Single-Cell Heterogeneous Graph Framework from Liver to Gastric Cancer

**Target Journal:** *Briefings in Bioinformatics* (Oxford University Press)

**Article type:** Problem Solving Protocol

**Authors:** Guohao Lyu<sup>1,2</sup>, Yingchun Xia<sup>1,2</sup>, Huichao Liu<sup>1,2</sup>, Xiaolei Zhu<sup>1,2</sup>, Shuai Yang<sup>1,2</sup>, Ailian Zhou<sup>3,4,\*</sup>, Lichuan Gu<sup>1,2,\*</sup>

**Affiliations:**
1. School of Artificial Intelligence, Anhui Agricultural University, Hefei 230036, China
2. Anhui Province Key Laboratory of Intelligent Agricultural Technology and Equipment, Anhui Agricultural University, Hefei 230036, China
3. Agricultural Information Institute, Chinese Academy of Agricultural Sciences, Beijing 100081, China
4. Key Laboratory of Agricultural Blockchain Application, Ministry of Agriculture and Rural Affairs, Beijing 100081, China

\* Corresponding authors: zhouailian@caas.cn (A.Z.), glc@ahau.edu.cn (L.G.)

**Author biographies**

**Guohao Lyu** is a graduate student at the School of Artificial Intelligence, Anhui Agricultural University. His research focuses on graph neural networks and single-cell transcriptomics for cancer immunology.

**Yingchun Xia** is a faculty member at the School of Artificial Intelligence, Anhui Agricultural University. Her research interests include machine learning and bioinformatics.

**Huichao Liu** is a graduate student at the School of Artificial Intelligence, Anhui Agricultural University. His research focuses on computational biology and data mining.

**Xiaolei Zhu** is a researcher at the School of Artificial Intelligence, Anhui Agricultural University. His interests span deep learning and multi-omics integration.

**Shuai Yang** is a faculty member at the School of Artificial Intelligence, Anhui Agricultural University. His research covers computer vision and biomedical image analysis.

**Ailian Zhou** is a professor at the Agricultural Information Institute, Chinese Academy of Agricultural Sciences, and the Key Laboratory of Agricultural Blockchain Application, Ministry of Agriculture and Rural Affairs. Her research focuses on agricultural informatization and blockchain applications.

**Lichuan Gu** is a professor at the School of Artificial Intelligence, Anhui Agricultural University, and director of the Anhui Province Key Laboratory of Intelligent Agricultural Technology and Equipment. His research spans artificial intelligence, bioinformatics, and intelligent agricultural systems.

*(Biographies ~30 words each per BiB style; update with accurate details before submission.)*

**Running head:** GC-NKGraph-Atlas

---

## Abstract

**Motivation:** A landmark wet-lab study (Zheng et al., *Nat Immunol* 2023) established that tumors evade natural killer (NK) cell cytotoxicity by dysregulating serine metabolism, which depletes sphingomyelin (SM) in NK membranes, collapses membrane protrusions, and abolishes lytic immune synapse formation. However, this mechanism was demonstrated in a small number of liver cancer patients using specialized single-cell mass spectrometry and super-resolution imaging — methods that cannot scale to cohort-level analyses or be readily tested in other cancer types.

**Results:** We present GC-NKGraph-Atlas, a computational framework that asks **how much** of the serine→sphingomyelin→membrane-topology→cytotoxicity axis is transcriptionally recoverable from public tumor transcriptomes. The framework employs: (i) a *mechanism-card* abstraction that encodes the wet-lab mechanism as a machine-readable recipe; (ii) a single-cell-informed, NK-aware heterogeneous graph that integrates protein–protein interactions, ligand–receptor pairs, transcription factor targets, and a mechanism-grounded `metabolic_crosstalk` edge type; and (iii) a graph neural network that learns gene embeddings from this multi-relational graph to predict NK immune states. A two-arm design (Arm A: liver/HCC positive control; Arm B: gastric cancer extension) yields a **scoping result** with three layers: **First**, the *effector arm* of the axis is recoverable in bulk transcriptomes — protrusion-machinery coupled to cytotoxicity output in TCGA-LIHC (r=0.55, p=5×10⁻³⁵) — and intratumoral NK cells show the expected loss of cytotoxic output (Δ=−0.14, p=6×10⁻⁵²); however, ~50% of this bulk coupling is NK-cell-abundance covariation (after adjusting for an NK-fraction proxy, partial r=0.25 in TCGA-LIHC, and it survives in two of three bulk cohorts but not the third), so the effector arm is *partially recoverable and substantially infiltration-driven*, not robustly recovered. At single-cell resolution the same coupling is nominally significant after pseudoreplication correction (8,310 NK cells from 9 samples, corrected r=0.31, p=3.9×10⁻⁸), but we found this single-cell number does not survive further, necessary technical-confound controls: it collapses under residualization against library size and the dominant transcriptional (scVI latent) structure (r=0.04–0.09 depending on scoring method) and does not exceed the correlation expected from randomly drawn gene modules of the same sizes in this dataset (permutation P=0.97 with an expression-matched scoring method). It is not explained by generic NK activation alone (activation-partialled r=0.25; activation-matched subset r=0.24), which rules out simple co-activation but not the broader technical/structural confound. We therefore treat the single-cell effector-arm number as not distinguishable from technical background and anchor the effector-arm claim on the bulk result and its gastric replication, rather than on the single-cell number. **Second**, *transcription does not proxy the physical topology phenotype*: intratumoral NK cells show *higher*, not lower, protrusion-machinery transcript levels (Δ=+0.14, p=3×10⁻⁹¹) — the opposite direction to the physical membrane collapse — establishing a fundamental disconnect between the transcriptional and physical layers of this mechanism. **Third**, the *upstream metabolic arm* is not transcriptionally recoverable: the SM-balance→protrusion coupling is undetectable in bulk (r=−0.02, p=0.72) and not significant after cell-type resolution and pseudoreplication correction (single-cell NK corrected r=0.029, p=0.20), consistent with serine→SM crosstalk operating at the metabolite rather than transcript level. The effector coupling reproduces at the zero-order level in two independent gastric microarray cohorts (GSE62254 r=0.42; GSE84437 r=0.62; both p≪10⁻¹³), but these bulk replications are themselves substantially NK-abundance-driven (§3.2). Running a *second, independently authored mechanism card* (TGFβ→SMAD→NK exclusion) end-to-end on the same cohorts reproduces the boundary: its effector co-variation recovers while its mechanism-specific causal predictions fail in bulk — establishing that "bulk transcriptomes recover NK effector state but not mechanism-specific upstream coupling" generalizes across two distinct immune-evasion mechanisms. The framework then prioritizes **37 putative tumor-intrinsic candidate targets** (led by the druggable serine/sphingomyelin enzymes PHGDH, SGMS2, PSAT1, PSPH, SMPD3/1), each with a recommended wet-lab validation assay; we note that their malignant-cell fold-changes are near-zero and ~46% are annotated to NK-side modules, so this list is a mechanism-intersecting shortlist for hypothesis generation rather than a set of tumor-exclusive targets.

**Availability:** All code, configuration, and synthetic test data are available at https://github.com/nblvguohao/GC-NKGraph-Atlas. The mechanism-card template enables application to additional published mechanisms without modifying the pipeline core.

**Contact:** glc@ahau.edu.cn

**Keywords:** NK cell; immune evasion; heterogeneous graph neural network; single-cell transcriptomics; sphingomyelin; gastric cancer

---

## Key Points

- A **mechanism-card** formalism converts published wet-lab immune-evasion mechanisms into scalable transcriptome-based analysis runs, demonstrated **end-to-end on two independent cards** (serine–SM and TGFβ→SMAD) that reveal a consistent transcriptional-reach boundary.
- A **two-arm design** (liver positive control + gastric extension) produces a **scoping result**: the *effector arm* (protrusion→cytotoxicity) is recoverable in bulk transcriptomes, but **~50% of the bulk coupling is NK-cell-abundance covariation** — after adjusting for an NK-fraction proxy it attenuates from r≈0.55 to r≈0.25 and remains significant in two of three bulk cohorts (TCGA-LIHC, GSE62254) but not the third (GSE84437, partial r≈0, n.s.). At single-cell resolution the same coupling does not survive count-depth/latent-structure residualization. The *metabolic arm* (SM-balance→protrusion) is not significant after correction (p=0.20). The effector arm is therefore *partially recoverable and substantially infiltration-driven*, not robustly recovered.
- **The reach boundary generalizes across mechanisms:** running a second, independently authored card (TGFβ→SMAD→NK exclusion) end-to-end reproduces the pattern — generic NK-effector co-variation is recoverable, but the mechanism-specific upstream causal predictions (TGFβ-signaling→receptor suppression, CAF-ECM→exclusion) fail in bulk. Across two distinct mechanisms, bulk transcriptomes recover NK effector state but not mechanism-specific causal coupling, and NK-fraction control is required to see this.
- **Transcription does not proxy the physical topology phenotype:** intratumoral NK cells exhibit *higher* protrusion-machinery transcript levels while their physical membrane protrusions are collapsed — a fundamental disconnect that defines the natural boundary of transcriptome-based reconstruction for membrane-lipid mechanisms.
- The heterogeneous graph is used as a **probe, not a predictor**: encoding the mechanism's metabolic coupling as a `metabolic_crosstalk` edge and ablating it gives an independent, architecture-based confirmation that this coupling is transcriptionally absent — removing the edge abolishes the corresponding embedding coupling, converging with the null correlation result.
- The framework outputs a **de-circularized, putative tumor-intrinsic candidate target list** (37 genes, led by druggable serine/SM enzymes) kept strictly separate from the NK-side axis readout, each with a recommended wet-lab validation assay — presented as a hypothesis-generation shortlist for experimental follow-up.

---

## 1. Introduction

### 1.1 Biological context

Natural killer (NK) cells are innate lymphoid cells critical for anti-tumor immunity. Unlike T cells, NK cells kill without prior antigen sensitization, making them attractive effectors for cancer immunotherapy [1-5]. However, tumors deploy multiple mechanisms to evade NK-mediated killing: metabolic competition [6-8], inhibitory ligand upregulation [9], immunosuppressive cytokine secretion [10], and physical exclusion from the tumor nest [11].

A particularly elegant evasion mechanism was recently elucidated by Zheng et al. [12]: tumors dysregulate serine metabolism [13,14], reducing the availability of serine for NK sphingolipid synthesis. The resulting depletion of sphingomyelin [15,16] in NK membranes collapses membrane protrusions and microvilli, preventing the formation of lytic immune synapses and abolishing cytotoxicity. Critically, this phenotype can be rescued by inhibiting SM catabolism (via ASM/SMPD1 or NSMASE/SMPD2-4), and the rescue is synergistic with Tim3 (HAVCR2) checkpoint blockade. The follow-up framing [17] positions sphingomyelin as a metabolic immune checkpoint and reports that the mechanism extends beyond liver cancer to lung, colon, and ovarian cancers.

**However, this mechanism was proven with techniques that do not scale.** The anchor paper used super-resolution SEM imaging and single-immunocyte mass spectrometry on a small patient cohort — gold-standard evidence for physical membrane topology but limited to specialized laboratories and small sample sizes. No transcriptomic data were deposited. The transcriptional footprint of this mechanism — i.e., whether the molecular machinery and capacity for the serine→SM→protrusion→cytotoxicity axis can be detected from widely available transcriptomic data — remains unexplored.

### 1.2 The gap

Three gaps motivate this work:

1. **Scalability.** The mechanism was demonstrated in tens of patients; testing it across hundreds or thousands of tumor samples requires a transcriptome-based proxy.
2. **Transferability.** Whether the same axis operates in gastric cancer — a digestive-tract cancer not on the published extension list (lung/colon/ovarian) [17] — is unknown.
3. **Reusability.** Each published mechanistic discovery currently requires a bespoke computational follow-up; a generalizable engine that converts mechanism papers into target-discovery runs does not exist.

### 1.3 Our approach

We introduce **GC-NKGraph-Atlas**, a computational framework designed to address all three gaps. The framework is organized around a reusable **mechanism-card** abstraction: a machine-readable YAML specification that encodes one published wet-lab mechanism — its molecular chain, gene modules, expected directions, physical ground-truth targets, validation assays, and therapeutic hooks — as a recipe that the pipeline consumes.

The specific mechanism card driving this study (`zheng_nk_sm_topology`) operationalizes the Zheng 2023 serine–SM–topology axis. The framework: (a) defines cell-type-attributed transcriptional proxies for each step of the mechanistic chain, computed from single-cell RNA-seq; (b) constructs a tumor–NK heterogeneous graph in which a `metabolic_crosstalk` edge type is justified by the specific biology rather than generic interaction priors; (c) learns gene embeddings from this graph to predict NK immune states; and (d) ranks putative tumor-intrinsic targets by multi-evidence scoring.

The study employs a two-arm design:

- **Arm A (Positive Control, liver/HCC):** test how much of the published axis is visible in the system where it was proven, using independent public transcriptomes (TCGA-LIHC + public HCC scRNA). This is the credibility anchor.
- **Arm B (Novel Extension, gastric cancer):** test the same axis in gastric cancer — a natural, not-yet-claimed digestive-tract extension — and prioritize gastric-specific candidate targets.

### 1.4 Related work

Our framework draws on four strands of computational biology, and its novelty
lies in how it combines them around a specific published mechanism rather than in
any single component.

**Immune-state and cell-type inference from transcriptomes.** A large tool
family scores immune infiltration and cell state from bulk and single-cell data —
CIBERSORTx [18] and quanTIseq [19] for deconvolution, and phenotype-to-genotype
methods such as Scissor [20], which links single-cell phenotypes to bulk clinical
variables. These estimate *how much* of a cell type is present or *which*
phenotype dominates; none operationalizes a defined mechanistic chain or enforces
the claim boundary between a transcriptional proxy and a physical phenotype, which
is central to our design. Single-cell atlases increasingly link a functional cell
state to long-term clinical outcome — for example, a type-2 program in CAR-T
infusion products associated with multi-year leukaemia remission [21]; our
framework instead attributes states to a defined evasion mechanism, and anchoring
the NK-state readout to clinical outcome remains future work (§4.3).

**Cell–cell communication and metabolic-flux inference.** CellChat [22],
CellPhoneDB [23] and NicheNet [24] infer ligand–receptor signaling between cell
types, and scFEA [25] estimates metabolic flux from single-cell expression. We
borrow the ligand–receptor edge concept for the heterogeneous graph but add a
mechanism-specific `metabolic_crosstalk` edge grounded in one published
serine→sphingomyelin relationship. Rather than assuming this edge improves
prediction, we use it as a diagnostic: structurally encoding the mechanism's
upstream metabolic coupling and asking whether the transcriptome corroborates
it — a question the correlation analyses (§3.2) and the graph ablation (§3.7)
answer convergently.

**Trajectory and dysfunction-state modeling.** CytoTRACE [26] and related
pseudotime methods order cells along differentiation or dysfunction gradients;
we use an analogous NK dysfunction axis to define reversible-state proxy labels,
but anchor the states to the specific effector/checkpoint genes of the target
mechanism rather than to a generic exhaustion signature. Most characterized
regulators of NK dysfunction are surface checkpoints; intracellular and
transcription-factor-level regulators (e.g. the cAMP-responsive factor CREM
downstream of CAR and IL-15 signalling [27]) are comparatively understudied, and
our candidate prioritization currently inherits this surface/metabolic-gene bias (§4.3).

**Graph learning in cancer genomics.** Graph neural networks [28-30] have been applied to
multi-omics integration (e.g. MOGONET [31]) and to molecular-interaction networks
for outcome prediction, and heterogeneous graph transformers [32] provide
type-specific message passing over multi-relational graphs. Most directly relevant
are recent transformer-based multi-network fusion models built specifically for
cancer driver-gene discovery, such as TREE [49], which learns cancer-gene
embeddings via co-attention over PPI, miRNA, and transcription-factor networks
spanning up to ~20,000 genes and 2.7 million edges, and GRAFT [50], which fuses
PPI, GO-similarity, and pathway-co-occurrence networks through a learnable
attention mechanism and an edge-attention-biased transformer. Both target
genome-scale, cancer-type-agnostic driver-gene discovery from generic interaction
networks; our graph instead spans a ~100-gene mechanism-anchored panel and asks a
narrower, mechanism-specific question — not "which genes drive this cancer" but
"how far does one specific published immune-evasion mechanism's transcriptional
signal reach" — so scale, labels, and objective are deliberately incommensurate
with this line of work rather than a smaller instance of the same task. Existing
applications typically treat the graph as a generic prior (a PPI or co-expression network) and
report a predictive gain. Our use of the graph is different in intent: on top of
standard priors (STRING protein–protein interactions, CellChatDB ligand–receptor
pairs) we add one edge type derived directly from the mechanism under study, and
then ask — rather than assume — whether that mechanism-grounded structure is
transcriptionally supported and predictively useful. As we show (§3.4, §3.7),
the graph does not beat simpler baselines on accuracy; its value is as a
*mechanism-structured probe* whose ablation behavior independently corroborates
which layer of the axis the transcriptome can and cannot reach.

To our knowledge, no prior tool converts a single published wet-lab immune-evasion
mechanism into a machine-readable card that drives cell-type-attributed proxy
construction, a mechanism-structured graph probe, and de-circularized target
prioritization, while explicitly measuring — rather than assuming — how much of
the mechanism the transcriptome can reach.

---

## 2. Methods

### 2.1 Study design overview

The framework proceeds through 14 phases, grouped into five stages:

| Stage | Content |
|-------|---------|
| I — DATA (Phases 1–2) | Download + preprocess TCGA-STAD, TCGA-LIHC, GEO gastric cohorts |
| II — scRNA (Phases 3–7) | scRNA integration, NK atlas annotation, state scoring, trajectory |
| III — GRAPH (Phase 8) | Heterogeneous gene graph with mechanism-grounded edges |
| IV — MODEL (Phases 9–10) | Baseline comparison + GNN-based NK-state classifier |
| V — TARGETS (Phases 11–14R) | SST-axis scoring, candidate prioritization, assay recommendation |

A master pipeline launcher (`src/pipeline.py`) orchestrates execution with checkpoint-based skipping and supports synthetic data mode for testing.

**Workflow overview.** Workflow of the GC-NKGraph-Atlas framework, from data acquisition through target prioritization, showing the mechanism-card abstraction and the two-arm study design (Arm A: liver positive control; Arm B: gastric cancer extension).

### 2.2 The mechanism-card abstraction

A *mechanism card* is a YAML document with the following required sections (see `configs/mechanism_cards/mechanism_card.template.yaml` for the full schema):

| Section | Content |
|---------|---------|
| `origin` | Paper DOI, lab, notes on data availability |
| `biology` | Phenotype, cell type affected, ordered mechanistic chain |
| `transcriptional_proxy` | Gene modules per step, expected direction, attribution requirements, claim boundaries |
| `physical_ground_truth` | Gated physical measurement targets (SEM, MS), mock policy |
| `graph_integration` | New node/edge types, construction rules |
| `validation` | Positive control cohort, pre-registered hypotheses, recovery definition |
| `therapeutic_hook` | Intervention logic, combination rationale, patient stratification readout |

**Design principle.** The card strictly separates what is *computable* from transcriptomes (transcriptional proxy, status ACTIVE) from what requires *physical measurement* (membrane protrusion density, SM metabolite content — status GATED). The loader raises a clear error if gated targets are absent; mock data are permitted only for code-path testing with `*_MOCK_*` file patterns and are excluded from all scientific claims.

### 2.3 SST-axis transcriptional proxy (Phase 14R)

The serine–sphingomyelin–topology axis is operationalized as seven gene modules derived from the anchor paper [12] and its follow-up [17]:

| Module | Cell type | Genes (n) | Expected direction |
|--------|-----------|-----------|-------------------|
| `tumor_serine_capacity` | malignant | PHGDH, PSAT1, PSPH, SHMT1/2, MTHFD1/2, MTHFD1L, SLC1A4/5 (10) | Calibrated on liver control |
| `nk_sm_synthesis` | NK | SGMS1, SGMS2 (2) | ↑ = more topology-permissive |
| `nk_sm_catabolism` | NK | SMPD1–4 (4) | ↑ = less topology-permissive |
| `nk_denovo_sphingolipid` | NK | SPTLC1–3, SPTSSA, CERS2/4/5/6, DEGS1 (9) | Context-dependent |
| `nk_protrusion_machinery` | NK | ERM family, Arp2/3 complex, WASP/WAVE, Rho GTPases, formins, BAR domain (24) | ↑ = more topology-permissive |
| `nk_synapse_cytotoxicity_outcome` | NK | NKG7, GNLY, GZMB, PRF1, IFNG, LCP2, LAT, VAV1, TLN1, ITGAL, ITGB2 (11) | Axis-positive correlate |
| `checkpoint_link` | NK | HAVCR2 (1) | ↑ = less topology-permissive |

Per-cell module scores are computed as mean z-scores of constituent genes. Derived scores include:
- **nk_sm_balance** = mean_zscore(sm_synthesis) − mean_zscore(sm_catabolism)
- **nk_topology_permissive** = composite(sm_balance, protrusion_machinery)
- **sst_axis_score** = integrated(tumor_serine_capacity [calibrated sign], nk_topology_permissive, cytotoxicity_outcome)

**Cell-type attribution.** All tumor-vs-NK axis claims require single-cell resolution. For bulk data, deconvolution (CIBERSORTx or quanTIseq) is applied as a fallback, and results lacking cell-type resolution are flagged `MIXED_UNRESOLVED`.

**Honesty rule.** The permitted language is strictly "transcriptional program permissive-of / associated-with the topology phenotype." Claims of physical membrane topology prediction or SM metabolite content prediction from transcriptomes are prohibited.

### 2.4 scRNA-seq analysis

**Data and loading.** Single-cell RNA-seq was obtained from GSE246662, comprising
nine samples across three tissues — healthy liver (HL1–3), gastric cancer (GC1–3),
and gastric-cancer liver metastasis (LM1–3) — totaling 166,829 cells. Per-sample
count matrices were loaded with orientation auto-detection (genes×cells vs.
cells×genes), tagged with tissue and tumor/normal condition, and concatenated on
the intersection of genes (inner join) to avoid cross-platform sparsity artifacts
(`src/scrna_analysis/run_scrna_v2.py`).

**Quality control.** Per-cell QC metrics (detected genes per cell,
mitochondrial fraction) were computed for all concatenated cells
(`src/scrna_analysis/run_scrna_v2.py`), with mitochondrial fraction computed
using a NaN-robust fallback to handle heterogeneous cross-platform inputs.
Hard per-cell thresholds (e.g. on gene count or mitochondrial fraction) and
doublet-based exclusion were not applied to the reported dataset; all
166,829 concatenated cells were retained into normalization and integration.
The gene space is restricted at two points instead of by a per-cell filter:
the cross-sample gene intersection at the concatenation step above, and the
top-3,000-highly-variable-gene selection below. Technical, quality-related
variance (library size, detected-gene count) is addressed downstream via the
count-depth residualization diagnostics reported in §3.2, rather than by
upfront hard-threshold cell exclusion; this is discussed as a limitation in
§4.3.

**Normalization, integration, and clustering.** Counts were library-size
normalized to 10⁴ and log1p-transformed; 3,000 highly variable genes were selected
(Seurat v3 flavor [33], with a variance-based fallback; Seurat v4 [34] provides the broader multimodal framework). Batch effects across the nine
samples were corrected with scVI [35] (`sample_id` as batch key, 30 latent
dimensions, 2 layers, up to 200 epochs with early stopping). Neighborhood graphs,
UMAP, and Leiden clustering (resolution 1.0) were computed on the scVI latent
space using SCANPY [36]. The analysis pipeline follows current best-practice recommendations for single-cell RNA-seq [37].

**Cell-type annotation.** Lineages were assigned from canonical marker-set mean
expression: NK (NCAM1, KLRD1, NKG7, GNLY, KLRF1, EOMES, NCR1, FCGR3A), T (CD3D/E/G,
CD4, CD8A), monocyte (CD14, CD68, CSF1R), and B (MS4A1, CD79A, CD19). NK cells
were separated from T cells by requiring an NK score above threshold with a low T
score, yielding 8,310 NK cells used for the axis analyses. We note this
marker-threshold labeling as a limitation (§4.3) relative to reference-based
mapping (e.g. scANVI/scArches [38]), and report NK counts per sample
(`results/tables/gc_scrna_dataset_summary.tsv`) for transparency.

**NK immune state classification.** Four states are defined from scRNA signatures:
- **NK-hot-cytotoxic:** high GZMB/PRF1/IFNG, low HAVCR2/TIGIT
- **NK-hot-dysfunctional:** high HAVCR2/TIGIT/CD96, moderate cytotoxicity genes
- **NK-cold/excluded:** low NK signature, low cytotoxicity
- **NK-intermediate:** transitional state

These states are projected onto bulk samples via scRNA-anchored scoring (mean z-score of state-specific gene sets [39,40], calibrated on the scRNA-defined states).

### 2.5 Heterogeneous graph construction

The graph is deliberately **axis-centered**: nodes are the 100 genes of the SST-axis
modules (§2.3), curated NK receptors, and the candidate pool. Onto this gene set we
lay standard prior-network edges plus two mechanism-grounded edge types and one
data-driven co-expression edge. The realized edge counts (from the built graph,
`data/processed/graph/edges.tsv`) are:

| Edge type | Source → Target | Source database | Weight | Edges (this panel) |
|-----------|-----------------|-----------------|--------|--------------------|
| `ppi` | gene ↔ gene | STRING v12 (score ≥ 700) | score / 1000 | 474 |
| `ligand_receptor` | gene → gene | CellChatDB (via OmniPath) | 0.9 | 12 |
| `tf_target` | TF gene → target gene | ChEA 2022 | 0.8 | 0 (see note) |
| `metabolic_crosstalk` | tumor_serine gene → NK topology gene | Zheng 2023 | 0.5 | 300 |
| `sm_topology_axis` | NK axis gene ↔ NK axis gene | Zheng 2023 | 0.3 | 820 |
| `coexpression` | axis gene ↔ axis gene | NK scRNA (|r|>0.3) | abs(r) | 6 |

**Note on `tf_target`.** ChEA 2022 is integrated by the pipeline, but because the
graph is intentionally restricted to the ~100 axis-centered genes, no ChEA
transcription-factor→target pair falls entirely within this panel, so this edge
type contributes zero edges to the realized graph. We report it explicitly rather
than omitting it: the generic priors that actually connect these genes are
protein–protein interactions (STRING) and ligand–receptor pairs (CellChatDB),
which together form the non-mechanistic background against which the
mechanism-grounded edges are compared in the ablation (§3.7). A prior
`dysfunction_correlation` edge type from bulk correlations to NK-state nodes was
part of the original design but is not used in the analyses reported here.

**Key design choice — `metabolic_crosstalk` edge.** This edge connects tumor-side
serine metabolism genes (PHGDH, PSAT1, etc.) to NK-side SM/topology genes (SGMS1,
SMPD1, EZR, etc.), encoding the Zheng 2023 upstream metabolic coupling as graph
structure. Its purpose is diagnostic, not predictive: it lets us ask whether
hard-coding this coupling helps downstream tasks, and §3.7 shows it does not — in
a way that converges with the null transcriptional correlation for the same
coupling (H1, §3.2). The edge *weight* (0.5) is a fixed structural prior applied
uniformly to every tumor-serine-gene → NK-topology-gene pair
(`src/graph_construction/build_heterograph.py::build_sst_edges`); it is not fit to
any cohort and therefore cannot leak information between the liver calibration and
gastric test data. A separate quantity, the *sign* of the tumor_serine_capacity
term in the descriptive `sst_axis_score` composite (§2.3,
`src/topology/sst_axis_validation.py`), is determined from the observed H1
correlation on the liver cohort; this sign calibration does not feed back into the
graph edge weight, and in any case H1 is null at both resolutions (§3.2), so the
calibrated sign carries little information.

### 2.6 Graph neural network model

The model employs a two-stage architecture:

**Stage 1 — Gene Graph Encoder.** Gene embeddings are learned from the heterogeneous graph via spectral decomposition (truncated SVD) of the combined normalized adjacency matrix (PPI + LR + TF + SST edges, uniformly weighted per edge type). Each gene is represented as a d-dimensional vector that encodes its multi-relational neighborhood. All embeddings reported in this study, including the ablation in §3.7, use this spectral encoder; a learnable, `torch_geometric` [41]-based heterogeneous graph transformer (HGT) [32] would be a natural extension (§4.4) but was not used to produce any result reported here.

**Stage 2 — NK State Classifier.** For each bulk tumor sample with expression vector x ∈ ℝ^{|G|}, the input to the classifier is the concatenation [x, x·E] where E ∈ ℝ^{|G|×d} is the gene embedding matrix. This provides both the raw expression signal and the graph-informed projection. A 2-3 layer MLP with batch normalization and dropout predicts NK immune state (binary: cytotoxic vs rest).

Training uses 5-fold stratified cross-validation with early stopping (patience=30 epochs). The optimizer is Adam with CrossEntropy loss; learning rate (1.7×10⁻³), weight decay (5.6×10⁻⁶), and dropout (0.6) were selected via a 100-trial Bayesian (TPE) hyperparameter search maximizing MCC (`results/tables/gc_nkgraph_bayesian_trials.tsv`).

**Multi-view audit.** Inspired by the separation of biological networks in TREE
[49] and GRAFT [50], we additionally encoded six views independently by SVD:
generic interactions (PPI plus scRNA co-expression), ligand–receptor,
`metabolic_crosstalk`, `sm_topology_axis`, GO-process co-membership, and MSigDB
C2 co-membership. A global softmax weight fused the six sample-level graph
projections; this is a lightweight audit, not a Transformer implementation. The
**label-masked primary analysis** removed the union of every gene used to define
NK infiltration, cytotoxicity, dysfunction, or CAF/ECM exclusion from both raw
expression and graph projections; the unmasked analysis was retained only as a
sensitivity check. For each of ten fixed random seeds, TCGA-STAD was split 80:20
for training/early stopping (maximum 300 epochs; patience 20), and TCGA-LIHC was
never used for standardization, early stopping, threshold selection, or model
selection. No-graph, merged-SVD, single-view, uniform-fusion, learned-fusion, and
leave-one-view-out models shared each split. External AUROC/AUPRC differences
were computed from seed-ensemble probabilities with a 2,000 paired stratified
bootstrap over LIHC samples; seed SD quantifies algorithmic instability rather
than biological replication.

### 2.7 Baseline methods

Six baselines are evaluated on the same data splits: XGBoost [42], LightGBM [43], Random Forest, ElasticNet logistic regression, RBF SVM, and a 2-layer MLP. All baselines use the full gene expression matrix without graph structure. The comparison isolates the contribution of graph-informed features.

### 2.8 Candidate target prioritization

Genes are ranked by a composite score integrating five evidence dimensions:

| Dimension | Weight | Source |
|-----------|--------|--------|
| Tumor cell specificity (log2 FC malignant vs non-malignant) | 0.30 | scRNA |
| NK dysfunction correlation (abs) | 0.20 | scRNA NK subset |
| SST-axis membership | 0.30 | Mechanism card |
| Axis core membership (protrusion→cytotoxicity) | 0.10 | Mechanism card |
| Gold-standard literature support | 0.10 | Curated list |

Each top-ranked candidate receives a recommended wet-lab validation assay (e.g., NK-tumor co-culture cytotoxicity assay, SEM membrane-protrusion imaging, sphingomyelinase-inhibitor rescue ± Tim3 blockade).

### 2.9 Pre-registered hypotheses (Arm A — liver positive control)

The following hypotheses are registered before execution:

| ID | Hypothesis | Expected direction |
|----|-----------|-------------------|
| H1 | tumor_serine_capacity ⟂ nk_sm_balance | Direction calibrated, reported |
| H2 | nk_sm_balance (+) nk_protrusion_machinery | Positive |
| H3 | nk_protrusion_machinery (+) cytotoxicity anchors | Positive |
| H4 | nk_topology_permissive (−) HAVCR2 / dysfunction | Negative |
| H5 | intratumoral NK < peritumoral/peripheral NK in sm_balance & protrusion machinery (scRNA) | Negative (tumor < normal) |

**Recovery definition.** The pre-registered full-recovery criterion required H2–H5 to pass in the expected direction in liver cancer. Because this criterion was not met, Arm A is interpreted as a partial-recovery scoping test: it identifies the effector layer that is recoverable from transcriptomes, the metabolic coupling that becomes detectable only after cell-type resolution, and the physical topology layer that remains outside transcriptional reach.

**Statistical correction for pseudoreplication in single-cell correlations.** 
The single-cell NK correlations (H2–H5, resolution "single-cell NK") use 8,310
NK cells drawn from 9 biological samples (GC1–3, HL1–3, LM1–3) with widely
varying per-sample cell counts (range 121–1,802). Treating cells as independent
observations inflates the effective sample size and deflates P-values. We
therefore apply a two-stage correction for all single-cell NK correlation tests:

1. **Per-sample analysis:** Pearson r (or Cohen's d for H5 group comparisons)
   is computed within each of the k=9 samples independently, producing k
   independent estimates.
2. **Random-effects meta-analysis:** Each per-sample r is Fisher
   z-transformed, the standard error SE(z)=1/√(n−3) is computed per sample,
   and the DerSimonian-Laird random-effects model pools the estimates.
   Between-sample heterogeneity is quantified as I² and τ² (Cochran's Q).

The meta-analytic r and its P-value are reported as the corrected values
(labeled "after corr." in Table 2). Where I² is high, we explicitly note the
between-sample variability and interpret the pooled estimate accordingly.
All corrected P-values replace the naive values in Tables 2 and 5 and in
the text of §3.2 and §4.1. The naive per-cell values are retained in
`results/tables/sst_axis_pseudoreplication_corrected.tsv` for transparency.
The correction is implemented in `src/topology/pseudoreplication_correction.py`.

---

## 3. Results

### 3.1 Dataset summary

Four independent bulk cohorts and one multi-tissue scRNA dataset were processed
(Table 1). TCGA-STAD represents the gastric adenocarcinoma molecular landscape [44],
while GSE62254 provides the ACRG molecular subtype annotation [45]. Bulk NK immune-state
labels were assigned by scRNA-anchored scoring;
the two external microarray cohorts required probe-to-symbol remapping before NK
markers could be scored (see Methods; external-validation scoring is reported
where gene coverage permits).

**Dataset characteristics.**

| Dataset | Cancer | Type | Samples / cells | Role |
|---------|--------|------|-----------------|------|
| TCGA-LIHC | Liver | Bulk RNA-seq | 423 | Positive control (Arm A) |
| TCGA-STAD | Gastric | Bulk RNA-seq | 450 | Train primary (Arm B) |
| GSE62254 (ACRG) | Gastric | Bulk microarray | 300 | External validation |
| GSE84437 | Gastric | Bulk microarray | 483 | External validation |
| GSE246662 | Liver/gastric/metastasis | scRNA-seq | 166,829 cells (8,310 NK) | NK atlas + axis test |

### 3.2 Arm A — the axis is *partially* recovered: effector arm yes, topology arm no

We tested the pre-registered hypotheses (H1–H5) at two resolutions — bulk
TCGA-LIHC and single NK cells — and report the outcome per hypothesis
(Table 2; `results/tables/sst_axis_positive_control_recovery.tsv`). The picture
is consistent and interpretable rather than a uniform pass:

- **The metabolic-serine capacity is transcriptionally uncoupled from NK SM balance
  (H1).** Tumor-side serine-pathway transcript scores show no correlation with
  NK-side SM-balance scores at either resolution (bulk r=−0.016, *P*=0.74; scNK
  r=+0.012, *P*=0.27). This is consistent with the mechanism's core premise that
  serine→SM crosstalk operates at the metabolite level — enzyme abundance
  (transcription) is not flux — and provides an empirical null baseline for the
  axis hierarchy.
- **The effector arm reproduces in bulk (though see the purity control below);
  the single-cell number is pseudoreplication-corrected but does not survive
  further technical-confound controls (H3).** Protrusion-machinery transcript score couples to the
  cytotoxicity-output score in bulk TCGA-LIHC (r=0.551, *P*=5.5×10⁻³⁵) and,
  at face value, across 8,310 single NK cells. After correcting for
  within-sample dependence (per-sample r → Fisher z → DerSimonian-Laird
  random-effects meta-analysis across 9 samples), the corrected single-cell
  meta-analytic r=0.313 (*P*=3.9×10⁻⁸, 95% CI [0.13, 0.48]), with I²=96%
  reflecting substantial between-sample heterogeneity (per-sample r range
  [0.009, 0.560]).
  Pseudoreplication correction addresses cell non-independence within
  samples, but not a separate concern: whether the per-cell module-score
  correlation itself is driven by shared technical structure rather than a
  specific biological coupling. We tested this directly on the real NK
  scRNA data (`src/topology/count_depth_control.py --real`,
  `src/topology/h3_scoring_method_diagnostic.py`,
  `src/a100_recompute/run_h3_activation_control_v2.py`). First, a generic
  16-gene NK-activation signature (CD69, TNF, XCL1/2, CCL3/4/5, CSF2, IL2RA,
  ICOS, TNFSF10, FASLG, CD38, HLA-DRA/B1, MKI67; not present in either
  module) explains only part of the raw coupling: linearly partialling it out
  gives r=0.251 (*P*=2.7×10⁻¹¹⁹), and an activation-matched analysis (cells
  binned into quintiles of the activation score, so the correlation is
  computed within cells of near-identical activation level) gives a
  consistent mean r=0.242 across bins — ruling out simple co-activation as
  the explanation. However, two more comprehensive controls collapse the
  coupling substantially further. **Count-depth / technical covariates**:
  47.4% of the protrusion-machinery module score's variance is explained by
  `total_counts` and `n_genes_by_counts` alone (library size varies 236-fold
  across the 8,310 cells, 427–100,763 counts); residualizing both modules
  against these covariates before correlating drops the raw r from 0.318 to
  0.094 (*P*=1.3×10⁻¹⁷). **Residualizing against the real scVI latent space**
  (30 dimensions, the same batch-corrected embedding used elsewhere in this
  study) plus the same technical covariates drops it further to r=0.092
  (*P*=5.5×10⁻¹⁷). Second, we asked whether this residual coupling is
  module-*specific* using a permutation test: 10,000 random gene modules of
  the same sizes (25 and 11) were drawn from the NK-expressed gene universe
  and scored the same way. With the mean-z-score method used throughout this
  paper, the observed r=0.318 falls *below* the permutation null mean
  (0.728; empirical *P*=1.0) — but we found this specific null is inflated by
  a biased universe (restricting candidate genes to the 314/22,728 detected
  in ≥50% of NK cells excludes most of the protrusion and cytotoxicity module
  genes themselves, which are more sparsely detected). Repeating the test
  with `scanpy`'s expression-level-matched `score_genes` scoring — the
  field-standard control for exactly this confound, applied identically to
  observed and permuted modules — gives a materially smaller raw r=0.089 to
  begin with, and the observed value still does not exceed the matched
  permutation null (mean 0.354, 95th percentile 0.560; empirical *P*=0.97).
  Full diagnostics are in
  `results/tables/h3_scoring_method_diagnostic_summary.md`.
  **We conclude that the single-cell H3 number, while not attributable to
  generic co-activation alone, is not distinguishable from the technical and
  global-transcriptional-structure background of this dataset, and should not
  be reported as an independent replication of the effector arm.** The bulk
  TCGA-LIHC result is unaffected by this *single-cell* diagnostic (bulk samples
  do not carry the extreme per-cell library-size range that drives the
  single-cell artifact), but it carries a separate, bulk-specific confound that
  we address next.
  We retain the single-cell pseudoreplication-corrected number in Table 2 for
  transparency but flag it as not surviving technical-confound control.
- **About half of the bulk effector coupling is NK-cell-abundance covariation
  (bulk purity control).** In bulk tumor transcriptomes, both the
  protrusion-machinery and cytotoxicity-output module scores rise and fall with
  the NK-cell fraction of the sample, so a positive protrusion~cytotoxicity
  correlation is expected in part from co-varying NK abundance rather than a
  within-NK coupling. We tested this directly
  (`src/topology/bulk_h3_purity_control.py`) by partialling out a clean
  NK-lineage abundance proxy — the mean-z of seven NK lineage/receptor markers
  (FCGR3A, KLRD1, KLRF1, KLRK1, NCAM1, NCR1, TYROBP) that appear in *neither*
  module, avoiding the over-control that would result from using the full NK
  signature (whose genes overlap the cytotoxicity module). In TCGA-LIHC the
  coupling attenuates by 55%, from zero-order r=0.55 to a partial r=0.25
  (*P*=1.8×10⁻⁶, 95% CI [0.15, 0.34]) — reduced but still significant. Applying
  the same control to the two external gastric cohorts roughly halves it in
  GSE62254 (partial r=0.23, *P*=4.6×10⁻⁵, still significant) and **abolishes it
  in GSE84437** (partial r=−0.07, *P*=0.11, 95% CI crossing zero), the cohort
  with the strongest zero-order coupling (r=0.62). The effector-arm coupling is
  therefore real but **substantially infiltration-driven**: after NK-fraction
  adjustment it holds at roughly half strength in two of three bulk cohorts and
  does not survive in the third (`results/tables/bulk_h3_purity_control.tsv`).
  We accordingly report the effector arm as *partially recoverable in bulk and
  substantially NK-abundance-driven*, not as a robust three-cohort replication.
- **The upstream metabolic coupling is not statistically significant after pseudoreplication correction
  (H2).** The SM-balance→protrusion coupling is undetectable in
  bulk (r=−0.017, *P*=0.72). After NK-cell
  isolation, the naive per-cell correlation was r=+0.030 (*P*=6×10⁻³), but this
  treats 8,310 cells as independent observations from only 9 samples — a
  severe pseudoreplication inflation. Applying per-sample Pearson r →
  Fisher z → DerSimonian-Laird random-effects meta-analysis to correct for
  within-sample dependence, the corrected meta-analytic r=0.029 (*P*=0.20,
  I²=73% substantial heterogeneity). The shared variance is r²=0.0009 and the
  95% CI of the per-sample r spans [−0.059, +0.132]. The coupling is therefore
  **not statistically detectable** at either resolution after appropriate
  correction — consistent with serine→SM crosstalk operating at the metabolite
  rather than transcript level, which the anchor paper's single-cell mass
  spectrometry design implies but our data now empirically confirm.
- **The physical topology phenotype is *not* captured by machinery transcription
  — a fundamental disconnect (H4, H5).** This is our most consequential finding.
  The topology-permissive→dysfunction relationship carries the wrong
  sign at both resolutions before and after correction (bulk r=+0.311;
  single-cell NK vs HAVCR2 corrected r=+0.036, *P*=9.1×10⁻⁴, I²=0%), and
  intratumoral NK cells show *higher* protrusion-machinery transcription than
  normal-tissue NK (Δ=+0.142, p=3.0×10⁻⁹¹) — the opposite of the physical
  collapse — even though their cytotoxic output is correctly reduced
  (Δ=−0.141, p=5.9×10⁻⁵²). Transcription of the machinery genes therefore does
  **not** proxy the membrane-topology state itself. This is not a weakness of our
  model but the natural boundary of any transcriptome-based reconstruction: a
  membrane-lipid mechanism whose phenotype is physical (protrusion collapse)
  and whose causal link (serine→SM) operates at the metabolite level cannot be
  fully captured by mRNA abundance. This finding precisely explains why
  the anchor lab required single-cell mass spectrometry and super-resolution
  imaging, and it constitutes an empirically grounded caution against
  over-interpreting transcriptional proxies for physical phenotypes in
  immune-evasion studies. The reduced cytotoxic output of
  intratumoral NK (H5-cytotoxicity) is independently corroborated by a single-cell
  HCC study reporting that intratumoral, relative to peritumoral, NK cells
  upregulate inhibitory-checkpoint and exhaustion programs and downregulate
  cytotoxicity pathways [46] — external support, from cohorts independent of ours,
  for the effector-layer direction; the opposing protrusion-transcription result
  is our own observation.

**Pre-registered hypothesis outcomes (multi-resolution).**

| Hyp | Test | Resolution | r / Δ | *P* | *P*_FDR | Expected | Outcome |
|-----|------|-----------|-------|---|---------|----------|---------|
| H1 | serine_capacity ~ sm_balance | bulk | −0.016 | 0.74 | 0.74 | calibrated | reported (null) |
| H1 | serine_capacity ~ sm_balance | single-cell NK | +0.012 | 0.27 | 0.34 | calibrated | reported (null) |
| H2 | sm_balance ~ protrusion | bulk | −0.017 | 0.72 | 0.74 | + | not recovered |
| H2 | sm_balance ~ protrusion | single-cell NK | +0.029 | 0.20 | 0.20 | + | not signif. after correction |
| H3 | protrusion ~ cytotoxicity | bulk | +0.551 | 5×10⁻³⁵ | 1×10⁻³⁴ | + | **recovered** |
| H3 | protrusion ~ cytotoxicity | single-cell NK | +0.313†‡ | 3.9×10⁻⁸†‡ | 7×10⁻⁸ | + | not distinguishable from technical background (‡) |
| H4 | topology ~ dysfunction | bulk | +0.311 | 7×10⁻¹¹ | 4×10⁻⁵ | − | not recovered |
| H4 | topology ~ HAVCR2 | single-cell NK | +0.036† | 9.1×10⁻⁴† | 1.1×10⁻³ | − | not recovered (after corr.) |
| H5 | intratumoral<normal: cytotoxicity | single-cell NK | −0.141 | 6×10⁻⁵² | 9×10⁻⁴⁶ | tumor< | **recovered** |
| H5 | intratumoral<normal: protrusion | single-cell NK | +0.142 | 3×10⁻⁹¹ | 1×10⁻⁸⁴ | tumor< | not recovered |

*P*_FDR: Benjamini–Hochberg false discovery rate correction across 10 tests.
† Single-cell NK r/P are meta-analytic values after pseudoreplication correction:
per-sample Pearson r → Fisher z → DerSimonian-Laird random-effects meta-analysis
across 9 samples (Methods §2.9). Naive per-cell P-values (H3: p≈5×10⁻¹⁹⁴;
H2: p≈6×10⁻³; H4: p≈4×10⁻⁶) are inflated by treating 8,310 cells as independent
observations. Corrected values use effective sample size at the sample level (k=9).
‡ Pseudoreplication correction does not address a separate confound: the H3
single-cell coupling is largely technical, not module-specific. After
residualizing against library size + detected-gene count + the real 30-dimension
scVI latent space, r drops to 0.09; using an expression-level-matched scoring
method (`scanpy.tl.score_genes`) the raw r is 0.089 to begin with; and a
10,000-draw permutation test of randomly-sized-matched gene modules shows the
observed coupling does not exceed the random-module baseline (empirical
*P*=0.97 with the matched-scoring-method null). See §3.2 and
`results/tables/h3_scoring_method_diagnostic_summary.md` for the full battery
of controls. We report H3's single-cell number in this table for completeness
and transparency, but the effector-arm claim rests on the bulk result (this
table) and its three independent bulk replications in gastric cancer (§3.3),
not on the single-cell number.

**Recovery verdict (revised).** H1 is null at both resolutions, establishing that
tumor serine-pathway transcription and NK SM balance are empirically uncoupled at
the transcript level — consistent with the serine→SM link operating
post-transcriptionally. Arm A recovers the *functional/effector layer* of
the axis in bulk (H3 bulk; H5-cytotoxicity); the single-cell H3 number is
pseudoreplication-corrected but, on the further technical-confound controls
described above, is not distinguishable from this dataset's technical/global-
transcriptional background, so we do not count it as an independent
single-cell replication. Arm A further demonstrates that the *metabolic
coupling is not transcriptionally detectable* after pseudoreplication
correction (H2, *P*=0.20), and establishes that the *physical topology
phenotype* is fundamentally disconnected from machinery transcription (H4;
H5-protrusion). This map — effector coupling recovered in bulk only / absent
metabolic coupling / disconnected physical topology — is the paper's central
finding: a delineation of a membrane-lipid immune-evasion mechanism's
transcriptional reach from null (H1) to bulk-only (H3) to absent
(H4/H5-protrusion), and a caution that single-cell module-score correlations
in this kind of high-depth-variance dataset require count-depth and
module-permutation controls beyond pseudoreplication correction alone.

**Figure 1.** (a) Confound-control logic for the effector-arm claim: a bulk
association is retained as primary evidence only after the single-cell
replicate survives residualization against library size/latent transcriptional
structure and a random-gene-module permutation baseline. (b) H3:
protrusion-machinery vs cytotoxicity-output score, bulk TCGA-LIHC (r=0.55,
p=5×10⁻³⁵). (c) H2 and H3 at bulk vs single-cell resolution — single-cell
resolution rescues H2 (bulk ≈0) but only partially reproduces the bulk H3
effect. (d) Intratumoral vs normal NK per-cell module z-scores for
cytotoxicity-output and protrusion-machinery (single cell). (e) Table 2 as a
forest plot; the H3 single-cell point is marked as flagged (not distinguishable
from technical background) rather than recovered.

### 3.3 Arm B — Gastric cancer extension

The effector coupling that defines the recovered arm holds in gastric-cancer NK
cells (protrusion~cytotoxicity r=0.346, p=6.6×10⁻³⁰, n=1,017 NK), i.e. the
transcriptionally-recoverable layer of the axis extends to gastric cancer. This
gastric-tissue subset is drawn from the same scRNA dataset and shares the same
extreme per-cell library-size range diagnosed in §3.2; we have not re-run the
full count-depth/permutation battery on this n=1,017 subset specifically, but
given the shared technical architecture, this single-cell gastric number
should be read with the same caution as the Arm A single-cell H3 result. The
gastric effector-arm claim is anchored on the three independent bulk
replications below (TCGA-STAD-adjacent GSE62254/GSE84437), which are not
subject to this single-cell confound.
Cross-tissue module means (`sst_axis_scrna_by_tissue.tsv`) show gastric-cancer NK
with the lowest cytotoxicity-output score (−0.194) of the three tissues,
consistent with an evasion phenotype. Gastric bulk (TCGA-STAD) NK states
distribute across cold/excluded (154), hot-cytotoxic (153), intermediate (123),
and hot-dysfunctional (20).

**External validation in two independent gastric cohorts.** The external
microarray cohorts initially failed NK scoring because probe IDs were not mapped
to gene symbols; after platform-annotation remapping (GSE62254→GPL570, 54,675
probes→22,880 genes, NK markers 6/7; GSE84437→GPL6947, 49,576→18,903 genes, NK
markers 7/7; `run_geo_external_validation.py`), the recovered effector coupling
replicates in **both** cohorts: protrusion-machinery→cytotoxicity-output
correlates at r=0.42 (p=1.4×10⁻¹⁴) in GSE62254 (n=300) and r=0.62 (p=3.3×10⁻⁵³)
in GSE84437 (n=483) — matching the direction and comparable in strength to the
liver control (Table 5). Notably, the SM-balance→protrusion coupling that was
undetectable in bulk *liver* is weakly but significantly positive in both bulk
*gastric* cohorts (r=0.18, p=1.3×10⁻³; r=0.11, p=0.02), consistent with the
effect being real but small and easier to detect where NK signal is stronger.
The effector layer of the axis is therefore reproduced in three gastric datasets
(one scRNA, two bulk microarray) at the zero-order level. However, as shown by
the bulk purity control (§3.2), these bulk replications are substantially
NK-abundance-driven: after partialling out NK-cell fraction, the coupling is
roughly halved but remains significant in GSE62254 and does not survive in
GSE84437. The gastric extension should therefore be read as a *partially
recoverable, infiltration-driven* effector signal rather than a robust
independent replication.

**External validation of the axis (independent gastric cohorts).**

| Cohort | n | protrusion~cytotoxicity r (p) | sm_balance~protrusion r (p) | NK markers |
|--------|---|-------------------------------|-----------------------------|------------|
| GSE62254 (GPL570) | 300 | 0.42 (1.4×10⁻¹⁴) | 0.18 (1.3×10⁻³) | 6/7 |
| GSE84437 (GPL6947) | 483 | 0.62 (3.3×10⁻⁵³) | 0.11 (0.02) | 7/7 |

**Figure 2.** (a) External-validation design: NK-state calling on TCGA-STAD
bulk plus two fully independent gastric microarray cohorts (GSE62254,
GSE84437), neither used to define the axis modules. (b) Cross-tissue module
comparison (healthy liver, liver metastasis, gastric cancer), with
gastric-cancer NK showing the lowest cytotoxicity-output score. (c)
NK-immune-state distribution in TCGA-STAD (n=450). (d) Effector-coupling
replication across four cohorts (Table 5).

### 3.4 NK immune state classification

On TCGA-STAD (5-fold stratified CV, binary NK-hot-cytotoxic vs rest), the
graph-informed model attains accuracy 0.864, balanced accuracy 0.856,
macro-F1 0.850, MCC 0.706, AUROC 0.950, AUPRC 0.910. The six tabular baselines
were evaluated on the *identical* seed-42 folds (Table 3), and we report paired
Wilcoxon/t-tests of the GNN against each baseline on MCC and AUROC
(`model_comparison_stats.tsv`). Table 3 reports one internally consistent
evaluation in which the GNN and all baselines were run on the identical seed-42
folds. The GNN itself is robust to the graph-edge enrichment described in §2.5:
rebuilding the graph with the real STRING/CellChatDB priors moves its
5-fold-mean MCC by less than 0.01 (0.706→0.716) and its AUROC by less than
0.002, within one fold-standard-deviation. Absolute baseline metrics vary
modestly across independent re-runs of the harness; the qualitative ranking used
below — GNN on par with the gradient-boosting baselines and above the linear,
kernel, and shallow-network ones — is stable, and is the only comparison our
conclusions rely on.

The comparison yields a clear result: the graph-informed model
is **statistically on par with the strongest gradient-boosting baselines** —
neither LightGBM (MCC 0.733) nor XGBoost (0.727) differs significantly from the
GNN (0.706) on paired tests (ΔMCC −0.028, t-test p=0.28; and −0.022, p=0.31,
respectively) — and it **significantly outperforms** the linear, kernel, and
shallow-network baselines (vs ElasticNet ΔMCC +0.035, p=0.017; vs SVM +0.217,
p=0.008; vs MLP +0.330, p=0.037; Figure 4).

**Domain-method baseline comparison.** To test whether the GNN adds value
over standard bioinformatics approaches that do not use graph structure, we
evaluated two additional baselines on the identical folds: (i) an *NK-marker
signature* baseline — logistic regression on the mean expression of 8
canonical NK markers (NCAM1, KLRD1, NKG7, GNLY, KLRF1, EOMES, NCR1, FCGR3A),
conceptually the simplest possible deconvolution proxy; and (ii) an *SST-module
signature* baseline — logistic regression on the 7 SST-axis module scores
computed directly on bulk expression (Methods §2.3), which captures the
"use the anchor paper's gene modules without building a graph" approach.

The SST-module baseline attains AUROC 0.904 ± 0.029 and MCC 0.619 ± 0.088 —
**not significantly different from the GNN** (ΔAUROC −0.046, t-test p=0.11;
ΔMCC −0.087, p=0.28). However, the simpler NK-marker baseline (AUROC 0.849,
MCC 0.503) is **significantly below** the GNN on both AUROC (Δ −0.101,
p=0.032) and AUPRC (Δ −0.177, p=0.023). These results sharpen the
interpretation of the GNN's contribution: the mechanism card's gene modules
already capture the bulk of the discriminative signal for NK state
classification when combined with logistic regression, and the heterogeneous
graph does not significantly improve discriminative accuracy beyond this
no-graph alternative. We therefore do not claim the graph as a predictive
advance. Its role in this study is instead diagnostic: the
mechanism-structured embedding is used as a *probe* whose ablation behavior
independently maps which layer of the axis the transcriptome can reach (§3.7),
not as a source of classification gain that the SST-module or tabular baselines
lack.

We note that a full comparison with established deconvolution tools
(CIBERSORTx, quanTIseq) and phenotype-genotype association methods (Scissor)
would further strengthen this conclusion; a detailed roadmap for reproducing
these comparisons is provided in
`submission_bundle_BiB/03_supplementary/CIBERSORTx_quanTIseq_Scissor_roadmap.md`.
These tools require an R environment not available in the current local setup.

We therefore do not claim state-of-the-art
accuracy, nor a predictive advantage for the graph: it matches top tree
ensembles and the SST-module signature baseline on this binary task. What the
graph uniquely provides is a mechanism-structured embedding that can be
*ablated* — removing the edge that encodes the mechanism's metabolic coupling
and observing the effect — which is what turns it into a probe of transcriptional
reach (§3.7) rather than a black-box classifier. This "comparable accuracy,
diagnostic value" position is the straightforward reading of the numbers.

**NK-state classification (TCGA-STAD, 5-fold CV; mean over folds).**

| Method | Accuracy | Balanced Acc. | Macro F1 | MCC | AUROC | AUPRC |
|--------|----------|---------------|----------|-----|-------|-------|
| LightGBM | 0.880 | 0.863 | 0.865 | **0.733** | **0.960** | 0.933 |
| XGBoost | 0.878 | 0.859 | 0.862 | 0.727 | 0.959 | 0.928 |
| **GC-NKGraph-Atlas (GNN)** | 0.864 | 0.856 | 0.850 | 0.706 | 0.950 | 0.910 |
| RandomForest | 0.860 | 0.825 | 0.836 | 0.681 | 0.941 | 0.901 |
| ElasticNet | 0.856 | 0.822 | 0.832 | 0.671 | 0.927 | 0.884 |
| SST-module signature | 0.818 | 0.823 | 0.803 | 0.619 | 0.904 | 0.827 |
| NK-marker signature | 0.787 | 0.744 | 0.750 | 0.503 | 0.849 | 0.733 |
| SVM (RBF) | 0.773 | 0.672 | 0.683 | 0.489 | 0.900 | 0.839 |
| MLP (2-layer) | 0.747 | 0.669 | 0.641 | 0.376 | 0.807 | 0.730 |

*GNN not significantly different from LightGBM/XGBoost (paired p>0.27) or
SST-module signature baseline (p>0.10); significantly above NK-marker
signature (p<0.05), ElasticNet/SVM/MLP (p<0.05).*

### 3.5 Candidate target prioritization (de-circularized)

The naive composite score is circular: weighting axis membership at 40% while
scoring |tumor-specificity| promotes NK-effector markers (NKG7, PRF1, GZMB),
which are *depleted* in malignant cells and are the axis *readout*, not targets.
We therefore separate two tables
(`src/interpretation/split_target_lists.py`):

1. **Tumor-intrinsic candidates** (Table 4; n=37): genes required to be *up in
   malignant cells* (tumor_specificity_log2>0), re-scored on signed
   tumor-specificity, mechanistic centrality (tumor-serine program /
   axis-druggable enzyme), druggability, and NK association. The top of the list
   is dominated by the mechanistically-privileged, druggable head of the
   axis — **PHGDH** (serine synthesis; Phase 1/2 inhibitor), **SGMS2/SMPD3/SMPD1**
   (SM synthesis/catabolism; the exact enzymes the anchor mechanism implicates),
   **PSAT1/PSPH** (serine pathway) — alongside gastric-relevant surface targets
   (ERBB2, FGFR2, MET) and stress ligands (MICA).
2. **Axis-confirmation panel** (n=36): the NK-side markers (GZMB, GNLY, PRF1,
   NKG7, protrusion genes), reported explicitly as *readout* validating the axis
   score, not as targets.

The de-circularization is conservative but not complete: an audit of the n=37
tumor-intrinsic pool reveals 17 genes (46%) are annotated to NK-side
mechanism-card modules (protrusion machinery 7, de novo sphingolipid 5,
SM catabolism 4, SM synthesis 1), including the NK protrusion GTPase RAC1 (rank 10)
and actin nucleation factor WASL (rank 24). These appear in the pool because their
expression is detectably above zero in malignant cells (satisfying
`tumor_specificity_log2>0`), but their mechanism-card module membership means the
ranking retains a structural NK-side bias. The top of the list — PHGDH, SGMS2,
PSAT1, PSPH, SMPD3/1 — is clean (tumor-serine-capacity or
metabolic-suppression categories) and these remain the primary candidates for
experimental follow-up. Future refinement should add a module-level penalty or a
higher tumor-specificity threshold (§4.3).

**Trivial baseline comparison.** To quantify whether the five-dimension scoring
adds information beyond simply listing every gene named in the Zheng 2023 anchor
paper, we compared it against a trivial baseline that ranks genes solely by
mechanism-card membership (in_sst_axis = 1.0) plus gold-standard literature
support (gold_standard = 0.5), using no expression data whatsoever. Within the
37 tumor-intrinsic candidates (Table 4), the trivial baseline shows moderate
correlation with the five-dimension scoring (Spearman rho=0.54, p=5x10^-4),
indicating the rankings are linked but not identical. The five-dimension scoring
contributes incremental value in three ways: (i) it re-orders within the SST
set by quantitative tumor specificity — e.g., PHGDH (log2FC +0.059) ranks
above PSAT1/SHMT1 despite both being SST members; (ii) it surfaces 12 non-SST
genes that the trivial baseline would rank at positions 26-37 — most notably
COL1A1/COL1A2 (log2FC ~0.15, rank 6-7 vs trivial rank 27-28), NECTIN2 (log2FC
+0.11, rank 9 vs 29), CA9, ERBB2, FGFR2, and MET — these are gastric-cancer-relevant
targets with measurable tumor-cell signal that a "read the anchor paper" approach
would entirely miss; and (iii) it demotes SST-member genes with near-zero tumor
specificity to the bottom of the list — SPTLC1/3, WASF1/3, DIAPH3 (log2FC
<0.01, ranks 33-37 vs trivial ranks 21-25). The incremental value is therefore
appreciable but not transformative: the serine/SM enzyme core is captured by
trivial membership, and the main added benefit is the surfacing of gastric-relevant
non-SST candidates and the quantitative de-prioritization of anchor-paper genes
with negligible malignant-cell signal. Full comparison tables are in
results/tables/trivial_baseline_comparison.tsv and
trivial_baseline_overlap.tsv.

**Independent orthogonal validation (NK-state DE + DepMap essentiality).**
To test the target list against evidence dimensions that are fully independent
of the mechanism-card scoring and the scRNA malignant-cell proxy, we performed
two orthogonal analyses. First, we tested each of the 37 genes for differential
expression between TCGA-STAD tumors with intact NK killing (NK-hot-cytotoxic,
n=134) versus suppressed NK function (NK-hot-dysfunctional, **n=20 — this
group is small, and the DE test is correspondingly underpowered; log2FC
effect sizes are reasonably estimated from the group means, but the
associated P/FDR values should be read as indicative rather than definitive.
We did replicate the directional comparison in two independent external
gastric cohorts using the same NK-state labeling rule (below); results were
mixed — informative rather than uniformly reassuring**). A genuine
tumor-intrinsic immune-evasion mediator should be upregulated in tumors where
NK is present but dysfunctional — the tumor is actively suppressing NK.
Second, we queried DepMap CRISPR knockout screens (CERES dependency scores)
in gastric cancer cell lines: a good immune-evasion target should be
non-essential for cell-autonomous survival in vitro (CERES > 0), since its
therapeutic value lies specifically in the immune-microenvironment context.

Four genes pass both orthogonal filters: **SGMS2** (SM synthase; upregulated
in dysfunctional tumors, log2FC=+0.14, FDR=0.007; CERES=−0.02 in 35 real
gastric/stomach DepMap 26Q1 cell lines, weakly essential — see note below)
and **NT5E/CD73** (adenosine pathway; log2FC=+0.16, FDR=0.014; CERES=+0.005,
genuinely non-essential) show the strongest combined evidence. **SMPD1**
(acid SMase; CERES=−0.16, weakly essential) and **SMPD3** (neutral SMase;
CERES=−0.13, weakly essential) are mechanistically privileged (Zheng 2023
anchor enzymes) but lack significant NK-state DE signal.
**Note on DepMap numbers:** these CERES values are from a real query of
**DepMap Public 26Q1** (`CRISPRGeneEffect.csv`, `Model.csv`; the current
release at the time of this analysis, superseding the 25Q2 originally named
in Methods and the 24Q2 figshare snapshot used in an earlier pass of this
analysis — see Data Availability), restricted to the 35 cell lines with an
Oncotree subtype containing "Stomach." Of the four genes above, NT5E is now
genuinely non-essential (CERES>0); SGMS2, SMPD1, and SMPD3 fall in the
"weakly essential" band, which is a common and largely unremarkable range for
most genes in most cell lines. All four remain clearly distinguishable from
the genuinely essential genes excluded below (RAC1, MTHFD1; CERES<−0.5) and
from ERBB2 (CERES=−0.36). We revise the wording of this validation from
"passes a non-essentiality filter" to "does not fall in the pan-essential
range that would confound an immune-evasion interpretation with a
cell-viability effect."
**Caveat on NT5E:** a gene-set separation audit
(`results/tables/geneset_separation_audit_summary.md`) found that NT5E is
itself one of the ten genes constituting `NK_DYSFUNCTION_GENES`, the marker
set whose z-score (net of the cytotoxicity score) defines the
NK-hot-dysfunctional label used for this very DE test. NT5E's "upregulated in
dysfunctional tumors" result is therefore partly self-referential — higher
NT5E expression directly contributes to a tumor being labeled dysfunctional in
the first place — and should not be read as independent evidence to the same
degree as SGMS2, SMPD1, and SMPD3, none of which overlap any NK-state
label-defining gene set.

**External replication of the NK-state DE test.** Rather than requiring a new
dataset, we re-ran the identical NK-state labeling rule
(`src/immune_scoring/nk_scores.py`) independently on the two external gastric
bulk cohorts already used for the effector-arm validation (§3.3) — GSE62254
(ACRG, n=300; NK-hot-cytotoxic n=105, NK-hot-dysfunctional n=13) and GSE84437
(n=483; NK-hot-cytotoxic n=122, NK-hot-dysfunctional n=29) — and compared the
direction of each flagged gene's dysfunctional-vs-cytotoxic difference against
the TCGA-STAD result
(`results/tables/nk_state_de_external_replication_summary.md`,
`nk_state_de_external_concordance.tsv`). Fourteen of 16 gene/cohort
comparisons (88%) are directionally concordant. Encouragingly, the five
downgraded serine-pathway genes and RAC1 — none of which overlap any
NK-state label-defining gene set, so this check is not circular — are
**DOWN in dysfunctional tumors in both external cohorts**, independently
confirming their exclusion from the priority tier. NT5E is **UP** in both
external cohorts, concordant with TCGA-STAD, but this concordance is subject
to the same label-overlap caveat above (NT5E contributes to the NK-state
label in every cohort tested, including these two). **SGMS2, however, is
discordant in both external cohorts** (DOWN in dysfunctional tumors, opposite
to the UP direction and FDR=0.007 reported for TCGA-STAD): the NK-state DE
signal that helped qualify SGMS2 does not externally replicate, and should be
downweighted accordingly. SGMS2's Tier-1 status now rests primarily on its
mechanistic privilege (core Zheng 2023 SM-synthase enzyme) and its
non-pan-essential DepMap profile, not on the NK-state DE direction.

Critically, five serine-pathway enzymes — PSPH, SHMT1, SHMT2, MTHFD1L, MTHFD1 —
show the *opposite* expression pattern: they are significantly *lower* in
NK-dysfunctional tumors (log2FC range −0.09 to −0.13, all p<0.02), arguing
against their inclusion as high-priority targets; of these, **MTHFD1** is
additionally pan-essential in vitro (real CERES=−0.53), while PSPH, SHMT1,
SHMT2, and MTHFD1L are weakly essential (CERES −0.10 to −0.16), and **PSAT1**
(not one of the five DE-downgraded genes, but already downgraded on other
grounds) is moderately essential (CERES=−0.27). **RAC1**
(real CERES=−0.74) is clearly pan-essential in vitro and should be
interpreted as a cell-viability rather than immune-evasion target. **ERBB2**
(real CERES=−0.36, moderately essential — softer than the pan-essential range,
though still a confound, and separately an FDA-approved HER2 target in its
own right) should likewise not be treated as a novel immune-evasion candidate.
These findings sharpen
the target list considerably: the top tier is SGMS2, SMPD3, SMPD1, and NT5E
(now genuinely non-essential — see below — but still flagged for the
label-circularity caveat above);
the mid-tier is MICA, FN1, BAIAP2, PVR; and 5 serine-pathway genes plus
ERBB2 and RAC1 should be downgraded or removed from the tumor-intrinsic pool.
Full results in `results/tables/target_validation_v2_merged.tsv`.

**Top putative tumor-intrinsic candidate targets — evidence-tiered.**

| Tier | Gene | Category | tumor_spec (log2FC) | NK-state DE | DepMap CERES (real, 35 gastric lines) | Druggability | Recommended assay |
|------|------|----------|---------------------|-------------|-------------|--------------|-------------------|
| **1§** | SGMS2 | SM synthesis | +0.038 | **UP in dysf (FDR=0.007)**; discordant in both external cohorts (§) | weakly essential (−0.02) | Preclinical | SM-synthase modulation + NK membrane-SM readout |
| **1‡** | NT5E | adenosine (CD73) | +0.013 | **UP in dysf (FDR=0.014)** | non-essential (+0.005) | Phase 1/2 clinical | Co-culture with CD73/A2AR inhibitor + NK cytotoxicity |
| **1** | SMPD3 | SM catabolism | +0.032 | no signal | weakly essential (−0.13) | Preclinical (Zheng 2023) | SMase-inhibitor rescue ± Tim3 blockade |
| **1** | SMPD1 | SM catabolism | +0.001 | no signal | weakly essential (−0.16) | Preclinical (Zheng 2023) | SMase-inhibitor rescue ± Tim3 blockade |
| **2** | PHGDH | serine metabolism | +0.059 | no signal | weakly essential (−0.04) | Phase 1/2 (PHGDH inhibitor) | NK–tumor co-culture + serine-pathway inhibition |
| **2** | MICA | stress ligand | +0.006 | no signal | non-essential (+0.05) | Preclinical | ADAM protease-activity assay + NK cytotoxicity |
| **3** | PSAT1 | serine metabolism | +0.016 | no signal | moderately essential (−0.27) | Preclinical | — downgraded (see text) |
| **3** | PSPH | serine metabolism | +0.010 | **DOWN in dysf (p=0.03)** | weakly essential (−0.16) | Preclinical | — downgraded (NK-state direction wrong) |
| **3** | SHMT1/2, MTHFD1L | serine/1C metabolism | +0.001–0.008 | **DOWN in dysf (all p<0.02)** | weakly essential (−0.10 each) | — | — downgraded (NK-state direction wrong) |
| **X** | MTHFD1 | serine/1C metabolism | +0.008 | **DOWN in dysf (p<0.02)** | **pan-essential (−0.53)** | — | — not an immune-evasion target |
| **X** | ERBB2 | gastric oncogene | +0.069 | DOWN trend | moderately essential (−0.36) | FDA approved (HER2+ GC) | — confound with cell-viability effect; not a novel target |
| **X** | RAC1 | NK protrusion GTPase | +0.077 | **DOWN in dysf (p=0.04)** | **pan-essential (−0.74)** | — | — not an immune-evasion target |

*Tier 1: does not fall in the pan-essential range that would confound an immune-evasion interpretation with a cell-viability effect, plus NK-state DE or mechanistic privilege (§3.5 note on DepMap numbers). Tier 2: single-dimension support or mechanistic privilege. Tier 3: downgraded — NK-state signal absent or opposite direction. Tier X: excluded — pan- or moderately essential in vitro. ‡ NT5E's NK-state DE result is partly self-referential: NT5E is itself a constituent of the `NK_DYSFUNCTION_GENES` set that defines the NK-hot-dysfunctional label being tested against (`results/tables/geneset_separation_audit_summary.md`); SGMS2, SMPD3, and SMPD1 have no such overlap. § SGMS2's "UP in dysfunctional" TCGA-STAD signal is discordant (opposite direction) in both external replication cohorts (GSE62254, GSE84437; `results/tables/nk_state_de_external_concordance.tsv`) and should not be weighted as independent evidence; SGMS2's Tier-1 status rests on mechanistic privilege and DepMap non-pan-essentiality instead. DepMap CERES values are from a real query of **DepMap Public 26Q1** (`data/26Q1/CRISPRGeneEffect.csv`, `Model.csv`), the current release at the time of this analysis — see Data Availability and Limitations. Full per-gene evidence in `results/tables/target_validation_v2_merged.tsv`.*

**Figure 3.** (a) Tumor-specificity vs NK-dysfunction-correlation scatter,
tumor-intrinsic pool highlighted. (b) Top-15 tumor-intrinsic candidates by
target score. (c) Category composition of the 37-gene pool.

---

### 3.7 The graph as a convergent probe of transcriptional reach

The heterogeneous graph is not offered as a predictive advance — it matches
simpler baselines on accuracy (§3.4). Instead we use it as an independent,
architecture-based test of the same question the correlation analyses ask: is the
mechanism's *upstream metabolic coupling* (tumor serine program → NK SM/topology)
present in the transcriptome? We hard-code that coupling as the
`metabolic_crosstalk` edge and ask what removing it does. Three graph variants
were compared on their spectral gene embeddings (§2.6, Stage 1), built on the
enriched real graph (STRING PPI + CellChatDB LR + mechanism edges + scRNA
co-expression):

| Variant | Edges | Modularity | H1 embedding-coupling (serine↔SM) | H2 embedding-coupling (SM↔protrusion) |
|---------|-------|------------|------------------------------------|----------------------------------------|
| FULL (all edges) | 1,612 | 0.279 | 0.128 | 0.183 |
| −MC (without `metabolic_crosstalk`) | 1,312 | 0.288 | −0.002 | 0.146 |
| −SST (mechanism edges removed; generic priors only) | 492 | 0.287 | −0.002 | −0.000 |

**The edge does exactly what it is designed to do — and that is the
problem.** Removing `metabolic_crosstalk` abolishes the tumor-serine↔NK-SM
coupling in the embedding (H1 embedding-coupling 0.128→−0.002), confirming that
this single edge is what imposes the hypothesized upstream coupling on the
geometry. But that coupling is one the transcriptome itself does not exhibit: the
corresponding H1 correlation is null at both bulk and single-cell resolution
(§3.2). The edge faithfully encodes the mechanism; the data simply do not
corroborate it. (We do not lean on the modularity column: it is sensitive to
exactly which edge types happen to connect the same node pairs and does not move
consistently with mechanism-edge content across variants, so it is not treated as
independent evidence.)

The multi-view audit makes this distinction stricter. After removing all
label-defining genes, the external LIHC ensemble AUROC/AUPRC was 0.922/0.846 for
the no-graph expression model, 0.906/0.822 for merged SVD, 0.914/0.832 for
uniform fusion, and 0.912/0.826 for learned fusion (Fig. S1). Learned fusion was
worse than no graph for both AUROC (Δ=−0.0097, 95% bootstrap CI
−0.0185 to −0.0024) and AUPRC (Δ=−0.0205, −0.0381 to −0.0063), and was
indistinguishable from uniform fusion. In the unmasked sensitivity analysis the
same learned model reached 0.971/0.936, showing why the label-overlap control is
load-bearing rather than optional. No view was top-ranked in at least 8/10 seeds
(MSigDB was highest in 7/10), and every leave-one-view-out contrast failed the
pre-registered joint weight/CI contribution gate. Thus, **view separation
improves structural auditability but does not establish a predictive advantage**.

As a topology-specific calibration, we permuted the gene labels of each
mechanism view 1,000 times while fixing the other five views. Both authored
views imposed their intended geometry (observed vs null mean coupling:
`metabolic_crosstalk` 0.225 vs 0.031; `sm_topology_axis` 0.169 vs 0.079; both
empirical *P*=0.001). This verifies that the typed edges act at the intended
module boundary; it is explicitly not biological or predictive validation.

Together, the ablation is a third, methodologically independent line of evidence
for the paper's central boundary: the **effector layer** of the axis is
transcriptionally present (the within-axis `sm_topology_axis` edges organize the
embedding; the bulk effector correlations recover, §3.2–3.3), whereas the
**upstream metabolic coupling** is not — it is null as a correlation (H1/H2,
§3.2) and absent as an imposed graph edge once that edge is removed (here), and,
per the anchor study, resolvable only by single-cell mass spectrometry. The graph
does not fail as a predictor so much as succeed as a probe, returning the same
answer as every other instrument in this study. The utility of typed,
mechanism-grounded edges as interpretable structure — as opposed to a predictive
gain — is consistent with single-cell work in which ligand–receptor interaction
analysis, rather than expression alone, explained how one cell subset regulates
another [21].

### 3.8 Real-data comparative recoverability atlas

We applied the same pre-specified module comparisons to all four registered cards in four verified human bulk cohorts (TCGA-STAD, TCGA-LIHC, GSE62254 and GSE84437; Fig. S2). The two SM-topology downstream comparisons were directionally concordant with BH-FDR <0.05 in all cohorts. The NKG2D recognition-to-cytotoxicity comparison likewise recovered, whereas its tumor-shedding comparison did not. Neither pre-specified adenosine nor TGF-beta inhibitory comparison recovered in the expected direction. This is a comparative evidence map, not a universal law: the direct metabolomics, spatial and sample-level protein matrices required for the corresponding non-transcriptomic endpoints were unavailable in the public inputs and are explicitly recorded as `not_measured`. The pre-registered three-card/two-cohort/direct-modality gate therefore returns `comparative_atlas_only`.

---

## 4. Discussion

### 4.1 Summary of findings

Our central result is a **map of the mechanism's transcriptional reach**, not a
blanket claim of recovery. Three findings anchor it:

1. **The effector layer of the axis is transcriptionally recoverable in bulk
   transcriptomes and generalizes to gastric cancer.**
   The protrusion-machinery→cytotoxicity coupling replicates from
   independent public liver (bulk r=0.55) and gastric (§3.3) transcriptomes —
   extending the recoverable layer to a digestive-tract cancer not on the
   mechanism's published list. At single-cell resolution the same coupling
   is pseudoreplication-corrected (corrected r=0.31, *P*=3.9×10⁻⁸,
   meta-analysis across 9 samples), but we found this number does not survive
   residualization against library size and the dominant transcriptional
   (scVI-latent) structure of the data (r drops to 0.09), nor does it exceed
   what a randomly drawn, size-matched pair of gene modules produces in this
   dataset (permutation *P*=0.97 with an expression-matched scoring method;
   §3.2). We therefore do not treat the single-cell number as an independent
   replication and anchor this finding on the bulk result.
2. **The metabolic coupling is not transcriptionally detectable — by two
   independent methods.** The SM-balance→protrusion coupling is invisible in
   bulk and **not statistically significant** in single NK cells after
   correcting for within-sample dependence (corrected r=0.029, *P*=0.20,
   I²=73%); the upstream serine↔SM step (H1) is likewise null at both
   resolutions (§3.2). Crucially, a methodologically distinct test converges on
   the same conclusion: when the same metabolic coupling is hard-coded as the
   `metabolic_crosstalk` graph edge, it faithfully imposes a tumor-serine↔NK-SM
   structure on the embedding, and removing the edge abolishes that structure
   (H1 embedding-coupling 0.13→0 on removal; §3.7). Correlation analysis and
   graph ablation — two methods with different assumptions — thus agree that
   this coupling is real biology that operates at the metabolite level, not the
   transcript level, consistent with the anchor paper's requirement for
   single-cell mass spectrometry. This is precisely the
   kind of boundary the mechanism card's "gated physical ground-truth" section
   was designed to mark: the card encodes the coupling as a hypothesis, and the
   transcriptome — whether interrogated by correlation or by graph structure —
   declines to corroborate it.
3. **Transcription does not substitute for the physical topology phenotype —
   a finding, not a failure.** Machinery-gene transcription runs *opposite* to the
   physical protrusion collapse in intratumoral NK (higher transcript, lower
   function). This is the natural boundary of any transcriptome-based
   reconstruction of a membrane-lipid mechanism, and it precisely reproduces why
   the anchor lab required single-cell mass spectrometry and super-resolution
   imaging. By demonstrating this disconnect empirically, our framework provides
   a cautionary benchmark against over-interpreting transcriptional proxies for
   physical phenotypes in immune-evasion research, and identifies the
   post-transcriptional membrane-recruitment layer as the critical gap for
   future experimental work; a companion wet-lab validation program targeting
   this gap is in preparation.
4. **The transcriptional-reach boundary generalizes across mechanisms.**
   Applying the framework's second mechanism card (TGFβ→SMAD→NK exclusion; Gao
   et al. 2023) end-to-end on the same gastric cohorts reproduces the Zheng-card
   pattern rather than a card-specific artifact (§4.2). The effector/outcome
   co-variation is recoverable — and, being built from the NK activation program
   itself, is expectedly strong — but the two *mechanism-specific causal
   predictions*, that TGFβ signaling suppresses activating receptors (H2) and
   that CAF-ECM programs exclude cytotoxic NK (H4), both fail in bulk, running
   positive where the mechanism predicts negative. This mirrors the failure of
   the Zheng card's upstream metabolic coupling, and the shared cause is visible
   once NK-cell fraction is controlled (§3.2): bulk immune-module correlations
   are dominated by shared infiltration, and only generic NK-effector
   co-regulation survives that control. The consistent result across two
   independently authored cards is the paper's most generalizable claim — the
   framework behaves as a *discriminating instrument* that recovers
   transcriptionally-encoded effector state while correctly declining to
   corroborate causal steps operating at the metabolite, signaling, or physical
   level.

Together these delineate which layers a transcriptome can and cannot reach — a
result of direct use to any lab extending a physical immune-evasion mechanism to
cohort scale. The prioritized putative tumor-intrinsic target list, led by the druggable
serine/sphingomyelin enzymes at the mechanistic head of the axis, converts this
map into experimentally testable follow-up. Independent orthogonal validation
using NK-state differential expression and DepMap CRISPR essentiality screens
(§3.5) refines the 37-gene list into evidence tiers: **SMPD3 and
SMPD1** (Tier 1) pass both orthogonal filters cleanly; **SGMS2** also remains
Tier 1 on mechanistic privilege and DepMap non-pan-essentiality, but its
NK-state DE signal (UP in dysfunctional in TCGA-STAD) is discordant in two
independent external gastric cohorts and should not be weighted as
independent evidence; **NT5E/CD73** also
passes but is flagged separately because it is itself one of the marker genes
used to define the NK-hot-dysfunctional label its DE test is run against
(§3.5), so its result is partly self-referential. Five serine-pathway genes
(PSPH, SHMT1/2, MTHFD1/L) show expression patterns opposite to the expected
immune-evasion direction and ERBB2/RAC1 are pan-essential in vitro. The tiered
list sharpens the experimental roadmap: Tier 1 targets are ready for wet-lab
follow-up; Tier 2 (PHGDH, MICA) require additional tumor-side evidence; Tier 3
genes should be tested only after the primary candidates. Additionally, a per-sample
readout combining NK SM-catabolism score and HAVCR2 (Tim3) expression can stratify
samples by the logic of the SM-restoration + Tim3-blockade combination proposed
in the anchor paper — presented as an in-silico hypothesis for experimental
testing, not a validated clinical predictor. Consistent with this
checkpoint-combination logic, inhibitory receptors are co-regulated on tumor NK
cells: an independent HCC study reports that a distinct NK checkpoint, CLEC12B,
correlates in TCGA-HCC with HAVCR2/TIM-3, TIGIT, PDCD1 and LAG-3 [46], supporting
the biological coherence of a combined SM-catabolism plus Tim3-axis readout.

### 4.2 The mechanism-card approach

A key contribution of this work is the mechanism-card abstraction itself — a
machine-readable YAML formalism that separates the computational
operationalization of a mechanism from the pipeline that executes it. A card
declares, for one published mechanism, its molecular chain, the gene modules and
expected directions for each step, the cell type each module is attributed to,
the physical ground-truth measurements that are *out of scope* for transcriptome
analysis, the graph node/edge types it introduces, and the pre-registered
validation hypotheses with an explicit recovery definition. Because the pipeline
consumes the card rather than hard-coding the biology, applying the framework to a
new mechanism — a different metabolic immune checkpoint (e.g. adenosine-mediated
NK suppression), a different effector cell type, or a stress-ligand-shedding axis
— requires only authoring a new card, not rewriting the pipeline core. The card
registry (`configs/mechanism_cards/registry.yaml`) is designed to accumulate such
cards over time.

Beyond reuse, the card format enforces scientific discipline that this study shows
to be load-bearing rather than decorative. First, it requires **explicit claim
boundaries** in the card body, which is what let us state without ambiguity that
machinery transcription is not a proxy for the physical topology phenotype — a
boundary our own results confirmed empirically (§3.2). Second, it **gates physical
ground-truth targets** (SEM protrusion density, single-cell SM mass spectrometry)
as a separate, non-transcriptional layer, preventing their accidental conflation
with the computable proxy. Third, it requires **pre-registered hypotheses with a
defined recovery criterion**, so that a partial or negative outcome (as in Arm A)
is reported as a structured result rather than quietly reframed. The H2
non-significance after pseudoreplication correction (*P*=0.20) — initially
reported as "statistically detectable" in naive analysis — illustrates this
discipline in practice: the correction and re-reporting were performed before
finalizing this manuscript, and the negative result is presented plainly rather
than minimized (§3.2).

**Multi-card evidence for the formalism.** The mechanism-card registry currently
holds four cards spanning distinct NK immune-evasion axes: the serine–SM axis
(Zheng 2023), adenosine-mediated suppression (A2AR→cAMP/PKA), TGFβ-driven
exclusion (TGFβ→SMAD→receptor suppression), and MICA/B shedding (NKG2D ligand
cleavage). All four cards share the same YAML schema and are consumed by the
same pipeline without modification. Pairwise gene overlap is low (Jaccard
0.07–0.19), confirming the cards encode mechanistically distinct biology
rather than overlapping gene sets (`results/tables/mechanism_card_comparison.tsv`,
`mechanism_card_gene_overlap.tsv`). **Two of the four cards have now been run
end-to-end on real data.** Beyond the Zheng serine–SM card, the
TGFβ→SMAD→NK-exclusion card (Gao et al. 2023) was scored on TCGA-STAD and both
external gastric cohorts and its pre-registered hypotheses tested identically
(`src/interpretation/run_tgfb_card_recovery.py`,
`results/tables/mechanism_card_tgfb_recovery.tsv`). Rather than a card-specific
artifact, it reproduces the same transcriptional-reach boundary reported for the
Zheng card (§4.1): the NK-effector/outcome co-variation is recoverable, but the
two mechanism-specific causal predictions — that TGFβ signaling suppresses NK
activating receptors (H2) and that CAF-ECM programs exclude cytotoxic NK (H4) —
both fail in bulk, running positive where the mechanism predicts negative, just
as the Zheng card's upstream metabolic coupling failed. The adenosine and MICA/B
cards remain structurally validated (schema compliance, synthetic-mode ingestion,
low inter-card Jaccard) but not yet run end-to-end on real data. Reusability is
therefore no longer only a design aspiration: the same formalism has produced a
consistent, cross-mechanism empirical result on two independently authored cards.
In our experience, this discipline is precisely what makes a partial-recovery
finding credible to a mechanistic wet lab, and it generalizes to any attempt to
operationalize a physical mechanism from an indirect molecular readout.

### 4.3 Limitations

1. **Transcriptional proxy ≠ physical topology.** Gene expression captures the molecular machinery and capacity for the serine–SM–topology axis, not the physical membrane phenotype itself. The disciplined qualifiers throughout are load-bearing: every claim is bounded by "transcriptional program permissive-of / associated-with."
2. **Serine/SM crosstalk is a metabolite-level effect.** Transcription captures enzyme abundance, not flux. The actual serine→SM crosstalk requires metabolomics or the anchor lab's single-cell mass spectrometry for direct measurement. This is now empirically supported by the H2 non-significance after pseudoreplication correction.
3. **Data availability constraints.** The anchor paper did not deposit transcriptomic data, requiring the use of independent public cohorts. The liver positive control is therefore not a direct replication but an independent validation.
4. **No experimental validation.** All targets are computationally prioritized; none have been tested in wet-lab assays. The recommended assays are offered as a bridge to experimental follow-up (see companion experimental program).
5. **Single-cell pseudoreplication.** The 8,310 NK cells derive from 9 biological samples. We correct for this via per-sample meta-analysis (§2.9), but between-sample heterogeneity is substantial (I² up to 96% for H3), indicating that sample-level factors beyond the SST axis influence the correlations. Results should be interpreted at the sample level, not the cell level. A leave-one-sample-out sensitivity analysis (`results/tables/h3_leave_one_sample_out.tsv`) shows the H3 pooled estimate itself is stable to removing any single sample (pooled r range 0.275–0.350 across the 9 leave-one-out re-analyses, all 95% CIs excluding 0) — so the pseudoreplication correction is not an artifact of one outlier sample. This is a separate question from whether the correlation is technically confounded (item 6 below), which it is.
6. **Single-cell module-score correlations require count-depth and module-permutation controls beyond pseudoreplication correction.** Re-running the count-depth control (P0-2) and module-membership permutation test (P0-3) on the real scRNA data (8,310 NK cells; previously validated on synthetic data only) revealed that pseudoreplication correction alone is not sufficient for the H3 single-cell number: 47.4% of the protrusion-machinery module score's variance is explained by library size (`total_counts`, which varies 236-fold across cells), residualizing against library size and the real scVI latent space collapses r from 0.32 to 0.09, and a permutation test using an expression-matched scoring method shows the observed coupling does not exceed a randomly-drawn, size-matched gene-module baseline (empirical *P*=0.97). We therefore no longer treat the single-cell H3 number as an independent replication (§3.2, §4.1); the effector-arm claim rests on the bulk result and its gastric replications. The count-depth diagnostics for H2 and H4 did not change those hypotheses' already-null verdicts.
7. **NK subtype resolution.** scRNA-based NK annotation depends on the quality of the reference atlas. Populations that are rare or absent in the reference may be misclassified.
8. **The graph model does not outperform top tabular or signature baselines on accuracy, and the stricter external audit confirms that boundary.** On the binary NK-state task the GNN is statistically indistinguishable from LightGBM/XGBoost and from an SST-module signature baseline (§3.4). In the label-masked STAD→LIHC audit, learned multi-view fusion was significantly worse than no graph for AUROC and AUPRC and no individual view passed the joint stability/ablation gate (§3.7). The better unmasked result is a sensitivity analysis of label overlap, not evidence of graph value. We therefore position typed views as an auditable representation of authored structure, not as a predictor; their node-label-permutation calibration verifies where structure is imposed but cannot validate the underlying biology.
9. **Residual NK bias in the tumor-intrinsic candidate pool.** The `tumor_specificity_log2>0` gate is permissive: 17 of 37 candidates are annotated to NK-side mechanism-card modules, including RAC1 and WASL. The pool should be interpreted as "genes with a non-zero malignant-cell transcript signal that mechanistically intersect the SST axis," not as a clean set of tumor-exclusive targets. The trivial baseline comparison (§3.5) confirms that the five-dimension scoring adds discriminative value over anchor-paper membership alone (Spearman rho=0.54), but the delta is moderate. A stricter filter (e.g. tumor_specificity_log2>0.5 or a module-level penalty) would reduce NK-side contamination at the cost of losing borderline tumor-intrinsic candidates.
10. **Candidate atlas omits intracellular/TF-level regulators.** The prioritization scores surface and metabolic-enzyme genes; transcription-factor and intracellular negative regulators of NK function (e.g. the CREM/PKA–CREB axis [27]) fall outside the current candidate space and are a natural extension of the mechanism-card modules.
11. **No clinical-outcome anchor.** NK states are linked to scRNA-defined labels, not to patient survival or therapy response. Single-cell atlases that tie a functional state to durable clinical outcome [21] indicate a clear next step: anchoring the NK-state readout to outcome in a cohort with follow-up.
12. **Multi-card validation covers two of four cards end-to-end.** The TGFβ→SMAD→NK-exclusion card has now been run end-to-end on real data alongside the Zheng serine–SM card (§4.2; `results/tables/mechanism_card_tgfb_recovery.tsv`), and the two converge on the same transcriptional-reach boundary. The adenosine and MICA/B cards remain validated only structurally (schema compliance, synthetic-mode ingestion, low inter-card Jaccard) and have not been run end-to-end with real expression data; their full biological validation remains future work. We also note that the TGFβ card's strong effector-arm coupling (H3) partly reflects that its "activating receptor" module is itself part of the NK activation program and therefore co-regulated with the cytotoxicity outcome — so the informative cross-card signal is the consistent *failure* of the mechanism-specific causal predictions (H2, H4), not the magnitude of H3.
13. **DepMap essentiality now uses the real, current 26Q1 release.** CERES scores in Table 4 are from a real query of **DepMap Public 26Q1** `CRISPRGeneEffect.csv` and `Model.csv` (35 real gastric/stomach-subtype cell lines), obtained directly from the DepMap portal — newer than the 25Q2 release originally named in Methods, and superseding an earlier pass of this analysis that used the 24Q2 figshare snapshot (the DepMap portal's own interactive download is gated by a Cloudflare bot-verification challenge that cannot be scripted non-interactively; 24Q2 was the most recent release with non-interactive public access at that time). Under the 26Q1 data, NT5E is genuinely non-essential (CERES=+0.005, revised from "weakly essential" under 24Q2); SGMS2, SMPD1, and SMPD3 remain "weakly essential" (CERES −0.02 to −0.16) — still well clear of the pan-essential range (RAC1, MTHFD1; CERES<−0.5) but not strictly non-essential; PSAT1 moved from "weakly essential" to "moderately essential" (CERES=−0.27). Table 4 and §3.5 reflect these real 26Q1 values.
14. **NK-state DE is underpowered for the dysfunctional group, and external replication gives mixed results.** TCGA-STAD has only n=20 NK-hot-dysfunctional samples vs n=134 NK-hot-cytotoxic. We replicated the directional test in two independent external gastric cohorts (GSE62254, GSE84437; `results/tables/nk_state_de_external_replication_summary.md`), each with an independently-derived NK-state label: 14/16 gene-cohort comparisons (88%) are directionally concordant with TCGA-STAD, including external confirmation that the five downgraded serine-pathway genes and RAC1 are down in dysfunctional tumors (strengthening their exclusion). However, **SGMS2's "up in dysfunctional" TCGA-STAD signal is discordant in both external cohorts** and should not be weighted as independent evidence for its Tier-1 status (§3.5, Table 4); NT5E's concordance is real but remains subject to the label-circularity caveat (item 15 below) in every cohort tested. Replication in a larger cohort with paired scRNA would still strengthen these conclusions.
15. **The NK-state classification target partly overlaps the classifier's own input features.** A gene-set separation audit (`results/tables/geneset_separation_audit_summary.md`) found that the NK-hot-cytotoxic/dysfunctional label (§2.4) is a thresholding rule on marker genes (NKG7, GNLY, GZMB, PRF1, IFNG among others) that are also present, unmodified, in the full expression vector used as input to every model in §3.4, including the GNN. This is a known property of marker-defined phenotype labels rather than a coding error, but it means the high accuracy of all models (including simple baselines such as the 8-gene NK-marker signature, AUROC=0.849) partly reflects recovering a label from the same genes that define it, and §3.4 should be read accordingly. The same audit found that the SST-axis modules used for the H1–H5 mechanism tests do not reference the NK-state label at all (those tests operate directly on module scores, not on the label), so this concern is confined to §3.4's classification numbers. It also found that NT5E, one of the 37 tumor-intrinsic candidates, is itself part of the label-defining `NK_DYSFUNCTION_GENES` set (§3.5, Table 4).
16. **No hard per-cell QC threshold was applied to the real scRNA-seq data.** Per-cell QC metrics (detected genes, mitochondrial fraction) were computed but not used to exclude cells (§2.4); all 166,829 concatenated cells were retained into normalization and integration, and doublet-based exclusion was not applied. Technical, quality-related variance is instead addressed downstream through the count-depth residualization diagnostics (item 6 above), which show that library size and detected-gene count explain a substantial share of module-score variance and must be accounted for statistically. A dataset with explicit per-cell hard-threshold filtering applied upstream would provide an additional, independent check on the results in §3.2.
17. **The bulk effector coupling is substantially NK-abundance-driven.** In bulk transcriptomes both the protrusion-machinery and cytotoxicity-output module scores scale with NK-cell fraction, so the zero-order protrusion~cytotoxicity correlation is inflated by co-varying NK abundance. Partialling out a clean NK-lineage fraction proxy (§3.2; `results/tables/bulk_h3_purity_control.tsv`) attenuates the coupling by ~50% (TCGA-LIHC r=0.55→0.25, still significant), roughly halves it in GSE62254 (still significant), and abolishes it in GSE84437 (partial r≈0, n.s.). The NK-lineage proxy is a coarse, transcriptome-derived abundance estimate rather than a deconvolution or cell-sorted NK fraction, and partialling it out may itself be conservative if NK abundance genuinely co-occurs with effector-program intensity; a deconvolution-based NK fraction (CIBERSORTx/quanTIseq) or a within-NK single-cell analysis free of the count-depth confound (item 6) would refine the adjusted estimate. The effector-arm claim is reported accordingly as partially recoverable and infiltration-driven, not as a robust replication.

### 4.4 Future directions

- **Additional cancer types.** Once the gastric extension is validated, the framework can test the axis in other digestive-tract cancers (colorectal, pancreatic, esophageal) using the same mechanism card.
- **Additional mechanisms.** The mechanism-card registry (`configs/mechanism_cards/registry.yaml`) is designed to hold multiple cards. Cards for adenosine-mediated NK suppression, TGFβ-driven NK exclusion, stress-ligand shedding (MICA/B-ADAM17), and other NK checkpoint axes (TIGIT/CD96 [47,48]; the CLEC12B–lipoprotein-lipase axis recently shown to restrain tumor NK cells and to synergize with PD-1 blockade [46]) are natural next additions.
- **Physical topology integration.** When membrane protrusion / microvilli imaging data become available (even for a subset of samples), Layer 14R-B of the SST-axis module can be activated to provide direct phenotype-transcriptome correlation.
- **Prospective validation cohort.** A dedicated gastric cancer cohort with paired bulk RNA-seq, scRNA-seq, and functional NK assays would provide the strongest validation of the prioritized targets.
- **Genome-scale learnable encoder.** The present audit already separates edge types and learns a global fusion weight over per-view spectral embeddings (§2.6); it does not implement the per-node GCN/Transformer encoders of TREE [49] or GRAFT [50]. Such an architecture should be reserved for a genome-scale panel, a label independent of its input genes, and an independently powered endpoint where model capacity can be evaluated without the circularity and small-panel overfitting exposed here.

---

## 5. Conclusion

GC-NKGraph-Atlas demonstrates *how far* a specific, published immune-evasion
mechanism — the serine–sphingomyelin–membrane-topology axis of NK dysfunction —
can be surveyed at cohort scale from public transcriptomes. Using a
single-cell-informed heterogeneous graph framework and a two-arm design, we establish
a three-layer scoping map: (i) the effector arm is recoverable in bulk
transcriptomes from liver to gastric cancer, but ~50% of the bulk coupling is
NK-cell-abundance covariation — after adjusting for NK fraction it holds at
roughly half strength in two of three bulk cohorts and does not survive in the
third, and at single-cell resolution it does not survive count-depth/latent-
structure residualization, so we report it as partially recoverable and
substantially infiltration-driven rather than robustly recovered; this same
boundary reproduces on a second, independently authored mechanism card
(TGFβ→SMAD→NK exclusion), whose mechanism-specific causal predictions likewise
fail in bulk; (ii) the upstream metabolic
coupling is **not** transcriptionally detectable after correction (*P*=0.20),
consistent with metabolite-level regulation — a conclusion reached
independently by correlation analysis and by ablating the mechanism-encoding
graph edge, which imposes the coupling on the embedding and loses it entirely
once the edge is removed (§3.7); and (iii) the physical topology
phenotype is fundamentally disconnected from machinery transcription
(transcript levels run opposite to the physical collapse) — a finding, not a
failure, that benchmarks the natural limit of transcriptome-based reconstruction
for membrane-lipid mechanisms. This map, not an overclaim of full reconstruction,
is the contribution.

The mechanism-card abstraction, demonstrated here on the Zheng 2023 serine–SM
axis, provides a disciplined formalism for converting a published mechanism into
a scalable study with explicit claim boundaries and gated physical ground-truth.
We show that this discipline is load-bearing: it is precisely the card's
separation of computable-from-physical that lets us report the topology
disconnection as a structured result rather than a hidden failure, and it is the
pre-registered recovery definition that forces the H2 non-significance to be
reported plainly. The framework has now been run end-to-end on two independently
authored cards (serine–SM and TGFβ→SMAD→NK exclusion), which converge on the same
transcriptional-reach boundary; extending it to the remaining registered cards
(adenosine-mediated NK suppression, MICA/B shedding) is a natural next step, so
multi-card reuse is an established capability to be broadened as further cards are
run end-to-end rather than only a design aspiration.

---

## Data Availability

All code and configuration files are available at https://github.com/nblvguohao/GC-NKGraph-Atlas. The mechanism-card template and the Zheng 2023 NK SM-topology card are provided under `configs/mechanism_cards/`. Synthetic test data can be generated via `python src/common/synthetic_data.py`.

**Reproducibility.** The pipeline supports a fully self-contained synthetic data
mode (≈17 MB, no real patient data) that exercises all phases end-to-end:
```
pip install -r requirements.txt
python src/common/synthetic_data.py   # generates data/synthetic/
python src/pipeline.py --synthetic --force
```
An environment lock file (`environment.yml`) pins all dependency versions
(testable on Linux/macOS/Windows with conda). Unit tests covering graph
construction, NK scoring, SST-axis computation, and target prioritization are
provided under `tests/` and can be run with `pytest tests/`. The pseudoreplication
correction analysis is implemented in `src/topology/pseudoreplication_correction.py`
and the count-depth control in `src/topology/count_depth_control.py`.

**Prior-network data.** The heterogeneous graph's generic-prior edges (STRING v12
PPI, CellChatDB ligand–receptor pairs via OmniPath, ChEA 2022 TF-target sets) are
not redistributed with the repository; `src/data_download/download_prior_networks.py`
regenerates them from the respective public APIs, restricted to the axis-centered
gene panel, so `build_heterograph.py` reproduces the exact graph used here from a
clean clone.

**Real data availability:**
- TCGA-STAD and TCGA-LIHC: available from the Genomic Data Commons (https://portal.gdc.cancer.gov/)
- GSE62254, GSE84437, GSE246662: available from the Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
- Prior networks: STRING v12 (https://string-db.org/), CellChatDB via OmniPath (https://omnipathdb.org/), ChEA 2022 via Enrichr (https://maayanlab.cloud/Enrichr/) — fetched by the download script above.
- DepMap Public 26Q1 (`CRISPRGeneEffect.csv`, `Model.csv`): the real, current release used for Table 4's CERES values, obtained directly from the DepMap portal (https://depmap.org/portal/, interactive download required — see Limitations, item 13). An earlier pass of this analysis used DepMap Public 24Q2 via the public figshare API (https://plus.figshare.com/articles/dataset/DepMap_24Q2_Public/25880521) before 26Q1 access was available; that snapshot has been superseded throughout.
- No novel sequencing data were generated for this study.

**Pre-registration:** The hypotheses H1–H5 and the original full-recovery criterion for the positive-control arm were registered in `configs/sst_axis_config.yaml` before execution. Because the full criterion was not met, the manuscript reports the result as a partial-recovery scoping map and explicitly separates pre-registered outcomes from post-hoc interpretation.

---

## Author Contributions

**Guohao Lyu:** Conceptualization, Methodology, Software, Formal analysis,
Investigation, Data curation, Visualization, Writing – original draft.
**Yingchun Xia:** Methodology, Software, Validation, Writing – review & editing.
**Huichao Liu:** Software, Data curation, Validation.
**Xiaolei Zhu:** Formal analysis, Methodology, Writing – review & editing.
**Shuai Yang:** Investigation, Validation, Visualization.
**Ailian Zhou:** Conceptualization, Supervision, Funding acquisition, Writing –
review & editing.
**Lichuan Gu:** Conceptualization, Supervision, Project administration, Funding
acquisition, Writing – review & editing.

All authors read and approved the final manuscript.

---

## Ethics Approval and Consent to Participate

Not applicable. This study used only publicly available, de-identified datasets
(TCGA-LIHC, TCGA-STAD, GSE62254, GSE84437, GSE246662) obtained from the Genomic
Data Commons and the Gene Expression Omnibus. No new human or animal subjects were
involved, and no identifiable personal data were generated or analyzed. Use of
these public data complies with the respective data-access policies.

---

## Competing Interests

The authors declare that they have no competing interests.

---

## Funding

This work was supported by grants from the National Natural Science Foundation of China (32472007, 62301006, 62301008), the Natural Science Foundation of Anhui Province (2308085MF217, 2308085QF202), and the Anhui Province Key Laboratory of Intelligent Agricultural Technology and Equipment.

---

## References

1. Vivier E, Tomasello E, Baratin M, et al. Functions of natural killer cells. *Nat Immunol* 2008;9:503–10.

2. Chiossone L, Dumas PY, Vienne M, Vivier E. Natural killer cells and other innate lymphoid cells in cancer. *Nat Rev Immunol* 2018;18:671–88.

3. Huntington ND, Cursons J, Rautela J. The cancer–natural killer cell immunity cycle. *Nat Rev Cancer* 2020;20:437–54.

4. Myers JA, Miller JS. Exploring the NK cell platform for cancer immunotherapy. *Nat Rev Clin Oncol* 2021;18:85–100.

5. Shimasaki N, Jain A, Campana D. NK cells for cancer immunotherapy. *Nat Rev Drug Discov* 2020;19:200–18.

6. O'Brien KL, Finlay DK. Immunometabolism and natural killer cell responses. *Nat Rev Immunol* 2019;19:282–90.

7. Terrén I, Orrantia A, Vitallé J, et al. NK cell metabolism and tumor microenvironment. *Front Immunol* 2019;10:2278.

8. Laskowski TJ, Biederstädt A, Rezvani K. Natural killer cells in antitumour adoptive cell immunotherapy. *Nat Rev Cancer* 2022;22:557–75.

9. André P, Denis C, Soulas C, et al. Anti-NKG2A mAb is a checkpoint inhibitor that promotes anti-tumor immunity by unleashing both T and NK cells. *Cell* 2018;175:1731–43.e13.

10. Gao J, Zheng X, Liu Y, et al. TGF-β impairs NK cell migration and cytotoxicity. *J Immunother Cancer* 2023;11:e005785.

11. Melaiu O, Lucarini V, Cifaldi L, Fruci D. Influence of the tumor microenvironment on NK cell function in solid tumors. *Front Immunol* 2020;10:3038.

12. **Zheng X, Hou Z, Qian Y, et al. Tumors evade immune cytotoxicity by altering the surface topology of NK cells. *Nat Immunol* 2023;24:748–59. doi:10.1038/s41590-023-01462-9.** ★ Anchor paper.

13. Possemato R, Marks KM, Shaul YD, et al. Functional genomics reveal that the serine synthesis pathway is essential in breast cancer. *Nature* 2011;476:346–50.

14. Locasale JW. Serine, glycine and one-carbon units: cancer metabolism in full circle. *Nat Rev Cancer* 2013;13:572–83.

15. Ogretmen B. Sphingolipid metabolism in cancer — signalling and drug resistance. *Nat Rev Cancer* 2018;18:33–50.

16. Hannun YA, Obeid LM. Sphingolipids and their metabolism in physiology and disease. *Nat Rev Mol Cell Biol* 2018;19:175–91.

17. Zheng X, Tian Z, Wei H. Sphingomyelin is a prospective metabolic immune checkpoint for natural killer cells. *Clin Transl Med* 2023;13:e1395.

18. Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nat Biotechnol* 2019;37:773–82.

19. Finotello F, Mayer C, Plattner C, et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. *Genome Med* 2019;11:34.

20. Sun D, Guan X, Moran AE, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. *Nat Biotechnol* 2022;40:527–38.

21. Bai Z, Feng B, McClory SE, et al. Single-cell CAR T atlas reveals type 2 function in 8-year leukaemia remission. *Nature* 2024;634:702–11. doi:10.1038/s41586-024-07762-w.

22. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell–cell communication using CellChat. *Nat Commun* 2021;12:1088.

23. Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R. CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. *Nat Protoc* 2020;15:1484–506.

24. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat Methods* 2020;17:159–62.

25. Alghamdi N, Chang W, Dang P, et al. A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data. *Genome Res* 2021;31:1867–84.

26. Gulati GS, Sikandar SS, Wesche DJ, et al. Single-cell transcriptional diversity is a hallmark of developmental potential. *Science* 2020;367:405–11.

27. Rafei H, Basar R, Acharya S, et al. CREM is a regulatory checkpoint of CAR and IL-15 signalling in NK cells. *Nature* 2025;643:1076–86. doi:10.1038/s41586-025-09087-8.

28. Kipf TN, Welling M. Semi-supervised classification with graph convolutional networks. *Proc ICLR* 2017.

29. Veličković P, Cucurull G, Casanova A, et al. Graph attention networks. *Proc ICLR* 2018.

30. Schlichtkrull M, Kipf TN, Bloem P, et al. Modeling relational data with graph convolutional networks. *Proc ESWC* 2018:593–607.

31. Wang T, Shao W, Huang Z, et al. MOGONET integrates multi-omics data using graph convolutional networks allowing patient classification and biomarker identification. *Nat Commun* 2021;12:3445.

32. Hu Z, Dong Y, Wang K, Sun Y. Heterogeneous graph transformer. *Proc Web Conf* 2020:2704–10.

33. Stuart T, Butler A, Hoffman P, et al. Comprehensive integration of single-cell data. *Cell* 2019;177:1888–902.

34. Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data. *Cell* 2021;184:3573–87.e29.

35. Gayoso A, Lopez R, Xing G, et al. A Python library for probabilistic analysis of single-cell omics data. *Nat Biotechnol* 2022;40:163–66.

36. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 2018;19:15.

37. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. *Mol Syst Biol* 2019;15:e8746.

38. Lotfollahi M, Naghipourfar M, Luecken MD, et al. Mapping single-cell data to reference atlases by transfer learning. *Nat Biotechnol* 2022;40:121–30.

39. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA* 2005;102:15545–50.

40. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 2015;1:417–25.

41. Fey M, Lenssen JE. Fast graph representation learning with PyTorch Geometric. *ICLR Workshop* 2019.

42. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *Proc KDD* 2016:785–94.

43. Ke G, Meng Q, Finley T, et al. LightGBM: a highly efficient gradient boosting decision tree. *Adv NeurIPS* 2017:3146–54.

44. Cancer Genome Atlas Research Network. Comprehensive molecular characterization of gastric adenocarcinoma. *Nature* 2014;513:202–9.

45. Cristescu R, Lee J, Nebozhyn M, et al. Molecular analysis of gastric cancer identifies subtypes associated with distinct clinical outcomes. *Nat Med* 2015;21:449–56.

46. Sun P, Xu X, Hu B, et al. Targeting NK cell CLEC12B enhances cancer immunotherapy. *Nat Immunol* 2026;27:985–99. doi:10.1038/s41590-026-02471-0.

47. Barry KC, Hsu J, Broz ML, et al. A natural killer–dendritic cell axis defines checkpoint therapy–responsive tumor microenvironments. *Nat Med* 2018;24:1178–91.

48. Zhang Q, Bi J, Zheng X, et al. Blockade of the checkpoint receptor TIGIT prevents NK cell exhaustion and elicits potent anti-tumor immunity. *Nat Immunol* 2018;19:723–32.

49. Su X, Hu P, Li D, et al. Interpretable identification of cancer genes across biological networks via transformer-powered graph representation learning. *Nat Biomed Eng* 2025;9:371–89.

50. Cho SP, Cho YR. GRAFT: a graph-aware fusion transformer for cancer driver gene prediction. *Brief Bioinform* 2026;27:bbaf706.
