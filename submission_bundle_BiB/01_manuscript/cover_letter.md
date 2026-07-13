**To:** The Editors, *Briefings in Bioinformatics*
**Re:** Submission of a Problem Solving Protocol

Date: July 12, 2026

Dear Editors,

We are pleased to submit our manuscript, **"Transcriptional Reach of the
Serine–Sphingomyelin–Membrane-Topology Axis in NK-Cell Immune Evasion: A
Single-Cell Heterogeneous Graph Framework from Liver to Gastric Cancer,"**
for consideration as a Problem Solving Protocol in *Briefings in Bioinformatics*.

**The problem.** A landmark wet-lab study (Zheng et al., *Nature Immunology* 2023)
showed that tumors evade natural killer (NK)–cell cytotoxicity through a physical
mechanism: dysregulated tumor serine metabolism depletes sphingomyelin in NK
membranes, collapsing the membrane protrusions required for lytic immune
synapses. That mechanism was established with single-cell mass spectrometry and
super-resolution imaging on a small liver-cancer cohort — gold-standard evidence
that, by construction, cannot scale to hundreds of tumors or be screened across
cancer types. No transcriptomic data were deposited. This leaves an open,
methodologically interesting question: **how much of a physical, metabolite-level
immune-evasion mechanism can actually be reconstructed from the transcriptome —
and where does transcription stop being a valid proxy?**

**What we contribute.** We introduce GC-NKGraph-Atlas, a framework built around a
reusable *mechanism-card* abstraction that encodes a published wet-lab mechanism
as a machine-readable recipe, and a single-cell-informed, NK-aware heterogeneous
graph with a biology-grounded `metabolic_crosstalk` edge. Rather than claiming a
blanket "recovery," our two-arm design (liver positive control; gastric-cancer
extension) yields a **precise scoping result** that we believe is the paper's most
useful and honest contribution:

- The **effector layer** of the axis (protrusion-machinery → cytotoxicity)
  reproduces robustly from independent public bulk transcriptomes — liver
  (r = 0.55) and two external gastric cohorts (r = 0.42, 0.62) — and
  generalizes to gastric cancer. The same coupling is nominally present in
  8,310 single NK cells after pseudoreplication correction, but we found this
  single-cell number does not survive further technical-confound controls
  (count-depth/latent-structure residualization; a random-gene-module
  permutation baseline), so we anchor the effector-arm claim on the bulk
  evidence alone, not on the single-cell number.
- The **upstream metabolic coupling** is undetectable in bulk and, after the
  same pseudoreplication correction, is **not** statistically significant in
  single NK cells either (corrected r = 0.029, P = 0.20) — an empirical
  demonstration that this coupling operates at the metabolite level rather
  than the transcript level, consistent with the anchor mechanism.
- The **physical membrane-topology phenotype** is demonstrably *not* captured by
  machinery-gene transcription, which we show and quantify rather than assume —
  delineating the honest boundary of any transcriptome-based reconstruction.

On NK-state classification the graph-informed model reaches AUROC 0.95 / MCC 0.71
under 5-fold cross-validation — statistically on par with the strongest
gradient-boosting baselines (LightGBM/XGBoost; paired p>0.27) and significantly
above linear, kernel, and shallow-network baselines (p<0.05) — while additionally
yielding the mechanism-structured gene embedding used for the axis analyses. The
framework then prioritizes 37 tumor-intrinsic candidate targets — led by the
druggable serine/sphingomyelin enzymes at the head of the mechanism (PHGDH,
SGMS2, SMPD3/1) — each paired with a recommended wet-lab validation assay, and
kept strictly separate from the NK-side axis readout. The recovered effector
coupling further replicates in two fully independent gastric microarray cohorts
(GSE62254, r=0.42; GSE84437, r=0.62; both p≪10⁻¹³), so the reproducible layer of
the axis is confirmed across three independent bulk datasets (liver plus two
gastric cohorts); the single-cell gastric subset shares the same technical
confound diagnosed in the liver arm and is not treated as an independent
confirmation.

**Why *Briefings in Bioinformatics*.** The manuscript is a methods-and-benchmarking
contribution aimed squarely at computational-biology readers: a transferable
formalism for operationalizing published mechanisms, a heterogeneous-graph design
with an explicit biological justification for each edge type, a rigorous
multi-resolution validation protocol, and — unusually — an honest map of a
method's limits. We think its combination of a reusable abstraction, a
benchmarked model, and disciplined claim boundaries fits the journal's readership
and its Problem Solving Protocol format.

**Declarations.** This work is original, has not been published previously, and is
not under consideration elsewhere. All authors have read and approved the
manuscript and agree to its submission. No novel sequencing data were generated;
all analyses use public cohorts (TCGA-LIHC, TCGA-STAD, GSE62254, GSE84437,
GSE246662), and all code and configuration are available for review at
https://github.com/nblvguohao/GC-NKGraph-Atlas. The authors declare no competing
interests.

**Suggested reviewers** (no conflicts of interest with the authors):

1. **Prof. Satu Mustjoki**, Translational Immunology Research Program, University of Helsinki, Helsinki, Finland. E-mail: satu.mustjoki@helsinki.fi. Expertise: natural killer cell functional genomics, single-cell transcriptomics of hematological and solid-tumor immune evasion, CRISPR-based dissection of NK resistance mechanisms.

2. **Prof. Fabian J. Theis**, Institute of Computational Biology, Helmholtz Zentrum München / Technical University of Munich, Munich, Germany. E-mail: fabian.theis@helmholtz-muenchen.de. Expertise: single-cell computational methods (scVI, scANVI, scArches), deep generative models for transcriptomics, transfer learning in single-cell atlases.

3. **Prof. Xia Li**, College of Bioinformatics Science and Technology, Harbin Medical University, Harbin, China. E-mail: lixia@hrbmu.edu.cn. Expertise: computational immune bioinformatics, transcription factor–mediated regulatory network analysis, single-cell data mining for tumor immunity. (Published in *Briefings in Bioinformatics*.)

We would be grateful for the opportunity to have this work considered, and we are
happy to provide any additional information. Correspondence should be addressed to
Prof. Lichuan Gu (glc@ahau.edu.cn) and Prof. Ailian Zhou (zhouailian@caas.cn).

Sincerely,

Guohao Lyu, on behalf of all authors
School of Artificial Intelligence, Anhui Agricultural University, Hefei 230036, China
lvguohao@stu.ahau.edu.cn
