# Supplementary Results — Extended Narrative

Manuscript: *Transcriptional Reach of the Serine–Sphingomyelin–Membrane-Topology
Axis in NK-Cell Immune Evasion* (GC-NKGraph-Atlas).

**Why this document exists.** The main text is condensed to meet *Briefings in
Bioinformatics*'s word-count limit for a Problem Solving Protocol (2,000–5,000
words). Every number, verdict, and conclusion in the main text is unchanged;
this document restores the full step-by-step statistical reasoning behind the
headline results, for readers and reviewers who want the complete derivation
rather than the verdict alone. All underlying values are also in the
Supplementary Tables (`SUPPLEMENTARY_INDEX.md`), cross-referenced below.

---

## S3.2 — Arm A confound-control battery (full narrative)

**H3 (effector arm) single-cell diagnostics.** The bulk TCGA-LIHC
protrusion-machinery↔cytotoxicity-output coupling (r=0.551, P=5.5×10⁻³⁵) also
appears, at face value, across 8,310 single NK cells. After correcting for
within-sample dependence (per-sample r → Fisher z → DerSimonian-Laird
random-effects meta-analysis across 9 samples), the corrected single-cell
meta-analytic r=0.313 (P=3.9×10⁻⁸, 95% CI [0.13, 0.48]), with I²=96% reflecting
substantial between-sample heterogeneity (per-sample r range [0.009, 0.560]).

Pseudoreplication correction addresses cell non-independence within samples,
but not a separate concern: whether the per-cell module-score correlation
itself is driven by shared technical structure rather than a specific
biological coupling. This was tested directly on the real NK scRNA data
(`src/topology/count_depth_control.py --real`,
`src/topology/h3_scoring_method_diagnostic.py`,
`src/a100_recompute/run_h3_activation_control_v2.py`). First, a generic
16-gene NK-activation signature (CD69, TNF, XCL1/2, CCL3/4/5, CSF2, IL2RA,
ICOS, TNFSF10, FASLG, CD38, HLA-DRA/B1, MKI67; not present in either module)
explains only part of the raw coupling: linearly partialling it out gives
r=0.251 (P=2.7×10⁻¹¹⁹), and an activation-matched analysis (cells binned into
quintiles of the activation score, so the correlation is computed within cells
of near-identical activation level) gives a consistent mean r=0.242 across
bins — ruling out simple co-activation as the explanation. However, two more
comprehensive controls collapse the coupling substantially further.
**Count-depth/technical covariates**: 47.4% of the protrusion-machinery module
score's variance is explained by `total_counts` and `n_genes_by_counts` alone
(library size varies 236-fold across the 8,310 cells, 427–100,763 counts);
residualizing both modules against these covariates before correlating drops
the raw r from 0.318 to 0.094 (P=1.3×10⁻¹⁷). **Residualizing against the real
scVI latent space** (30 dimensions, the same batch-corrected embedding used
elsewhere in this study) plus the same technical covariates drops it further
to r=0.092 (P=5.5×10⁻¹⁷). Second, whether this residual coupling is
module-*specific* was tested with a permutation test: 10,000 random gene
modules of the same sizes (25 and 11) were drawn from the NK-expressed gene
universe and scored the same way. With the mean-z-score method used
throughout, the observed r=0.318 falls *below* the permutation null mean
(0.728; empirical P=1.0) — but this specific null is inflated by a biased
universe (restricting candidate genes to the 314/22,728 detected in ≥50% of NK
cells excludes most of the protrusion and cytotoxicity module genes
themselves, which are more sparsely detected). Repeating the test with
`scanpy`'s expression-level-matched `score_genes` scoring — the field-standard
control for exactly this confound, applied identically to observed and
permuted modules — gives a materially smaller raw r=0.089 to begin with, and
the observed value still does not exceed the matched permutation null (mean
0.354, 95th percentile 0.560; empirical P=0.97). Full diagnostics are in
`results/tables/h3_scoring_method_diagnostic_summary.md`.

**Conclusion (H3 single-cell).** The single-cell H3 number, while not
attributable to generic co-activation alone, is not distinguishable from the
technical and global-transcriptional-structure background of this dataset,
and should not be reported as an independent replication of the effector arm.
The bulk TCGA-LIHC result is unaffected by this *single-cell* diagnostic (bulk
samples do not carry the extreme per-cell library-size range that drives the
single-cell artifact), but it carries a separate, bulk-specific confound
addressed next. The single-cell pseudoreplication-corrected number is retained
in the main-text table for transparency but flagged as not surviving
technical-confound control.

**Bulk purity control (H3).** In bulk tumor transcriptomes, both the
protrusion-machinery and cytotoxicity-output module scores rise and fall with
the NK-cell fraction of the sample, so a positive protrusion~cytotoxicity
correlation is expected in part from co-varying NK abundance rather than a
within-NK coupling. This was tested directly (`src/topology/bulk_h3_purity_control.py`)
by partialling out a clean NK-lineage abundance proxy — the mean-z of seven NK
lineage/receptor markers (FCGR3A, KLRD1, KLRF1, KLRK1, NCAM1, NCR1, TYROBP)
that appear in *neither* module, avoiding the over-control that would result
from using the full NK signature (whose genes overlap the cytotoxicity
module). In TCGA-LIHC the coupling attenuates by 55%, from zero-order r=0.55
to a partial r=0.25 (P=1.8×10⁻⁶, 95% CI [0.15, 0.34]) — reduced but still
significant. Applying the same control to the two external gastric cohorts
roughly halves it in GSE62254 (partial r=0.23, P=4.6×10⁻⁵, still significant)
and **abolishes it in GSE84437** (partial r=−0.07, P=0.11, 95% CI crossing
zero), the cohort with the strongest zero-order coupling (r=0.62). The
effector-arm coupling is therefore real but **substantially
infiltration-driven**: after NK-fraction adjustment it holds at roughly half
strength in two of three bulk cohorts and does not survive in the third
(`results/tables/bulk_h3_purity_control.tsv`).

**H2 (metabolic arm).** The SM-balance→protrusion coupling is undetectable in
bulk (r=−0.017, P=0.72). After NK-cell isolation, the naive per-cell
correlation was r=+0.030 (P=6×10⁻³), but this treats 8,310 cells as
independent observations from only 9 samples — a severe pseudoreplication
inflation. Applying per-sample Pearson r → Fisher z → DerSimonian-Laird
random-effects meta-analysis, the corrected meta-analytic r=0.029 (P=0.20,
I²=73% substantial heterogeneity). Shared variance is r²=0.0009 and the 95% CI
of the per-sample r spans [−0.059, +0.132].

**H4/H5 (physical topology).** The topology-permissive→dysfunction
relationship carries the wrong sign at both resolutions before and after
correction (bulk r=+0.311; single-cell NK vs HAVCR2 corrected r=+0.036,
P=9.1×10⁻⁴, I²=0%). Intratumoral NK cells show *higher* protrusion-machinery
transcription than normal-tissue NK (Δ=+0.142, p=3.0×10⁻⁹¹) — opposite to the
physical collapse — even though cytotoxic output is correctly reduced
(Δ=−0.141, p=5.9×10⁻⁵²). The reduced cytotoxic output of intratumoral NK is
independently corroborated by a single-cell HCC study reporting that
intratumoral, relative to peritumoral, NK cells upregulate
inhibitory-checkpoint/exhaustion programs and downregulate cytotoxicity
pathways [46]; the opposing protrusion-transcription result is this study's
own observation.

Leave-one-sample-out sensitivity analysis
(`results/tables/h3_leave_one_sample_out.tsv`) shows the H3 pooled estimate
itself is stable to removing any single sample (pooled r range 0.275–0.350
across the 9 leave-one-out re-analyses, all 95% CIs excluding 0) — so the
pseudoreplication correction is not an artifact of one outlier sample; this is
a separate question from the technical-confound diagnostics above.

## S3.4 — Domain-knowledge baseline comparison and tooling roadmap

To test whether the GNN adds value over standard bioinformatics approaches
that do not use graph structure, two additional baselines were evaluated on
the identical folds: (i) an *NK-marker signature* baseline — logistic
regression on the mean expression of 8 canonical NK markers (NCAM1, KLRD1,
NKG7, GNLY, KLRF1, EOMES, NCR1, FCGR3A), conceptually the simplest possible
deconvolution proxy; and (ii) an *SST-module signature* baseline — logistic
regression on the 7 SST-axis module scores computed directly on bulk
expression, capturing the "use the anchor paper's gene modules without
building a graph" approach.

The SST-module baseline attains AUROC 0.904±0.029 and MCC 0.619±0.088 — not
significantly different from the GNN (ΔAUROC −0.046, t-test p=0.11; ΔMCC
−0.087, p=0.28). The simpler NK-marker baseline (AUROC 0.849, MCC 0.503) is
significantly below the GNN on both AUROC (Δ−0.101, p=0.032) and AUPRC
(Δ−0.177, p=0.023). A full comparison with established deconvolution tools
(CIBERSORTx, quanTIseq) and phenotype-genotype association methods (Scissor)
would further strengthen this conclusion; a detailed reproduction roadmap is
provided in `03_supplementary/CIBERSORTx_quanTIseq_Scissor_roadmap.md` (these
tools require an R environment not available in the current local setup).

## S3.5 — Candidate-target orthogonal validation (full narrative)

**Trivial baseline comparison.** To quantify whether the five-dimension
scoring adds information beyond simply listing every gene named in the Zheng
2023 anchor paper, it was compared against a trivial baseline that ranks genes
solely by mechanism-card membership (`in_sst_axis`=1.0) plus gold-standard
literature support (`gold_standard`=0.5), using no expression data. Within the
37 tumor-intrinsic candidates, the trivial baseline shows moderate correlation
with the five-dimension scoring (Spearman ρ=0.54, p=5×10⁻⁴): the rankings are
linked but not identical. The five-dimension scoring contributes incremental
value in three ways: (i) it re-orders within the SST set by quantitative tumor
specificity — e.g., PHGDH (log2FC +0.059) ranks above PSAT1/SHMT1 despite both
being SST members; (ii) it surfaces 12 non-SST genes that the trivial baseline
would rank at positions 26–37 — most notably COL1A1/COL1A2 (log2FC ~0.15, rank
6–7 vs trivial rank 27–28), NECTIN2 (log2FC +0.11, rank 9 vs 29), CA9, ERBB2,
FGFR2, and MET — gastric-cancer-relevant targets with measurable tumor-cell
signal that a "read the anchor paper" approach would entirely miss; and (iii)
it demotes SST-member genes with near-zero tumor specificity to the bottom of
the list — SPTLC1/3, WASF1/3, DIAPH3 (log2FC <0.01, ranks 33–37 vs trivial
ranks 21–25). Full comparison tables:
`results/tables/trivial_baseline_comparison.tsv`,
`trivial_baseline_overlap.tsv`.

**DepMap essentiality note.** CERES values are from a real query of DepMap
Public 26Q1 (`CRISPRGeneEffect.csv`, `Model.csv`; the current release at the
time of analysis, superseding the 25Q2 named in Methods and an earlier 24Q2
figshare snapshot), restricted to 35 cell lines with an Oncotree subtype
containing "Stomach." Of the four genes qualifying on the NK-state DE filter,
NT5E is genuinely non-essential (CERES>0); SGMS2, SMPD1, and SMPD3 fall in the
"weakly essential" band, common and largely unremarkable for most genes in
most cell lines. All four are clearly distinguishable from the genuinely
essential genes excluded (RAC1, MTHFD1; CERES<−0.5) and from ERBB2
(CERES=−0.36). The validation language is accordingly "does not fall in the
pan-essential range that would confound an immune-evasion interpretation with
a cell-viability effect," not "passes a non-essentiality filter."

**NT5E self-referentiality caveat.** A gene-set separation audit
(`results/tables/geneset_separation_audit_summary.md`) found that NT5E is
itself one of the ten genes constituting `NK_DYSFUNCTION_GENES`, the marker
set whose z-score (net of the cytotoxicity score) defines the
NK-hot-dysfunctional label used for this DE test. NT5E's "upregulated in
dysfunctional tumors" result is therefore partly self-referential — higher
NT5E expression directly contributes to a tumor being labeled dysfunctional in
the first place — and should not be read as independent evidence to the same
degree as SGMS2, SMPD1, and SMPD3, none of which overlap any NK-state
label-defining gene set.

**External replication of the NK-state DE test.** The identical NK-state
labeling rule (`src/immune_scoring/nk_scores.py`) was re-run independently on
the two external gastric bulk cohorts — GSE62254 (ACRG, n=300;
NK-hot-cytotoxic n=105, NK-hot-dysfunctional n=13) and GSE84437 (n=483;
NK-hot-cytotoxic n=122, NK-hot-dysfunctional n=29) — and the direction of each
flagged gene's dysfunctional-vs-cytotoxic difference compared against the
TCGA-STAD result (`results/tables/nk_state_de_external_replication_summary.md`,
`nk_state_de_external_concordance.tsv`). Fourteen of 16 gene/cohort
comparisons (88%) are directionally concordant. The five downgraded
serine-pathway genes and RAC1 — none of which overlap any NK-state
label-defining gene set, so this check is not circular — are **DOWN in
dysfunctional tumors in both external cohorts**, independently confirming
their exclusion from the priority tier. NT5E is **UP** in both external
cohorts, concordant with TCGA-STAD, but subject to the same label-overlap
caveat above. **SGMS2 is discordant in both external cohorts** (DOWN in
dysfunctional tumors, opposite the UP direction and FDR=0.007 reported for
TCGA-STAD): the NK-state DE signal that helped qualify SGMS2 does not
externally replicate and should be downweighted; SGMS2's Tier-1 status rests
primarily on mechanistic privilege and DepMap non-pan-essentiality, not on the
NK-state DE direction.

**Serine-pathway downgrade.** Five serine-pathway enzymes — PSPH, SHMT1,
SHMT2, MTHFD1L, MTHFD1 — show the *opposite* expression pattern: significantly
*lower* in NK-dysfunctional tumors (log2FC range −0.09 to −0.13, all p<0.02),
arguing against high-priority inclusion. MTHFD1 is additionally pan-essential
in vitro (CERES=−0.53); PSPH, SHMT1, SHMT2, MTHFD1L are weakly essential
(CERES −0.10 to −0.16); PSAT1 is moderately essential (CERES=−0.27). RAC1
(CERES=−0.74) is clearly pan-essential and should be interpreted as a
cell-viability rather than immune-evasion target. ERBB2 (CERES=−0.36,
moderately essential, and separately an FDA-approved HER2 target) should
likewise not be treated as a novel immune-evasion candidate. Full results:
`results/tables/target_validation_v2_merged.tsv`.

## S3.7 — Multi-view audit detail

After removing all label-defining genes, the external LIHC ensemble
AUROC/AUPRC was 0.922/0.846 for the no-graph expression model, 0.906/0.822 for
merged SVD, 0.914/0.832 for uniform fusion, and 0.912/0.826 for learned
fusion. Learned fusion was worse than no graph for both AUROC (Δ=−0.0097, 95%
bootstrap CI −0.0185 to −0.0024) and AUPRC (Δ=−0.0205, −0.0381 to −0.0063),
and indistinguishable from uniform fusion. In the unmasked sensitivity
analysis the same learned model reached 0.971/0.936, showing why the
label-overlap control is load-bearing rather than optional. No view was
top-ranked in at least 8/10 seeds (MSigDB highest in 7/10), and every
leave-one-view-out contrast failed the pre-registered joint weight/CI
contribution gate. A topology-specific calibration permuted the gene labels
of each mechanism view 1,000 times while fixing the other five views: both
authored views imposed their intended geometry (observed vs null mean
coupling: `metabolic_crosstalk` 0.225 vs 0.031; `sm_topology_axis` 0.169 vs
0.079; both empirical P=0.001) — a structural check, not biological or
predictive validation.

## S3.8 — Real-data comparative recoverability atlas (exploratory, gate not met)

**Scope note.** This section is an exploratory scan, separate from the two-card
(Zheng serine-SM; TGFβ→SMAD) result that anchors the main text's
transcriptional-reach claim (§3.7, §4.1, §4.2). It applies the same style of
comparison to all four registered mechanism cards — including the adenosine
and MICA/B cards, which are validated only structurally (schema compliance,
synthetic-mode ingestion) and have not been run end-to-end elsewhere in this
study — across four verified bulk cohorts plus additional orthogonal-modality
sources, pre-registering a stricter three-card/two-cohort/direct-modality gate
for treating it as confirmatory evidence. That gate is **not met**, and the
verdict below should be read as a comparative map of what has and has not been
checked, not as evidence either for or against the reach boundary reported in
the main text.

Applying the same pre-specified module comparisons to all four registered
cards in four verified human bulk cohorts (TCGA-STAD, TCGA-LIHC, GSE62254 and
GSE84437; Fig. S2), the two SM-topology downstream comparisons were
directionally concordant with BH-FDR <0.05 in all cohorts. The NKG2D
recognition-to-cytotoxicity comparison likewise recovered, whereas its
tumor-shedding comparison did not. Neither the pre-specified adenosine nor the
TGFβ inhibitory comparison recovered in the expected direction. Direct
metabolomics and sample-level protein matrices remain unavailable and are
recorded as `not_measured`. Spatial transcriptomics is available only as an
explicitly exploratory subset of four independently verified GSE251950 GSM
sections: strict first-order Visium spot-grid CAF-NK module adjacency is
heterogeneous (three calibrated positive sections and one null section) and
is not a cell-contact measurement. It therefore does not satisfy the
pre-registered cross-mechanism direct-modality gate. The pre-registered
three-card/two-cohort/direct-modality gate consequently remains
`comparative_atlas_only`.

## S4.3 — Limitations (full text)

The main text lists a condensed set of limitations. The complete 17-item list
is preserved here for reference; each restates a point already summarized
above or in the corresponding main-text section.

1. **Transcriptional proxy ≠ physical topology.** Gene expression captures the
   molecular machinery and capacity for the axis, not the physical membrane
   phenotype itself; every claim is bounded by "transcriptional program
   permissive-of/associated-with."
2. **Serine/SM crosstalk is a metabolite-level effect.** Transcription
   captures enzyme abundance, not flux; direct measurement requires
   metabolomics or single-cell mass spectrometry. Now empirically supported by
   the H2 non-significance after correction.
3. **Data availability constraints.** The anchor paper did not deposit
   transcriptomic data; the liver positive control is an independent
   validation, not a direct replication.
4. **No experimental validation.** All targets are computationally
   prioritized; recommended assays bridge to a companion experimental program.
5. **Single-cell pseudoreplication.** 8,310 NK cells from 9 samples; corrected
   via per-sample meta-analysis, but I² up to 96% for H3. Leave-one-sample-out
   analysis shows the pooled H3 estimate is stable to removing any single
   sample (range 0.275–0.350, all CIs excluding 0) — a separate question from
   technical confounding (item 6).
6. **Single-cell module-score correlations require count-depth and
   module-permutation controls beyond pseudoreplication correction** — see
   S3.2 above; H2 and H4's already-null verdicts were unchanged by these
   diagnostics.
7. **NK subtype resolution.** Marker-threshold annotation depends on reference
   atlas quality; rare/absent populations may be misclassified.
8. **The graph model does not outperform top tabular or signature baselines
   on accuracy**, confirmed by two independent stricter audits (§3.4 domain
   baselines; §3.7 label-masked multi-view audit). An exploratory single-edge
   cross-cohort transfer test (`metabolic_crosstalk`, STAD→LIHC MCC) proved
   highly seed-count-sensitive (3 seeds: no difference; 10 seeds: edge
   significantly improves transfer) and no directional claim is drawn from it.
9. **Residual NK bias in the tumor-intrinsic candidate pool.** 17 of 37
   candidates retain NK-side mechanism-card module annotation; the pool is
   "genes with non-zero malignant-cell signal that mechanistically intersect
   the SST axis," not a clean tumor-exclusive set.
10. **Candidate atlas omits intracellular/TF-level regulators** (e.g. the
    CREM/PKA–CREB axis [27]) — a natural extension of the mechanism-card
    modules.
11. **No clinical-outcome anchor.** NK states are linked to scRNA-defined
    labels, not survival or therapy response.
12. **Multi-card validation covers two of four cards end-to-end**; the
    adenosine and MICA/B cards remain validated only structurally. The
    TGFβ card's strong H3 partly reflects that its "activating receptor"
    module is itself part of the NK activation program — the informative
    cross-card signal is the *failure* of H2/H4, not the magnitude of H3.
13. **The bulk effector coupling is substantially NK-abundance-driven** — see
    S3.2 bulk purity control above. The NK-lineage proxy is a coarse,
    transcriptome-derived estimate rather than a deconvolution or
    cell-sorted fraction; a deconvolution-based estimate would refine it.
14. **DepMap essentiality uses the real, current 26Q1 release** (see S3.5
    DepMap note above); the Cloudflare bot-check on DepMap's interactive
    portal prevented non-interactive access to intermediate releases.
15. **NK-state DE is underpowered for the dysfunctional group** (TCGA-STAD
    n=20 vs n=134); external replication (S3.5 above) gives 88% directional
    concordance but a discordant result for SGMS2.
16. **The NK-state classification target partly overlaps the classifier's own
    input features** — a gene-set separation audit found the classification
    label is a thresholding rule on marker genes (NKG7, GNLY, GZMB, PRF1,
    IFNG, among others) also present in the input expression vector, so
    high classification accuracy partly reflects recovering a label from the
    genes that define it (§3.4 should be read accordingly). The SST-axis
    modules used for H1–H5 do not reference this label at all.
17. **No hard per-cell QC threshold was applied to the real scRNA-seq data.**
    All 166,829 concatenated cells were retained; technical variance is
    addressed downstream via count-depth residualization (item 6) rather than
    upfront filtering.
