# H3 scoring-method diagnostic - is the real-data collapse an artifact?

## Question
On real scRNA data (8,310 NK cells, gc_nk_subset_remote.h5ad), the manuscript's
mean-zscore H3 correlation (protrusion~cytotoxicity, r=0.318) collapses under
count-depth residualization (r=0.09) and under a module-membership permutation
test (null mean r=0.73, observed r=0.318 below it, empirical P=1.0). Before
reporting this as a real finding, we checked whether it is an artifact of the
scoring method or the permutation test's null construction.

## Finding A: the permutation test's "universe" is biased
The permutation null was built by drawing random genes from those detected in
>=50% of NK cells - only 314/22728 genes, with median
detection rate 0.665. The real
protrusion module (median detection 0.380) and cytotoxicity
module (median detection 0.332) are both mostly *below*
this threshold (5/24 and
1/11 genes qualify). The narrow,
atypically highly-expressed 314-gene universe is not a fair null population for
these modules, and likely inflates the specific null-mean value of 0.73.

## Finding B: the qualitative result replicates with the field-standard method
Using scanpy's `sc.tl.score_genes` (expression-level-matched control-gene
scoring, the standard method for exactly this technical confound) instead of
mean-zscore:
- Raw H3: r=0.0886 (p=6.01e-16) - already far below the
  manuscript's mean-zscore r=0.318, indicating mean-zscore itself inflates the
  raw single-cell correlation.
- After count+n_genes residualization: r=0.0411 (p=1.80e-04),
  R2_tech = 0.280 (protrusion) / 0.120 (cytotoxicity).
- A matched permutation null (same score_genes method applied to both real and
  randomly-drawn modules, N=200): null mean = 0.3539,
  95th percentile = 0.5596. The observed r
  (0.0886) does **not** exceed this null (empirical P =
  0.9700).

## Verdict
The specific magnitude of the null distribution (0.73 vs 0.35) is
method-dependent and the naive-universe permutation test overstates it. But the
qualitative conclusion is **not** an artifact: with the field-standard,
expression-matched scoring method, the observed single-cell
protrusion~cytotoxicity correlation still does not exceed what a randomly drawn
gene-module pair of the same sizes would produce, and a large fraction of its
variance (28%/12%) is explained by technical library-size
covariates. This is a real methodological finding, not a scoring bug.

## Implication for the manuscript
The single-cell H3 pseudoreplication-corrected number (corrected r=0.313,
P=3.9e-8) addresses cell non-independence within samples, but does not address
this separate technical-confound / module-specificity problem. The single-cell
H3 result should not be reported as an independent, specific replication of the
effector arm. The bulk TCGA-LIHC result (r=0.55, deconvolution-based, not
subject to per-cell dropout/depth artifacts in the same way) is unaffected by
this diagnostic and remains the primary evidence for the effector-arm claim.
