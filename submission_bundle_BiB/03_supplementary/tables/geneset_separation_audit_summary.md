# Gene-set separation audit (M4)

## Question
Does the NK immune-state label (used as the classifier target in §3.4, and to
define the NK-hot-cytotoxic / NK-hot-dysfunctional groups for the NK-state DE
test in §3.5) share genes with (a) the SST-axis modules used for the H1-H5
mechanism tests, or (b) the 37 tumor-intrinsic candidate targets?

## Finding 1: NK-state label vs SST-axis cytotoxicity-outcome module
`NK_CYTOTOXICITY_GENES` (used both to build `NK_cytotoxicity_score` and,
subtracted, `NK_dysfunction_score` -- i.e. it defines the classifier's target
label) overlaps with the SST-axis `nk_synapse_cytotoxicity_outcome` module
(used for the H3/H5 mechanism tests) in 5 genes:
['GNLY', 'GZMB', 'IFNG', 'NKG7', 'PRF1'].

This overlap does **not** affect the H3/H5 correlation tests themselves (those
operate directly on expression-derived module scores and never reference the
NK-state label). It **does** mean that the NK-state classification task
(§3.4) has partial "label leakage" by construction: the target label is a
thresholding rule on a handful of marker genes (NKG7, GNLY, GZMB, PRF1, IFNG
among others) that are also present, unmodified, in the full expression vector
x used as classifier input for every baseline and the GNN. This is a known
and generally-accepted property of marker-defined phenotype classification
(the label is a deterministic function of a subset of the input features), not
a coding error -- but it explains why even simple baselines (e.g. the 8-gene
NK-marker signature baseline, AUROC=0.849) perform well, and it means §3.4
accuracy numbers should be read as "how well can a model recover a
marker-gene-defined label from the transcriptome that contains those same
markers," not as an unconstrained prediction task.

## Finding 2: NK-state label vs checkpoint_link module
`NK_DYSFUNCTION_GENES` (used to build `NK_dysfunction_score`, part of the
label rule) overlaps the SST-axis `checkpoint_link` module ({HAVCR2}) in
1 genes: ['HAVCR2'].

## Finding 3: 37 tumor-intrinsic candidates vs NK-state label genes
The 37 tumor-intrinsic candidate targets (Table 4, §3.5) overlap the NK-state
label-defining gene sets (infiltration + cytotoxicity + dysfunction) in
1 genes: ['NT5E'].
This is a non-trivial overlap and means the NK-state DE test for these genes is partly self-referential.

## Full pairwise overlap table
See `geneset_separation_audit.tsv` for all pairwise Jaccard overlaps among the
NK-state label components, the seven SST-axis modules, and the 37-gene target
list.

## Recommendation
State explicitly in Methods (§2.4/§2.8) that the NK-state classification
task's target label is partly constructed from marker genes present in the
classifier's own input (Finding 1/2). For the NK-state DE test (§3.5), flag
that NT5E -- one of the four Tier-1 orthogonally-validated candidates -- is
itself a constituent of `NK_DYSFUNCTION_GENES` (Finding 3): its "upregulated
in NK-hot-dysfunctional tumors" result is partly self-referential, since
higher NT5E expression directly contributes to a tumor being labeled
dysfunctional in the first place. The other three tested genes among the
37-candidate pool have no overlap with any label-defining gene set, so their
NK-state DE results are not subject to this concern; NT5E's Tier-1 status
should be downgraded to reflect the circularity, or corroborated with an
NK-state definition that excludes NT5E from the label rule.
