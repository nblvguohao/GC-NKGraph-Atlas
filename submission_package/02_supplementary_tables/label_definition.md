# NK Immune-State Label Definition

## Scores
- NK_infiltration_score: mean z-score of NCAM1, NCR1, KLRD1, KLRK1, KLRF1, NKG7, GNLY, GZMB, PRF1, FCGR3A, XCL1, XCL2, CCL5, IFNG, TYROBP
- NK_cytotoxicity_score: mean z-score of NKG7, GNLY, GZMB, PRF1, IFNG, XCL1, XCL2, CCL5
- NK_dysfunction_score: mean z-score(dysfunction) - mean z-score(cytotoxicity)
- NK_exclusion_score: mean z-score(CAF/ECM) - NK_infiltration_score

## Immune States
- NK-hot-cytotoxic: high infiltration, high cytotoxicity, low dysfunction
- NK-hot-dysfunctional: high infiltration, low cytotoxicity, high dysfunction
- NK-cold/excluded: low infiltration, high exclusion
- NK-intermediate: remaining samples

## Thresholds
Computed on training primary dataset (TCGA-STAD) as median values.
External validation uses saved thresholds from training.

## Caveats
- Gene signatures are literature-derived starting sets.
- scRNA-derived cell-type attribution may refine scores.
- These are transcriptional proxies, not functional measurements.
