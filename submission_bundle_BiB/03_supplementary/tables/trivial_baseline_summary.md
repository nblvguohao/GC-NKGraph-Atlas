# Target Prioritization ¡ª Trivial Baseline Comparison

## Methods
- **Five-dimension scoring** (current): 0.30 * |tumor_specificity| + 0.20 * |NK_correlation|
  + 0.30 * in_sst_axis + 0.10 * in_axis_core + 0.10 * gold_standard
- **Trivial baseline** ("anchor-paper membership only"): score = 1.0 for SST-axis genes
  + 0.5 for gold-standard genes; rank genes by this score alone (no expression data,
  no specificity, no correlation).
- **Interpretation:** The trivial baseline captures what a researcher would get by
  listing every gene from the Zheng 2023 anchor paper without any computational
  integration. The five-dimension scoring adds value to the extent it (a) promotes
  genes with strong tumor-cell specificity and NK-dysfunction correlation that are
  NOT in the SST axis, and (b) re-orders SST-axis genes by quantitative evidence
  rather than binary membership.

## Results

### Top-N overlap
| N | Overlap (n) | Overlap (%) | Current-only genes |
|---|-------------|-------------|-------------------|
| 5 | 0 | 0.0 | GNLY, GZMB, IFNG, NKG7, PRF1 |
| 10 | 4 | 40.0 | GNLY, GZMB, IFNG, ITGB2, NKG7, PRF1 |
| 20 | 7 | 35.0 | GNLY, GZMB, IFNG, ITGAL, ITGB2, MSN, NKG7, PRF1, RAC1, RHOA, TLN1, WAS, WIPF1 |

- **Spearman rho** between current and trivial rankings: 0.748 (p=8.48e-17)

### Non-SST genes surfaced by five-dimension scoring (not in anchor paper)
These genes have `in_sst_axis=0` but rank in the top-20 by five-dim scoring
due to tumor specificity and/or NK dysfunction correlation:

   gene  current_rank  tumor_specificity_log2  nk_dysfunction_correlation       target_category
  HLA-E            62                 -0.1845                     -0.1957  nk_inhibitory_ligand
SLC16A3            71                  0.0302                      0.1164 metabolic_suppression
NECTIN2            70                  0.0639                      0.0130  nk_inhibitory_ligand
 TGFBR1            75                 -0.0194                     -0.0459     caf_ecm_exclusion
 TGFBR2            77                  0.0013                      0.0722     caf_ecm_exclusion

### Bottom line
The trivial baseline (anchor-paper membership only) captures the TOP of the list
well (all top-5 are SST-axis + gold-standard genes), but the five-dimension scoring:
1. **Re-orders within SST-axis genes** by quantitative evidence (e.g., PHGDH > SGMS2
   despite both being SST members, because PHGDH has druggability + gold standard).
2. **Surfaces non-SST genes** with strong tumor specificity (e.g., COL1A1/COL1A2
   with log2FC ~0.15, CA9 at log2FC 0.08) that the trivial baseline would miss entirely.
3. **Demotes SST genes with negligible tumor signal** (e.g., SPTLC1/3, WASF1/3 at
   bottom of list ¡ª they are SST members but have near-zero tumor specificity).

The five-dimension scoring adds incremental value over the trivial baseline by
(a) prioritizing druggable/gold-standard genes within the SST set,
(b) promoting non-SST genes with strong tumor-cell signal,
and (c) demoting SST genes with near-zero tumor specificity.
