# NK-state DE external replication (S3 optional item)

## Method
The TCGA-STAD NK-state DE test (Sec 3.5, Table 4) compares the 37
tumor-intrinsic candidate genes between NK-hot-cytotoxic (n=134) and
NK-hot-dysfunctional (n=20, underpowered per Limitations) tumors. Rather
than acquiring a new external dataset, we re-ran the identical labeling rule
(`src/immune_scoring/nk_scores.py`) on the two bulk gastric cohorts already
used for the effector-arm external validation (Sec 3.3): GSE62254 (ACRG,
n=300, GPL570) and GSE84437 (n=483, GPL6947).

## Cohort NK-state distributions
- **GSE62254**: NK-hot-cytotoxic n=105, NK-hot-dysfunctional n=13
- **GSE84437**: NK-hot-cytotoxic n=122, NK-hot-dysfunctional n=29

## Directional concordance with the TCGA-STAD internal test
For the genes flagged in Table 4 (Tier 1 UP genes SGMS2/NT5E; Tier 3/X DOWN
genes PSPH/SHMT1/SHMT2/MTHFD1L/MTHFD1/RAC1), we compare the direction of the
dysfunctional-vs-cytotoxic difference in each external cohort against the
TCGA-STAD direction. Concordance rate: **14/16**
(88% if n_total>0 else 'n/a').

See `nk_state_de_external_concordance.tsv` for the full per-gene, per-cohort
comparison, and `nk_state_de_external_replication.tsv` for DE results on all
37 genes in both cohorts.

## Caveats
- These external cohorts are microarray (not RNA-seq), and NK-state labels
  are derived independently per cohort from cohort-specific score
  distributions (median-based thresholds), not transferred from TCGA-STAD --
  so this is a test of directional consistency, not an exact replication of
  the same patients or platform.
- Sample sizes in the dysfunctional group vary by cohort and may still be
  small; see the per-cohort n's above.
