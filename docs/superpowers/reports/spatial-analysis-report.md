# GSE251950 Visium spot-module adjacency report

**Date:** 2026-07-17 (Asia/Shanghai)  
**Data scope:** exactly four independently verified official GEO per-GSM archives
from GSE251950: GSM7990474, GSM7990476, GSM7990479, and GSM7990481. This is
not a claim about the six unrecovered/unverified sections or the full study.

## Method

The analysis reads the supplied Matrix Market expression matrix, barcode table,
feature table, and `tissue_positions_list.csv` directly from each verified
archive. It calculates the `caf_ecm_program` and `nk_cytolytic_machinery`
scores specified in the active TGF-beta mechanism card. All genes in both
modules were present in every archive (feature coverage 1.0).

The endpoint is the mean cross-module score over an undirected six-nearest
neighbour graph built from the supplied Visium `array_row`/`array_col` spot-grid
labels. It is a **spot-module spatial adjacency** metric. It does not identify
CAF or NK cells, estimate cell-to-cell distances, or establish cellular contact.

For each section, a one-sided coordinate-label permutation p value uses 1,000
random permutations of the NK module score over the fixed supplied spot grid;
CAF scores and grid edges remain fixed. This is structural calibration only.

## Real-data results

| GSM | In-tissue spots | Grid edges | Observed CAF--NK spot-module adjacency | Permutation p value |
|---|---:|---:|---:|---:|
| GSM7990474 | 4,274 | 14,313 | 0.04925 | 0.0010 |
| GSM7990476 | 4,244 | 14,250 | 0.05481 | 0.0010 |
| GSM7990479 | 2,483 | 8,396 | 0.04707 | 1.0000 |
| GSM7990481 | 3,037 | 10,194 | 0.03969 | 0.0010 |

The direct-modality submission table is
`submission_bundle_BiB/03_supplementary/tables/recoverability_direct_modality.tsv`.
It records all four spatial rows as `measured`, while preserving unavailable
metabolomics and protein endpoints as `not_measured`.

## Tests and reproducibility

- `python -m pytest tests/test_recoverability_modalities.py -q` — 4 passed.
- Tests used the real four archives where available; they require exactly that
  accession set and check measured status, nonzero feature coverage, spot-grid
  edges, and a bounded permutation p value.
- Production command:

```powershell
python src/interpretation/run_recoverability_modalities.py --n-permutations 1000 --seed 20260717
```

## Interpretation limits and concerns

- Three sections have positive adjacency relative to the fixed-grid label null,
  whereas GSM7990479 does not. The heterogeneity precludes a universal spatial
  TGF-beta/CAF--NK exclusion claim from this four-section subset.
- The module scores are expression programs at mixed Visium spot resolution;
  they are not cell-type assignments, TGF-beta protein measurements, or
  phospho-SMAD evidence.
- Six GSE251950 sections have not been independently recovered and validated.
  The scope label must remain `four_verified_per_gsm_subset` in manuscript and
  supplementary material until that changes.
- No synthetic data, RNA proxy for an unavailable modality, or inferred
  cell-contact metric was used.
