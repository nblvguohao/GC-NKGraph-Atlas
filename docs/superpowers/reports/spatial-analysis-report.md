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

The endpoint is the mean cross-module score over an undirected, strict
first-order Visium hex-grid graph built from the supplied
`array_row`/`array_col` labels. Legal offsets are `(0, +/-2)` and
`(+/-1, +/-1)` only: missing spots stay missing and are never bridged with kNN
edges. It is a **spot-module spatial adjacency** metric. It does not identify
CAF or NK cells, estimate cell-to-cell distances, or establish cellular contact.

For each section, a one-sided coordinate-label permutation p value uses 1,000
random permutations of the NK module score over the fixed supplied spot grid;
CAF scores and grid edges remain fixed. This is structural calibration only.

## Real-data results

| GSM | In-tissue spots | Grid edges | Observed CAF--NK spot-module adjacency | Permutation p value |
|---|---:|---:|---:|---:|
| GSM7990474 | 4,274 | 12,474 | 0.04922 | 0.0010 |
| GSM7990476 | 4,244 | 12,294 | 0.05505 | 0.0010 |
| GSM7990479 | 2,483 | 7,182 | 0.04738 | 0.8442 |
| GSM7990481 | 3,037 | 8,768 | 0.03990 | 0.0010 |

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
  The scope label must remain `exploratory_four_verified_per_gsm_subset` in manuscript and
  supplementary material until that changes.
- No synthetic data, RNA proxy for an unavailable modality, or inferred
  cell-contact metric was used.

## Reviewer-remediation addendum

The formal run now validates every archive against the real-data manifest and
its SHA-256 before opening the TAR. The manifest contains the official URL,
local path, content length, digest, `available` status, and explicit
`exploratory_four_verified_per_gsm_subset` scope for each GSM. It also rejects
duplicate expression barcodes and duplicate in-tissue array-grid coordinates.

The 1,000-permutation null summaries and BH-FDR across the four exploratory
sections are:

| GSM | Null mean | Null SD | Raw p | BH-FDR |
|---|---:|---:|---:|---:|
| GSM7990474 | 0.048925 | 0.0000382 | 0.0010 | 0.0013 |
| GSM7990476 | 0.053993 | 0.0000522 | 0.0010 | 0.0013 |
| GSM7990479 | 0.047491 | 0.0001162 | 0.8442 | 0.8442 |
| GSM7990481 | 0.039298 | 0.0000744 | 0.0010 | 0.0013 |

These are exploratory subset-calibration results, not a replication-level
claim. The three positive sections and one null section reinforce the need to
report section-level heterogeneity rather than a universal CAF--NK spatial
exclusion effect.

## Final manuscript and atlas scope

Section 3.8 of both the submission Markdown and TeX manuscript now separates
the modalities explicitly: metabolomics and sample-level protein matrices are
unavailable (`not_measured`); spatial evidence exists only for the exploratory
four-verified-GSM subset above. The latter is heterogeneous spot-module
adjacency rather than a cell-contact measurement and cannot satisfy the
pre-registered cross-mechanism gate. The cross-mechanism verdict is therefore
unchanged: `comparative_atlas_only`.

The source manifest was regenerated from the YAML configuration and now retains
the pending outer GSE251950 archive alongside the four `available` per-GSM
assets. Targeted tests passed (13 total), and the updated TeX manuscript was
compiled twice; the rendered section 3.8 continuation was visually checked in
the resulting 23-page PDF.
