# Target Validation v2 ¡ª Integrated Evidence

## Dimension 1: NK-state differential expression
Tests whether each of the 37 genes is differentially expressed between
TCGA-STAD tumors with active NK killing (NK-hot-cytotoxic) vs suppressed
NK function (NK-hot-dysfunctional, NK-cold/excluded).

**A genuine immune-evasion target should be UP in dysfunctional/cold tumors**
(the tumor expresses it to suppress NK).

## Dimension 2: DepMap CRISPR essentiality
CERES dependency scores from genome-wide CRISPR KO screens in gastric cancer
cell lines grown IN VITRO (no immune cells).

**A good immune-evasion target should be NON-ESSENTIAL in vitro**
(its value is specifically in the immune-microenvironment context).
A gene that is pan-essential (CERES < -0.5) is NOT a good immune-evasion
target: inhibiting it would kill tumor cells regardless of NK.

## Integrated evidence score
- +3: non-essential in vitro (strong immune-target signal)
- +2: significantly UP in dysfunctional tumors (evasion pattern)
- +1: trending UP or weakly non-essential
-  0: no clear evidence either way
- -1: weakly essential or DOWN in dysfunctional
- -2: pan-essential in vitro (NOT an immune target)

## Interpretation
Genes with score >= 3 have independent, orthogonal evidence (NK-state DE +
DepMap) supporting their role as tumor-intrinsic immune-evasion targets.
Genes with score <= -1 may be false positives in the target list.
