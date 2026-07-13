# Domain-Method Baseline Comparison ¡ª Roadmap for External Tools

## Purpose
This document describes how to reproduce the domain-method baseline comparison
using CIBERSORTx, quanTIseq, and Scissor, extending the analysis in
`src/baselines/run_domain_baselines.py`.

## Environment requirements
- R >= 4.0
- R packages: `immunedeconv` (quanTIseq), `Scissor`
- CIBERSORTx: Docker image `cibersortx/fractions` or web portal
  (https://cibersortx.stanford.edu)

## Baseline 1: CIBERSORTx NK-fraction baseline

### Rationale
CIBERSORTx deconvolves bulk RNA-seq into immune cell-type fractions using a
signature matrix (LM22). The NK-cell fraction estimate is used as a single
feature to predict NK-hot-cytotoxic state.

### Procedure
```r
# Step 1: Prepare mixture file (TCGA-STAD expression, genes x samples)
# Step 2: Run CIBERSORTx
#   docker run -v $PWD:/data cibersortx/fractions #     --username <token> --token <token> #     --refsample LM22.txt --mixture tcga_stad_mixture.txt #     --perm 100 --QN FALSE
# Step 3: Extract NK fraction column (NK cells resting + NK cells activated)
# Step 4: 5-fold stratified CV, logistic regression on NK fraction -> state
```

### Expected metrics (approximate, based on known TCGA-STAD NK content)
- AUROC: 0.65¨C0.75 (NK fraction alone is a weak predictor of NK activation state)
- This baseline tests whether the GNN adds information beyond simple NK abundance

## Baseline 2: quanTIseq immune-deconvolution baseline

### Rationale
quanTIseq (Finotello et al., *Genome Med* 2019) provides absolute immune cell
fractions using a constrained least-squares approach. Like CIBERSORTx, the NK
fraction serves as the predictor.

### Procedure
```r
library(immunedeconv)
# Load TCGA-STAD expression (TPM-normalized, genes x samples)
expr <- read.table("tcga_stad_expression.tsv", row.names=1, header=TRUE)
# Run quanTIseq
res <- deconvolute(expr, method="quantiseq")
# Extract NK cell fraction
nk_fraction <- as.numeric(res["NK cell", ])
```

### Expected comparison
- NK fraction per sample -> logistic regression -> 5-fold CV
- AUROC expected 0.65¨C0.75 (similar to CIBERSORTx)
- Also compute: compare with SST-module baseline (which uses mechanism-specific
  modules rather than generic immune markers)

## Baseline 3: Scissor phenotype-to-genotype baseline

### Rationale
Scissor (Sun et al., *Nat Biotechnol* 2022) links single-cell phenotypes to bulk
clinical variables. The scRNA-defined NK states (hot-cytotoxic / dysfunctional /
cold) are treated as phenotypes, and Scissor identifies which bulk samples
associate with each phenotype.

### Procedure
```r
library(Scissor)
# Load scRNA data: scRNA NK cells (8,310 cells, 9 samples)
#   with phenotype labels (NK-hot-cytotoxic, NK-hot-dysfunctional,
#   NK-cold/excluded, NK-intermediate)
# Load bulk data: TCGA-STAD expression (450 samples)
# Run Scissor:
scissor_output <- Scissor(
  bulk_dataset = tcga_stad_expr,
  sc_dataset = nk_scrna_expr,
  sc_phenotype = nk_phenotype_labels,  # "NK-hot-cytotoxic" = phenotype of interest
  alpha = 0.05,
  family = "binomial"
)
# Scissor identifies "Scissor+" samples (associated with NK-hot-cytotoxic)
# Use Scissor+/- label as predictor
```

### Expected comparison
- Scissor+/- label -> classification performance vs GNN
- Key question: does the GNN's mechanism-structured embedding capture
  information beyond phenotype-genotype association?

## What the current repository already provides

`src/baselines/run_domain_baselines.py` implements two directly runnable
baselines that require no R:
1. **NK-marker signature baseline** ¡ª logistic regression on mean NK marker
   expression (conceptually the simplest possible deconvolution proxy)
2. **SST-module signature baseline** ¡ª logistic regression on 7 SST-axis
   module scores (the closest no-graph alternative to the GNN)

These two baselines capture the core question: "Does the graph add anything
beyond the anchor paper's gene modules scored directly on bulk expression?"

The CIBERSORTx/quanTIseq baselines answer a related but distinct question:
"Does the GNN add anything beyond generic NK abundance estimates from
state-of-the-art deconvolution?"

The Scissor baseline answers: "Does the GNN add anything beyond
phenotype-genotype association from single-cell data?"

All three external comparisons are complementary and would strengthen the
manuscript's "our framework adds value over existing tools" claim. They
require an R environment not available in the current local setup.
