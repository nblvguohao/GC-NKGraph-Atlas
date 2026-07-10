# Nature-style pre-submission reviewer assessment (EN)

> Generated 2026-07-10 via the `nature-reviewer` skill against
> `manuscript/main_manuscript.md` (v0.5). Domain: computational immuno-oncology /
> bioinformatics. None of the skill's domain gates cover this domain, so technical
> soundness is assessed from the manuscript's own evidence chain plus general axes.
> Stated target journal is *Briefings in Bioinformatics*; evaluated against
> Nature-style axes as the skill requires, with fit flagged honestly and **no**
> editorial decision asserted.

## Review setup
- **Input scope:** Full main manuscript (abstract, key points, intro, methods, results, discussion, conclusion, tables, figure legends, references). Figures as legends only; supplementary tables and code referenced but not inspected; no code run.
- **Assessment boundary:** Grounded in manuscript text and reported numerics. Cannot verify underlying computations, figure panels, supplementary tables, or that statistics match deposited outputs. Single-cell claims rest on one scRNA dataset (GSE246662).
- **Shared claim summary:** How much of Zheng et al. 2023's serine→sphingomyelin→membrane-topology→cytotoxicity NK-evasion axis is recoverable from public transcriptomes. Via a reusable "mechanism-card" YAML abstraction, a mechanism-grounded heterogeneous gene graph, and a GNN, the authors report a **scoping map**: the protrusion→cytotoxicity effector arm recovers robustly (bulk LIHC r=0.55; single-cell r=0.32) and generalizes to gastric cancer (external cohorts r=0.42/0.62); the upstream metabolic arm is negligible (single-cell r=0.03, r²<0.001); the physical topology phenotype is NOT captured by machinery transcription (wrong sign; intratumoral NK show higher, not lower, protrusion transcripts). Output: 37 putative tumor-intrinsic candidate targets led by druggable serine/SM enzymes.
- **Visible evidence base:** Four bulk cohorts (TCGA-LIHC 423, TCGA-STAD 450, GSE62254 300, GSE84437 483) + one scRNA (166,829 cells; 8,310 NK). Pre-registered H1–H5, multi-resolution per-hypothesis outcomes, FDR correction, 95% CIs, an activation partial-correlation control, a de-circularization audit, an edge-type ablation, six-method baseline comparison with paired tests.
- **Missing materials affecting confidence:** Figure panels; supplementary tables; independent single-cell replication; any wet-lab validation (explicitly none); ORCID/author placeholders; the multi-card (adenosine/TGFβ) results alluded to but not shown.

## Reviewer 1 — emphasis: originality, significance, breadth of interest
- **Overall:** Methodologically careful and unusually honest; the "map of transcriptional reach" framing is appealing, but the positive findings largely re-confirm known NK biology while the genuinely novel tests return negative/negligible results. Credible as a bioinformatics methods+scoping contribution; below a Nature-style bar for a decisive, broadly-interesting conceptual advance.
- **Who would be interested:** Computational immuno-oncology and NK-biology groups extending physical/metabolic mechanisms to cohort scale; single-cell methodologists; target-discovery teams. Not obviously of broad cross-disciplinary Nature interest.
- **Major strengths:** (1) Honest reframing of a would-be validation into an explicit scoping result; (2) pre-registration with an explicit recovery criterion, reporting a *failed* criterion as a structured result; (3) consistent claim-boundary discipline (transcription ≠ physical topology).
- **Major concerns:** (1) The advance may re-derive known biology — protrusion~cytotoxicity co-expression (H3) and reduced intratumoral NK cytotoxicity (H5) are established properties of NK activation/dysfunction; (2) the novel tests (metabolic arm, topology) failed; (3) framework novelty is largely an engineering (YAML) abstraction with one card exercised.
- **Technical failings to address:** Prove the effector coupling is mechanism-specific rather than generic activation covariance (negative-control module pairs); show a second card yields a non-trivial result or drop the "engine" claim.
- **Recommendation posture:** Below broad-significance threshold as framed; strong fit for a specialist computational-biology venue (consistent with the BiB target). Not a settled editorial matter.

## Reviewer 2 — emphasis: technical soundness, statistics, de-circularization, ablation
- **Overall:** Statistically disciplined and self-critical, but several load-bearing conclusions rest on effect sizes so small their biological meaning is questionable, and two supporting analyses (target ranking, edge ablation) are closer to tautological than acknowledged.
- **Major strengths:** FDR correction; 95% CIs; explicit r²; the H3 activation partial-correlation control; paired Wilcoxon/t-tests; a self-reported 46% NK-side de-circularization audit.
- **Major concerns:** (1) Near-zero fold-changes drive the target list (log2FC +0.059/+0.038/+0.010/+0.001); (2) the edge ablation is near-tautological (SST edges are constructed from the same mechanism that defines H1/H2) and shows no external predictive value; (3) H2 (r²=0.0009) is honestly labeled negligible yet still cited in §4.1 as "validating the edge's design premise"; (4) all single-cell claims rest on one scRNA dataset with marker-threshold annotation.
- **Technical failings to address:** Report DE statistics/FDR/magnitude for the 37 candidates; replace the ablation with an out-of-sample predictive test; replicate key single-cell couplings in an independent scRNA cohort; reconcile the §4.1 "validates" statement with r²<0.001.
- **Recommendation posture:** Major technical revision before the target-list and edge-value claims stand; transparency strengths are real.

## Reviewer 3 — emphasis: biological interpretation, single-cell methodology, readability
- **Overall:** Clearly written and biologically literate. Main risk: the most-quoted numbers (the "recovered" arm) may measure NK activation covariance rather than the specific axis; the target list mixes NK-side and tumor-side genes despite de-circularization.
- **Major strengths:** Transparent seven-module operationalization; explicit separation of the axis-confirmation panel (readout) from tumor-intrinsic candidates (targets); thorough limitations.
- **Major concerns:** (1) Effector-coupling interpretation — both modules are hallmarks of activated cytotoxic NK; the 16-gene activation control may not span the shared manifold (partial r still 0.286); (2) marker-threshold NK annotation; NKT/tissue-resident contamination could drive the "higher intratumoral protrusion transcription" that anchors the topology conclusion — run reference-based mapping (scANVI/scArches); (3) residual contamination — 17/37 candidates are NK-side; the abstract should carry this caveat.
- **Technical failings to address:** Specificity control for the effector coupling; reference-mapping sensitivity for the NK compartment (especially H5); align the abstract's target-list claim with the audited contamination.
- **Recommendation posture:** Major revision on interpretation and annotation controls; presentation quality high.

## Cross-review synthesis
- **Consensus strengths:** Rare statistical/scientific honesty (pre-registration, reporting a failed recovery criterion, FDR/CIs, activation control, self-reported de-circularization audit); clear writing; open code.
- **Consensus technical risks:** (1) The effector arm may reflect generic NK-activation covariance rather than the specific axis (no specificity control); (2) the 37-gene list is driven by near-zero fold-changes and ~46% NK-side; (3) the ablation is largely tautological with no external predictive value, and the GNN adds no accuracy over LightGBM/XGBoost; (4) single-cell claims rest on one scRNA dataset with marker-threshold annotation.
- **Where emphasis differs:** R1 weights significance/novelty; R2 weights effect-size and self-referential analyses; R3 weights biological specificity and annotation robustness. All converge on the specificity-control and target-list-magnitude issues.
- **Broad-interest readout:** Useful to computational immuno-oncology (consistent with the BiB target); net new biological knowledge and breadth below a flagship bar; the honest framing underlines this.
- **Most important issues before a strong case:** (1) specificity control for the effector coupling; (2) effect-size justification/reframing for the target list + move caveat into abstract; (3) out-of-sample test replacing the ablation; (4) independent single-cell replication + reference-mapping sensitivity (H5); (5) a second card result or downgrade the "engine" claim.

## Risk / unsupported claims
- "Independent replication of the axis's functional endpoint" (§3.2) — weak/not established without a specificity control; partial r after the 16-marker control is 0.286.
- "Validates the `metabolic_crosstalk` edge's design premise" (§4.1) — overreach relative to r²<0.001.
- "37 putative tumor-intrinsic candidate targets" (abstract) — partly unsupported as headline (near-zero fold-changes, ~46% NK-side).
- Edge-ablation "confirmable structural effect beyond generic edges" (§3.7) — not assessable as biological validity (internal to the constructed edges).
- Multi-card reusability — not assessable (no in-text results).
- Figure-level claims (Figs 1–4) — not assessable (legends only).
- Fit to *Nature* — not asserted.
