# BiB Submission Bundle — Guide

**Target journal:** *Briefings in Bioinformatics* (Oxford University Press)
**Article type:** Problem Solving Protocol (original method + benchmarking)
**Prepared:** 2026-07-12
**Manuscript title:** *Transcriptional Reach of the
Serine–Sphingomyelin–Membrane-Topology Axis in NK-Cell Immune Evasion: A
Single-Cell Heterogeneous Graph Framework from Liver to Gastric Cancer*

This bundle is a **clean, self-contained snapshot** of everything needed to submit,
drawn from the canonical repository (`temp_clone`, GitHub `nblvguohao/GC-NKGraph-Atlas`).
Internal development notes (TDD logs, reviewer-simulation reports, memos, backups)
are deliberately **excluded** — they live in the repo's `submission_package/` for
your reference but are not part of a journal submission.

---

## 1. What's in this bundle

```
00_SUBMISSION_GUIDE.md          ← you are here
01_manuscript/
    main_manuscript.md          ← manuscript source (Markdown; human-readable)
    main.tex                    ← manuscript source (LaTeX; compile to PDF)
    cover_letter.md             ← cover letter (title synced, dated 2026-07-12)
02_figures/
    fig0_workflow               ← Workflow overview     (PDF + PNG) — embedded in main.tex via \includegraphics, Fig. label fig:workflow
    fig1_armA_positive_control  ← Arm A recovery        (PDF vector + PNG 300dpi)
    fig2_armB_extension         ← Arm B + external val.  (PDF + PNG)
    fig3_targets                ← candidate targets      (PDF + PNG)
    fig4_model_comparison       ← GNN vs baselines       (PDF + PNG)
03_supplementary/
    SUPPLEMENTARY_INDEX.md      ← index + Supplementary Methods prose (S.M.1–S.M.4)
    tables/*.tsv, *.md, *.json  ← 45 supplementary tables/summaries (2026-07-15 count; all 45 are referenced in SUPPLEMENTARY_INDEX.md — verified in sync)
04_reproducibility/
    environment.yml, requirements.txt, LICENSE (MIT)
```

## 2. ScholarOne upload map (BiB uses Oxford ScholarOne)

| Submission-system slot | File(s) from this bundle |
|---|---|
| Cover letter | `01_manuscript/cover_letter.md` (paste as text) |
| Main document | Compiled **`main.pdf`** (see §4) — title page, abstract, Key Points, body, references |
| Figures | `02_figures/fig0–4` — upload the **PDF** (vector) or 300-dpi PNG per portal rules; fig0 (workflow overview) is embedded in the compiled PDF via `\includegraphics` but still needs its own source file uploaded like fig1–4; each has a self-contained legend/alt text in the manuscript |
| Supplementary files | `03_supplementary/` — tables + `SUPPLEMENTARY_INDEX.md` (now includes the Supplementary Methods prose, S.M.1–S.M.4; see §6 for one file cited in-text that could not be located) |
| Data/code availability | Already stated in the manuscript; repo: https://github.com/nblvguohao/GC-NKGraph-Atlas |
| Suggested reviewers | Listed in the cover letter (Mustjoki, Theis, Li) |

## 3. Readiness status

**Done / verified (2026-07-12):**
- ✅ Manuscript complete: full IMRaD + Key Points + author biographies + running head + all declarations (Data Availability, Author Contributions/CRediT, Ethics, Competing Interests, Funding).
- ✅ All headline numbers cross-checked against result tables — consistent. Note the H3 single-cell number (r=0.32, pseudoreplication-corrected) is **not** a headline recovery result: it is reported for transparency only and is explicitly flagged in the manuscript as not surviving count-depth/latent-structure residualization or a random-module permutation baseline (collapses to r≈0.09; P=0.97). The effector-arm headline evidence is the three independent bulk results (TCGA-LIHC r=0.55; GSE62254 r=0.42; GSE84437 r=0.62). H5 Δ=−0.14; 37 candidates — consistent.
- ✅ Figures 1–4 present (vector PDF + 300-dpi PNG).
- ✅ Supplementary tables assembled (18) + index written; ablation table added.
- ✅ Cover letter now fully aligned with the manuscript body (2026-07-13 pass): title synced, and the effector-arm/metabolic-arm bullets rewritten so the cover letter no longer describes the single-cell H3 number as an independent replication or the single-cell H2 number as significant — both now match the manuscript's corrected, confound-controlled conclusions.
- ✅ Reproducibility verified: 120/120 unit tests pass; `python src/pipeline.py --synthetic` runs end-to-end (exit 0).
- ✅ `main.tex` and `main_manuscript.md` are content-synced (refs [1]–[50] in both; see §7 for the 2026-07-16 citation/retraction pass).

## 4. Remaining action items BEFORE you submit (author tasks)

1. **✅ PDF recompiled and verified (2026-07-16).** `main.pdf` was rebuilt from the current `main.tex` with a local MiKTeX pdflatex (two passes for cross-refs): **exit 0, 0 undefined references, no undefined citations** (all [1]–[48] resolve), **21 pages** (the modern/large OUP layout with figures runs longer than the earlier ~15-page estimate). Remaining warnings are cosmetic only (a few Overfull \hbox from a wide table + the template's per-page footer vbox). No further action needed unless you edit the source; if you do, re-run `pdflatex main.tex` twice before upload.
2. **⬜ Commit and push the repository.** `master` is even with `origin/master`; the doc-polish edits from the 2026-07-15/16 pass (Supplementary Methods, index sync, article-type line) are uncommitted in the working tree. Commit and push before submission so the code-availability link matches what reviewers read. See §7 for the branch-cleanup question before this push.
3. **⏸ ORCID iDs — deferred by author decision.** Placeholders `[0000-…]` remain in the manuscript by choice; fill the 7 authors' 16-digit iDs whenever ready (https://orcid.org). Not blocking the current pass.
4. **✅ CRediT roles & grant numbers verified (2026-07-16)** against `作者信息和基金.txt`. Authors, affiliations, corresponding authors (A. Zhou, L. Gu), and all grant numbers (NSFC 32472007, 62301006, 62301008; Anhui Prov. NSF 2308085MF217, 2308085QF202; Anhui Prov. Key Lab of Intelligent Agricultural Technology and Equipment) match across `main_manuscript.md`, `main.tex`, and the cover letter — the manuscript even corrects a stray-period typo in the source file's Anhui grant IDs. CRediT role assignments are internally consistent and not contradicted by the source file.
5. **✅ Article type confirmed: Problem Solving Protocol.** Declared in the cover letter (3×), recorded on the manuscript title page (`main_manuscript.md`) and as a `main.tex` comment; **select this same type in ScholarOne at submission.** Re-check against BiB's current *Instructions to Authors* at submission time in case categories change.
6. **✅ Supplementary Methods prose written** (S.M.1–S.M.4 in `SUPPLEMENTARY_INDEX.md`: NK-state labels/thresholds, SST-axis modules/scoring, scRNA QC/scVI/Leiden settings, model-comparison paired-test protocol). Skim once against BiB's supplementary format at submission time.
7. **⬜ Verify** BiB's live author guidelines (Oxford updates them): title-length, abstract structure, reference style (Vancouver numbered — already used), figure format/DPI.

## 5. Doc consistency note (resolved)

The mechanism-card count is now consistent at **4** across the YAML files, registry,
web playground, and README (`zheng_nk_sm_topology`, `adenosine_nk_suppression`,
`tgfb_nk_exclusion`, `nkg2d_micab_shedding`), with the README noting that only
`zheng` is executed end-to-end and the others are authored to demonstrate the format.
The manuscript body commits to no count, so the paper text is unaffected either way.

## 6. Code/data/figure/table alignment audit (2026-07-13)

**Extended pass (same day):** beyond filename existence, checked whether the
*behavior* described in Methods for each cited script matches what the script
actually does. This caught two further, more serious mismatches (Methods
describing something the code does not do / describes different numbers than
the code's actual defaults), on top of the `scrna_qc_summary.tsv` gap above:

- **RESOLVED:** §2.6 stated "When `torch_geometric` is available, a
  heterogeneous graph transformer (HGT) is used instead" for the gene
  embedding. A repo-wide search found no `torch_geometric`/`HGTConv` usage
  anywhere in `src/`; `src/models/gc_nkgraph_atlas.py` unconditionally uses
  the spectral (SVD-based) `GeneGraphEncoder` for every real result in this
  paper, including the §3.7 ablation. A separate `GATEncoder` class exists in
  the same file but is never instantiated anywhere in the codebase — dead
  code, not a runtime fallback. §2.6 and §4.4 reworded to state the spectral
  encoder was used throughout, with HGT named as unimplemented future work.
- **RESOLVED:** §2.6 stated the classifier optimizer was "Adam (lr=1e-3,
  weight_decay=1e-5)." The actual `NKStateClassifier` defaults (used
  unmodified by `pipeline.py`) are `learning_rate=1.7e-3, weight_decay=5.6e-6,
  dropout=0.6`, explicitly commented in code as selected by a Bayesian (TPE)
  search. The saved trial log (`gc_nkgraph_bayesian_trials.tsv`) has **100**
  completed trials, not the "50-trial" figure in the code's own comment — the
  code comment itself was stale; the manuscript now cites the real trial
  count from the data file. §2.6 reworded accordingly and the two supporting
  tables (`gc_nkgraph_bayesian_trials.tsv`, `gc_nkgraph_best_hyperparams.tsv`)
  added to the supplementary bundle. Note for the authors: the search's
  recorded `best_mcc` (0.728) does not exactly match Table 3's reported GNN
  MCC (0.706), while `best_auroc` (0.951) is very close to Table 3's AUROC
  (0.950) — plausibly the search's best-trial metric versus the final
  5-fold-mean metric reported in Table 3, but this was not chased further and
  is worth a second look if you want full numeric reconciliation.

Cross-checked every table/figure filename cited by name in `main_manuscript.md`
against what actually exists in `results/tables/`, `results/figures/`, and the
`03_supplementary/` bundle.

- **43 of the 74 files in `results/tables/`** are internal/intermediate artifacts
  (hyperparameter-search logs, `_v2`/backup duplicates, bootstrap-CI internals,
  bulk pipeline internals). None of these 31 excluded files are cited by name
  anywhere in the manuscript text — the authors already correctly kept them out
  of the reader-facing supplementary bundle. No action needed.
- **RESOLVED (2026-07-13):** `results/tables/scrna_qc_summary.tsv` was cited in
  Methods §2.4 ("QC thresholds are logged and reported in
  `results/tables/scrna_qc_summary.tsv`") but did not exist anywhere in this
  repository's git history (checked across all branches, including the two
  since deleted from origin) — it was evidently never brought back from
  whatever remote/server environment ran the real QC step, and
  `results/tables/` is gitignored by policy in any case. Since it could not be
  retrieved, Methods §2.4 has been reworded to stop asserting the file exists:
  it now states the thresholds are fixed in `src/scrna_analysis/qc_filter.py`
  and that per-sample before/after counts were not preserved for this
  submission. If the original file later turns up, it can be added to
  `03_supplementary/tables/` and the sentence restored.
- **RESOLVED (2026-07-13):** the internal version-changelog block referencing
  `manuscript/notes/SUBMISSION_READINESS.md` has been removed from
  `main_manuscript.md` (it was also interrupting the numbered reference list).
  The Discussion §4.1 reference to the private
  `plan_NatImmunol_wetlab/WETLAB_PROGRAM.md` path has been reworded to
  describe the companion wet-lab program generically, since that file is
  intentionally staying out of the public repo.
- **Figures:** `01_manuscript/figures/` (LaTeX source) and `02_figures/`
  (upload bundle) now both contain `fig0_workflow` through `fig4` consistently
  (fig0's PDF was missing from `02_figures/` before this pass — added). Note
  `results/figures/` additionally contains three stale, superseded figures not
  used anywhere in the current manuscript (`fig9_sst_axis_scores.pdf`,
  `fig9_sst_axis_scrna.pdf`, `fig10_liver_positive_control.pdf`, dated 2026-07-07
  and 2026-07-12) — harmless (not part of the submission bundle) but candidates
  for repo cleanup if you want the pipeline output directory to reflect only
  current figures.
- **RESOLVED (2026-07-15):** `SUPPLEMENTARY_INDEX.md`'s table listing (previously
  written for an earlier ~18-table version) has been re-synced to the current
  45-file set — all 45 files in `03_supplementary/tables/` are now referenced in
  the index, with no dangling entries for files that do not exist (verified by
  filename diff). The Supplementary Methods prose (S.M.1–S.M.4) was also added.

## 7. 2026-07-16 (later pass): related-work citations added, one ablation claim retracted after a bug fix

- **Added two citations** (refs [49] TREE — Su et al., *Nat Biomed Eng* 2025 — and
  [50] GRAFT — Cho & Cho, *Brief Bioinform* 2026, the same journal as this
  submission) to §1.4 Related Work and §4.4 Future Directions, explicitly
  positioning this paper's ~100-gene mechanism-anchored graph against that
  line of genome-scale, multi-network driver-gene-discovery work rather than
  leaving the comparison unaddressed. Reference count is now **[1]–[50]**
  (`main.tex` bibitem count and `main_manuscript.md` numbered-reference count
  both verified at 50; highest in-text citation number in both files is 50).
- **Found and fixed a real bug**, independent of the citation work: in both
  `src/a100_recompute/run_t17_edge_external_value.py` (`build_adj_matrix`) and
  `src/a100_recompute/run_local_t14_t7.py` (`build_adj`), when two edge types
  connected the same node pair, whichever edge type's row came later in
  `edges.tsv` silently overwrote the earlier one's weight (last-write-wins)
  instead of the two being combined. Fixed to take `max(weight)` across
  colliding edge types in both scripts.
- **§3.7's H1/H2/modularity ablation table numbers changed** after the fix
  (e.g. FULL H1 0.140→0.128) but the qualitative conclusion is unchanged and
  was re-verified: removing `metabolic_crosstalk` still abolishes the H1
  embedding-coupling (0.128→−0.002). `results/tables/ablation_results.tsv`
  and the supplementary copy were both regenerated and updated in-text.
- **§3.7's cross-cohort transfer claim was retracted, not just corrected.**
  The original single-seed result (FULL MCC=0.416 vs −MC MCC=0.509, read as
  "the edge measurably degrades cross-cohort transfer") turned out to rest on
  both the bug above and a single training seed. Refitting with the bug fixed
  and averaged over 10 seeds gave a third, different answer each time the
  seed count changed (3 seeds: no detectable difference; 10 seeds: FULL
  significantly *better*, p=0.01–0.02) — evidence the test itself is too
  seed-unstable to support a directional claim either way, not evidence for a
  reversed conclusion. Rather than swap in a new number, the cross-cohort
  transfer sub-claim was removed entirely from Key Points, §3.7, §4.1 (item
  2), §4.3 (Limitations item 8), and §5 Conclusion, in both `main.tex` and
  `main_manuscript.md`. The paper's ablation evidence now rests solely on the
  H1 embedding-coupling result, which does not depend on this test.
  `t17_edge_external_value.tsv` was removed from `03_supplementary/tables/`
  and its `SUPPLEMENTARY_INDEX.md` row deleted (46 files, 46 index entries,
  verified in sync by filename diff — zero dangling or missing entries).
- **`main.pdf` recompiled** (two `pdflatex` passes, MiKTeX): exit 0, 0
  undefined references/citations, 22 pages. Figures 1–4 and fig0 (workflow)
  were checked and do **not** reference any of the removed content — no
  figure regeneration was needed for the retraction itself.
- **Separately, found and fixed a real Figure 4 / Table 3 mismatch** while
  visually reviewing the figures: `src/figures/make_figures.py`'s
  `figure4()` read `results/tables/model_comparison.tsv`, which only has 7
  of the 9 methods reported in Table 3 — it was missing the two
  domain-knowledge baselines (SST-module signature, NK-marker signature).
  Figure 4 and its caption ("significantly above the NK-marker signature...")
  therefore made a claim the figure itself didn't visually support.
  `results/tables/domain_baselines_per_fold.tsv` is a verified complete
  superset (identical values for all 7 shared methods, confirmed by exact
  comparison) with the 2 missing methods included; `figure4()` now reads
  from that file instead. Regenerated and redeployed to
  `results/figures/`, `02_figures/`, and `01_manuscript/figures/`; `main.pdf`
  recompiled again (still 22 pages, 0 undefined refs). All 9 methods now
  visible in Figure 4, matching Table 3 exactly.
- Net effect on the paper's central claims: **none** — the three-layer
  scoping map, the two-arm design, the 37 candidate targets, and the
  dual-mechanism-card generalization result are all unchanged. The only
  removed material was a secondary, ultimately unsupportable side-claim about
  the graph's downstream predictive value, which was always in tension with
  the paper's own "probe, not predictor" framing; removing it made that
  framing more internally consistent, not less.

## 8. 2026-07-17: merged a concurrent branch adding real multi-modal evidence (S3.7/S3.8, Fig. S1/S2)

A separate concurrent agent session (`codex/multiview-strengthening`, forked
from this repo's `705bf48`) independently built two substantive, real-data
additions and asked to have them merged into `master`:

1. **Label-masked multi-view external audit** (extends S2.6/S3.7, adds Fig.
   S1). Properly redoes this session's own earlier exploratory multi-view-
   fusion idea (informally "T19"): masks every label-defining gene out of
   both expression and graph projections, uses 10 fixed seeds with STAD
   80:20 train/early-stop and LIHC held out entirely from tuning. Finding:
   learned multi-view fusion is significantly *worse* than no graph once
   label overlap is controlled -- this session's own earlier, unmasked T19
   result was optimistic due to exactly that circularity, and was correctly
   never added to the manuscript at the time.
2. **Real-data comparative recoverability atlas** (new S3.8, Fig. S2). Tests
   all four registered mechanism cards across the four bulk cohorts plus
   three additional real orthogonal-modality sources (GSE122401 protein,
   MTBLS3303 metabolomics, GSE251950 Visium spatial). Verdict:
   `comparative_atlas_only` (pre-registered gate not met) -- protein/
   metabolomics correctly reported `not_measured` rather than faked; spatial
   evidence scoped to an explicitly exploratory 4-of-9/10-section subset
   after the full archive failed an integrity check and was excluded.

**Verification before merging:** ran the new test suite (43 passed, 3 skipped
for the still-pending metabolomics download) plus the full existing suite
(163 passed total, no regressions); cross-checked one headline number
(learned-fusion-vs-no-graph AUROC/AUPRC deltas) against its underlying TSV
and got an exact match; confirmed the no-synthetic-data guarantees are
enforced in code (`src/common/real_data.py` rejects `synthetic`/`mock`/`demo`
paths), not just documented.

**Merge conflicts resolved by hand:** Limitations item 8 (git's auto-merge
picked one side and silently dropped this session's T17
cross-cohort-transfer-instability clause -- manually recombined both
cautions); `main.pdf` (binary conflict, resolved by recompiling from the
merged `.tex`, 23 pages, 0 undefined refs); `SUPPLEMENTARY_INDEX.md` (merged
cleanly, then audited and fixed 7 table entries the source branch had
created but never indexed). The "Additional mechanisms" Future Directions
paragraph auto-merged correctly in favor of this session's newer 8-card text
(the other branch's version was stale, predating the mechanism-card commit).

Supplementary table count is now **63** (was 46); figure count is now **7**
(fig0-4 plus S1, S2). `data/external/` (1.1GB of real downloaded
GSE122401/GSE251950 data backing the new work) added to `.gitignore`, same
policy as `data/26Q1/`.

## 9. Open question before pushing to GitHub

The local repo has five branches, four of which exist on `origin` alongside
`master`: `nk-pre-submission`, `strengthen-paper`,
`system-optimization-2026-07-12`, and `gh-pages-deploy`/`gh-pages`. Before any
push or branch cleanup, confirm with the repo owner exactly which of these are
superseded and safe to remove from the public remote (with a local backup kept
first) versus which must stay (e.g. `gh-pages` likely serves the GitHub Pages
web playground referenced in the manuscript and should probably not be
deleted). See conversation for the pending clarification.
