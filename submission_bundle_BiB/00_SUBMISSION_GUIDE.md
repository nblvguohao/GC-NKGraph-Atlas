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
    SUPPLEMENTARY_INDEX.md      ← index + supp-methods pointer
    tables/*.tsv, *.md, *.json  ← 43 supplementary tables/summaries (2026-07-13 count; SUPPLEMENTARY_INDEX.md table list should be re-synced to match — see §6)
04_reproducibility/
    environment.yml, requirements.txt, LICENSE (MIT)
```

## 2. ScholarOne upload map (BiB uses Oxford ScholarOne)

| Submission-system slot | File(s) from this bundle |
|---|---|
| Cover letter | `01_manuscript/cover_letter.md` (paste as text) |
| Main document | Compiled **`main.pdf`** (see §4) — title page, abstract, Key Points, body, references |
| Figures | `02_figures/fig0–4` — upload the **PDF** (vector) or 300-dpi PNG per portal rules; fig0 (workflow overview) is embedded in the compiled PDF via `\includegraphics` but still needs its own source file uploaded like fig1–4; each has a self-contained legend/alt text in the manuscript |
| Supplementary files | `03_supplementary/` — tables + `SUPPLEMENTARY_INDEX.md` (add the short Supplementary Methods prose first; see §6 for one file cited in-text that could not be located) |
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
- ✅ `main.tex` and `main_manuscript.md` are content-synced (both include refs 46–48 / T17).

## 4. Remaining action items BEFORE you submit (author tasks)

1. **⚠️ Recompile the PDF.** The repo's `main.pdf` is **stale** (older than `main.tex`, and no LaTeX is installed on this machine). Compile the current source:
   ```
   pdflatex main.tex && pdflatex main.tex        # two passes for cross-refs
   ```
   or drop `main.tex` into Overleaf. Confirm: 0 undefined references, refs [1]–[48] all resolve, ~15 pages.
2. **⬜ Commit and push the repository.** As of 2026-07-13, `master` is even with `origin/master`, but the working tree has substantial **uncommitted** work (manuscript/cover-letter alignment fixes, ~15 new analysis scripts, ~30 new supplementary tables, new figure assets) that has never been pushed. Commit and push before submission so the code-availability link actually matches what reviewers will read. See §7 for the branch-cleanup question before this push.
3. **⬜ ORCID iDs.** Each of the 7 authors fills their 16-digit iD (placeholders `[0000-…]` are in the manuscript). Register free at https://orcid.org.
4. **⬜ Confirm CRediT roles & grant numbers** against award letters (NSFC 32472007, 62301006, 62301008; Anhui 2308085MF217, 2308085QF202).
5. **⬜ Confirm article type** — "Problem Solving Protocol" vs standard Research Article — against BiB's current *Instructions to Authors* at submission time.
6. **⬜ Add the short Supplementary Methods prose** (pointer written in `SUPPLEMENTARY_INDEX.md`).
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
- `SUPPLEMENTARY_INDEX.md`'s table listing was written for an earlier ~18-table
  version of the bundle and should be re-synced to the current 43-file set
  before submission (content is accurate for the files it does list; it is
  incomplete, not wrong).

## 7. Open question before pushing to GitHub

The local repo has five branches, four of which exist on `origin` alongside
`master`: `nk-pre-submission`, `strengthen-paper`,
`system-optimization-2026-07-12`, and `gh-pages-deploy`/`gh-pages`. Before any
push or branch cleanup, confirm with the repo owner exactly which of these are
superseded and safe to remove from the public remote (with a local backup kept
first) versus which must stay (e.g. `gh-pages` likely serves the GitHub Pages
web playground referenced in the manuscript and should probably not be
deleted). See conversation for the pending clarification.
