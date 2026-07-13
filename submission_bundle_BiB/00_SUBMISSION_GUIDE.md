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

Cross-checked every table/figure filename cited by name in `main_manuscript.md`
against what actually exists in `results/tables/`, `results/figures/`, and the
`03_supplementary/` bundle.

- **43 of the 74 files in `results/tables/`** are internal/intermediate artifacts
  (hyperparameter-search logs, `_v2`/backup duplicates, bootstrap-CI internals,
  bulk pipeline internals). None of these 31 excluded files are cited by name
  anywhere in the manuscript text — the authors already correctly kept them out
  of the reader-facing supplementary bundle. No action needed.
- **Every supplementary file cited by name in the manuscript is present in the
  bundle, with one exception:** `results/tables/scrna_qc_summary.tsv`, cited in
  Methods §2.4 ("QC thresholds are logged and reported in
  `results/tables/scrna_qc_summary.tsv`"), **does not exist anywhere in this
  repository.** The script that produces it
  (`src/scrna_analysis/qc_filter.py`) is explicitly marked in its own docstring
  as "untested on the author's workstation (scanpy not installed here)" and
  requires `data/processed/scrna/gc_integrated.h5ad` (the full 166,829-cell
  pre-QC integrated object) as input — that file is also not present locally;
  only the already-QC'd NK subset (`gc_nk_subset_remote.h5ad`) is. This table
  was evidently produced on a remote/server run and never brought back into
  this repo. **This is a genuine, unresolved gap between what the manuscript
  claims exists and what is currently reproducible/verifiable — flagged for the
  authors to either (a) retrieve the original file from wherever the remote QC
  run was performed, (b) re-run QC filtering remotely and commit the output, or
  (c) soften the Methods §2.4 sentence if the table cannot be recovered.**
- **Two other filenames cited in-text are not in `03_supplementary/` by
  design, not by omission:** `manuscript/notes/SUBMISSION_READINESS.md`
  (referenced only in the manuscript's own internal version-changelog block,
  which itself should be removed before formal submission — see the
  nature-reviewer assessment, Reviewer 3) and
  `plan_NatImmunol_wetlab/WETLAB_PROGRAM.md` (referenced in Discussion §4.1 as
  a companion future-work document, not a supplementary data file). Neither
  needs to be added to the supplementary bundle.
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
