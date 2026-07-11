# BiB Submission Bundle — Guide

**Target journal:** *Briefings in Bioinformatics* (Oxford University Press)
**Article type:** Problem Solving Protocol (original method + benchmarking)
**Prepared:** 2026-07-12
**Manuscript title:** *Mapping the Transcriptional Reach of the
Serine–Sphingomyelin–Membrane-Topology Axis of NK-Cell Immune Evasion: A
Single-Cell-Informed Heterogeneous Graph Framework from Liver to Gastric Cancer*

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
    fig1_armA_positive_control  ← Arm A recovery        (PDF vector + PNG 300dpi)
    fig2_armB_extension         ← Arm B + external val.  (PDF + PNG)
    fig3_targets                ← candidate targets      (PDF + PNG)
    fig4_model_comparison       ← GNN vs baselines       (PDF + PNG)
03_supplementary/
    SUPPLEMENTARY_INDEX.md      ← index + supp-methods pointer
    tables/*.tsv, *.md, *.json  ← 18 supplementary tables
04_reproducibility/
    environment.yml, requirements.txt, LICENSE (MIT)
```

## 2. ScholarOne upload map (BiB uses Oxford ScholarOne)

| Submission-system slot | File(s) from this bundle |
|---|---|
| Cover letter | `01_manuscript/cover_letter.md` (paste as text) |
| Main document | Compiled **`main.pdf`** (see §4) — title page, abstract, Key Points, body, references |
| Figures | `02_figures/fig1–4` — upload the **PDF** (vector) or 300-dpi PNG per portal rules; each has a self-contained legend in the manuscript |
| Supplementary files | `03_supplementary/` — tables + `SUPPLEMENTARY_INDEX.md` (add the short Supplementary Methods prose first) |
| Data/code availability | Already stated in the manuscript; repo: https://github.com/nblvguohao/GC-NKGraph-Atlas |
| Suggested reviewers | Listed in the cover letter (Mustjoki, Theis, Li) |

## 3. Readiness status

**Done / verified (2026-07-12):**
- ✅ Manuscript complete: full IMRaD + Key Points + author biographies + running head + all declarations (Data Availability, Author Contributions/CRediT, Ethics, Competing Interests, Funding).
- ✅ All headline numbers cross-checked against result tables (H3 r=0.55/0.32; H5 Δ=−0.14; external 0.42/0.62; 37 candidates) — consistent.
- ✅ Figures 1–4 present (vector PDF + 300-dpi PNG).
- ✅ Supplementary tables assembled (18) + index written; ablation table added.
- ✅ Cover letter title synced to the current manuscript title (was previously the old "Reconstructing…" title — now fixed).
- ✅ Reproducibility verified: 120/120 unit tests pass; `python src/pipeline.py --synthetic` runs end-to-end (exit 0).
- ✅ `main.tex` and `main_manuscript.md` are content-synced (both include refs 46–48 / T17).

## 4. Remaining action items BEFORE you submit (author tasks)

1. **⚠️ Recompile the PDF.** The repo's `main.pdf` is **stale** (older than `main.tex`, and no LaTeX is installed on this machine). Compile the current source:
   ```
   pdflatex main.tex && pdflatex main.tex        # two passes for cross-refs
   ```
   or drop `main.tex` into Overleaf. Confirm: 0 undefined references, refs [1]–[48] all resolve, ~15 pages.
2. **⬜ Push the repository.** The local branch is **2 commits ahead of GitHub** (the T17 analysis + the latest manuscript sync are not yet on the public repo). Push `master` so the code-availability link matches the submitted manuscript, and confirm the repo (and GitHub Pages playground) is public. *(Outward action — left for you to run.)*
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
