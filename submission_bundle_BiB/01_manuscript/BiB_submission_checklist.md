# Briefings in Bioinformatics — Submission Checklist (current)

> ✅ done · ⚠️ verify · ⬜ author action. Re-check every item against BiB's live
> *Instructions to Authors* at submission — Oxford updates them.
> Last updated: **2026-07-12** (supersedes the 2026-07-07 checklist).

## A. Article type & scope
- [x] ✅ Type: **Problem Solving Protocol** (original method + benchmarking). ⚠️ confirm vs. Research Article at submission.
- [x] ✅ Scope fit: computational method of broad interest; benchmarked; transferable abstraction.

## B. Manuscript structure
- [x] ✅ Title (present; cover letter title now **synced** to manuscript).
- [x] ✅ **Key Points** (3–5 bullets) — present.
- [x] ✅ Structured Abstract (Motivation/Results/Availability/Contact).
- [x] ✅ Intro / Methods / Results / Discussion / Conclusion + Related work (§1.4).
- [x] ✅ Running head (GC-NKGraph-Atlas); author biographies (~30 words × 7).

## C. Results completeness
- [x] ✅ Table 1 datasets; Table 2 positive-control recovery; Table 3 GNN + baselines (paired tests); Table 4 candidates (n=37); Table 5 external validation (0.42/0.62).
- [x] ✅ Figures 1–4 (vector PDF + 300-dpi PNG); numbers read from result tables.
- [x] ✅ **Ablation (§3.7)** — done. Edge-type ablation table included in supplementary (`ablation_edge_types.tsv`, `ablation_results.tsv`). *(This was the last "partial" item on the old checklist — now closed.)*
- [x] ✅ All headline stats re-verified against result tables (2026-07-12).

## D. Front matter & declarations
- [x] ✅ Authors + affiliations; corresponding authors (L. Gu, A. Zhou).
- [x] ✅ Author Contributions (CRediT), Funding, Competing Interests, Data Availability, Ethics — all present.
- [x] ✅ Code availability: repo + reproduction commands + synthetic mode + MIT LICENSE.
- [ ] ⬜ ORCID iDs — 7 authors each fill their 16-digit iD (placeholders in text).
- [ ] ⚠️ Verify grant numbers / CRediT roles against source documents.

## E. References & formatting
- [x] ✅ References [1]–[48], Vancouver numbered; refs 46–48 (T17 literature) integrated in both `.tex` and `.md`.
- [x] ✅ Figures vector + ≥300-dpi, self-contained legends.
- [x] ✅ Supplementary tables (18) + `SUPPLEMENTARY_INDEX.md`.
- [ ] ⬜ Add short **Supplementary Methods** prose (pointer in the index).

## F. Reproducibility (BiB values this)
- [x] ✅ 120/120 unit tests pass.
- [x] ✅ `python src/pipeline.py --synthetic` runs end-to-end (exit 0).
- [x] ✅ `environment.yml` + `requirements.txt` + Dockerfile in repo.

## G. Submission-system items
- [x] ✅ Cover letter finalized (dated 2026-07-12; 3 suggested reviewers w/ emails).
- [x] ✅ Keywords ≤ 6.
- [ ] ⚠️ **Recompile `main.pdf`** from current `main.tex` (repo PDF is stale; no LaTeX here). Two `pdflatex` passes; expect ~15 pp, 0 undefined refs.
- [ ] ⬜ **Push repo** — local is 2 commits ahead of GitHub; make sure public repo matches submitted version.
- [ ] ⬜ Register corresponding authors in Oxford ScholarOne; confirm no dual submission; all authors approve.
- [ ] ⬜ Graphical abstract / highlights only if the editorial office requests.

---

### The only things between you and "submit"
1. Recompile PDF (§G).
2. Push the 2 pending commits to GitHub (§G).
3. Fill ORCID iDs + confirm grants/CRediT (§D).
4. Add Supplementary Methods prose (§E).

Everything else (analysis, figures, tables, cover letter, declarations,
reproducibility) is complete and verified.
