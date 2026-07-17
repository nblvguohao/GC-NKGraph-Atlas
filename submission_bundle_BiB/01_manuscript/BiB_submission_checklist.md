# Briefings in Bioinformatics — Submission Checklist (current)

> ✅ done · ⚠️ verify · ⬜ author action. Re-check every item against BiB's live
> *Instructions to Authors* at submission — Oxford updates them.
> Last updated: **2026-07-16** (supersedes the 2026-07-12 checklist; see
> `00_SUBMISSION_GUIDE.md` §7 for what changed in this pass — two new
> citations added, one unsupportable ablation sub-claim retracted after a
> code bug was found and fixed).

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
- [x] ✅ Figures 1–4 + S1 (label-masked multi-view audit) + S2 (real-data recoverability atlas) (vector PDF + 300-dpi PNG); numbers read from result tables. S1/S2 added 2026-07-17 via merge from `codex/multiview-strengthening` (see guide §8).
- [x] ✅ **Ablation (§3.7)** — done. Edge-type ablation table included in supplementary (`ablation_results.tsv`). A collision bug in the adjacency-construction code was found and fixed on 2026-07-16 (see guide §7); the table numbers were regenerated and the qualitative conclusion re-verified as unchanged. A separate cross-cohort transfer sub-claim that had rested on the same bug proved seed-unstable after the fix and was removed from the manuscript rather than replaced with a new number — the ablation section now presents only the (robust) embedding-coupling result.
- [x] ✅ All headline stats re-verified against result tables (2026-07-12; ablation table re-verified again 2026-07-16 after the bug fix above).

## D. Front matter & declarations
- [x] ✅ Authors + affiliations; corresponding authors (L. Gu, A. Zhou).
- [x] ✅ Author Contributions (CRediT), Funding, Competing Interests, Data Availability, Ethics — all present.
- [x] ✅ Code availability: repo + reproduction commands + synthetic mode + MIT LICENSE.
- [x] ✅ ORCID iDs — by author decision, not included in the manuscript text (placeholders removed 2026-07-17); each author will enter their iD directly in Oxford ScholarOne at submission instead.
- [x] ✅ Grant numbers / CRediT roles verified against `作者信息和基金.txt` (2026-07-16) — see guide §4 item 4.

## E. References & formatting
- [x] ✅ References [1]–[50], Vancouver numbered; refs [49]–[50] (TREE, GRAFT — added 2026-07-16 to position the graph design against recent multi-network driver-gene-discovery work) integrated in both `.tex` and `.md`; bibitem count and highest in-text citation number both verified at 50 in each file.
- [x] ✅ Figures vector + ≥300-dpi, self-contained legends; verified 2026-07-16 that none of fig0–4 depend on the retracted cross-cohort transfer claim. Separately, Figure 4 was regenerated (2026-07-16) to fix a real gap: it previously plotted only 7 of Table 3's 9 methods (missing the SST-module-signature and NK-marker-signature baselines), so the caption's "significantly above the NK-marker signature" claim wasn't visually supported by the figure. All 9 methods now shown, matching Table 3.
- [x] ✅ Supplementary tables (63, per `03_supplementary/tables/`) + `SUPPLEMENTARY_INDEX.md`, verified in sync by filename diff (2026-07-17, re-verified after the multiview/recoverability merge): 63 files present, 63 referenced, zero dangling or missing entries.
- [x] ✅ Supplementary Methods prose present (S.M.1–S.M.4 in `SUPPLEMENTARY_INDEX.md`).

## F. Reproducibility (BiB values this)
- [x] ✅ 120/120 unit tests pass.
- [x] ✅ `python src/pipeline.py --synthetic` runs end-to-end (exit 0).
- [x] ✅ `environment.yml` + `requirements.txt` + Dockerfile in repo.

## G. Submission-system items
- [x] ✅ Cover letter finalized (dated 2026-07-12; 3 suggested reviewers w/ emails).
- [x] ✅ Keywords ≤ 6.
- [x] ✅ **`main.pdf` recompiled** (2026-07-16, MiKTeX, two `pdflatex` passes): exit 0, 0 undefined references/citations, 22 pages. Re-run this after any further edit to `main.tex`.
- [ ] ⬜ **Push repo** — commit and push before submission so the public code-availability link matches the submitted version (see guide §8 for the branch-cleanup question first).
- [ ] ⬜ Register corresponding authors in Oxford ScholarOne; confirm no dual submission; all authors approve.
- [ ] ⬜ Graphical abstract / highlights only if the editorial office requests.

---

### The only things between you and "submit"
1. ~~Push the repo to GitHub~~ — done 2026-07-17; `master` is even with `origin/master` (commit `523130a`). The three stale branches flagged in guide §8 were resolved: `nk-pre-submission` and two fully-merged local branches deleted, `origin/merged-paper` (fully superseded) deleted from the remote; `origin/gh-pages` left untouched (serves the web playground).
2. ~~ORCID~~ — resolved by author decision: removed from the manuscript entirely (2026-07-17); each author enters their iD directly in Oxford ScholarOne at submission instead.
3. Register in Oxford ScholarOne; confirm no dual submission; all authors approve (§G).
4. One last skim of BiB's live *Instructions to Authors* immediately before submission — a same-day web check (2026-07-17) could not confirm current word-limit/article-type figures reliably (see chat record); cross-check directly on academic.oup.com/bib/pages/msprep_submission at submission time, and note the manuscript body (~13,100 words, Intro-Conclusion) may be well over a "2,000-5,000 word" Problem Solving Protocol figure if that number is still current -- weigh against the fact that recently published BiB Problem Solving Protocols (e.g. GRAFT, bbaf706, 2026) are of comparable length to this manuscript.

Everything else (analysis, figures, tables, cover letter, declarations,
reproducibility, PDF compilation, reference list, supplementary-table sync,
mechanism-card registry) is complete and verified as of 2026-07-17.
