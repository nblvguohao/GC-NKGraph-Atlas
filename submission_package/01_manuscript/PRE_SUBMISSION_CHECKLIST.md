# GC-NKGraph-Atlas — Pre-Submission Verification Checklist

> **Target:** *Briefings in Bioinformatics* (Oxford University Press)
> **Article type:** Problem Solving Protocol
> **Date:** 2026-07-09
> **Status:** All data/analysis blockers cleared. Remaining: author actions + formatting.

This document is the single source of truth for every action needed before
clicking "Submit" in ScholarOne. Check off items as each author completes them.

---

## 🔴 BLOCKING — Must complete before submission

### A. Author ORCID iDs (each author, individually)

Each author must register at https://orcid.org and fill their 16-digit ORCID iD
below. The placeholder `[0000-0000-0000-0000]` in `main_manuscript.md` must be
replaced.

| # | Author | ORCID iD | Verified |
|---|--------|----------|----------|
| 1 | Guohao Lyu | `____-____-____-____` | ⬜ |
| 2 | Yingchun Xia | `____-____-____-____` | ⬜ |
| 3 | Huichao Liu | `____-____-____-____` | ⬜ |
| 4 | Xiaolei Zhu | `____-____-____-____` | ⬜ |
| 5 | Shuai Yang | `____-____-____-____` | ⬜ |
| 6 | Ailian Zhou (corresponding) | `____-____-____-____` | ⬜ |
| 7 | Lichuan Gu (corresponding) | `____-____-____-____` | ⬜ |

### B. Public code repository

- [ ] Push the complete repository to https://github.com/nblvguohao/GC-NKGraph-Atlas
- [ ] Verify the repo is **public** (not private)
- [ ] Add repository description: "Reconstructing the serine–sphingomyelin–membrane-topology axis of NK-cell immune evasion from tumor transcriptomes — a single-cell-informed heterogeneous graph framework"
- [ ] Add repository topics: `bioinformatics`, `single-cell`, `graph-neural-network`, `cancer-immunology`, `nk-cells`, `transcriptomics`
- [ ] On a **clean machine or fresh clone**, run: `pytest -q` (all tests pass)
- [ ] On a **clean machine or fresh clone**, run: `python src/pipeline.py --synthetic` (end-to-end succeeds)
- [ ] Verify LICENSE file is present (MIT) and visible on GitHub

### C. Manuscript formatting (BiB requirements)

- [ ] Replace all `[0000-0000-0000-0000]` placeholders with real ORCID iDs
- [ ] Delete the `*(ORCID iDs: fill in each author's ORCID...)*` note below author list
- [ ] Delete the `*(Biographies ~30 words each...)*` note below author biographies
- [ ] Delete the `*(CRediT taxonomy; adjust...)*` note in Author Contributions
- [ ] Author biographies: each author confirms their 30-word biography is accurate
- [ ] CRediT roles: each author confirms their contribution roles
- [ ] Check that all references [1]–[24] resolve correctly on PubMed/CrossRef
- [ ] Verify grant numbers against award letters (NSFC 32472007, 62301006, 62301008; Anhui 2308085MF217, 2308085QF202)
- [ ] Verify the running head "GC-NKGraph-Atlas" is present

### D. Figures

- [ ] Figure 1–4 regenerated from current result tables (`python src/figures/make_figures.py`)
- [ ] Figure 5 (mechanism-card concept) generated
- [ ] All figures: PDF (vector) + PNG (≥300 dpi), colour-blind-safe palette
- [ ] Each figure has a self-contained legend
- [ ] Figure files copied to `manuscript/figures/` directory
- [ ] All figures referenced in text in sequential order

### E. ScholarOne submission system

- [ ] Corresponding authors registered at https://mc.manuscriptcentral.com/bib
- [ ] ORCID iDs linked in ScholarOne profiles
- [ ] Cover letter uploaded
- [ ] Manuscript file uploaded (as .docx or LaTeX → PDF per BiB instructions)
- [ ] Supplementary tables uploaded as separate files or combined PDF
- [ ] Figures uploaded individually or embedded in manuscript per BiB format
- [ ] Suggested reviewers entered (Satu Mustjoki, Fabian J. Theis, Xia Li)
- [ ] Opposed reviewers (if any) entered
- [ ] All authors receive and approve the submission confirmation

---

## 🟡 IMPORTANT — Strongly recommended before submission

### F. Content verification

- [ ] All numeric claims in manuscript cross-checked against `results/tables/*.tsv`
- [ ] Table 2 (H1–H5 outcomes): r, p, and verdict match `sst_axis_positive_control_recovery.tsv`
- [ ] Table 3 (model comparison): all 7 methods match `model_comparison.tsv`
- [ ] Table 4 (top targets): ranks and scores match `tumor_intrinsic_candidates.tsv`
- [ ] Table 5 (external validation): r and p match `external_validation_results.tsv`
- [ ] No placeholder text remains: `grep -nE '\[PLACEHOLDER|\[TBD|To be (expanded|written|designed)' main_manuscript.md` returns empty
- [ ] No unmarked ORCID placeholders: `grep -c '0000-0000-0000-0000' main_manuscript.md` returns 0

### G. Reproducibility

- [ ] `README.md` contains clear reproduction instructions
- [ ] `environment.yml` lists all dependencies with versions
- [ ] `requirements.txt` minimal install works (`pip install -r requirements.txt`)
- [ ] Synthetic data mode (`--synthetic`) exercises full pipeline
- [ ] `.gitignore` excludes data files and results (but tracks `submission_package/`)
- [ ] All Python scripts have a `__main__` guard or are importable without side effects

### H. Honesty guardrails

- [ ] No claim of "predicting physical membrane topology" from transcriptomes
- [ ] All axis claims bounded by "transcriptional program permissive-of / associated-with"
- [ ] Positive-control outcome reported per-hypothesis (no cherry-picking)
- [ ] NK-side readout panel explicitly labeled as "not targets"
- [ ] GNN described as "on par with" / "comparable to" top baselines, not "outperforming"
- [ ] No survival association presented as mechanism
- [ ] In-silico SM-restoration stratified as "hypothesis for experimental testing"

---

## 🟢 NICE-TO-HAVE — If time permits

- [ ] Run `qc_filter.py` on scRNA data and update QC table in Methods §2.4
- [ ] Verify QC does not change H3/H2 conclusion direction (regression test)
- [ ] Generate a graphical abstract for the journal submission system
- [ ] Pre-registration calibration note logged in `configs/sst_axis_config.yaml`
- [ ] Add Supplementary Methods note explaining the supplementary tables
- [ ] Ask one colleague unfamiliar with the project to read the manuscript and flag unclear passages

---

## ⚡ Quick self-check commands

```bash
# Verify no placeholders remain
grep -nE '\[PLACEHOLDER|\[TBD|To be (expanded|written|designed)' manuscript/main_manuscript.md

# Verify ORCIDs are filled (should return 0)
grep -c '0000-0000-0000-0000' manuscript/main_manuscript.md

# Regenerate all figures
python src/figures/make_figures.py --dpi 300

# Run synthetic pipeline end-to-end
python src/pipeline.py --synthetic

# Run all tests
pytest -q

# Count words (rough)
wc -w manuscript/main_manuscript.md
```

---

## Submission timeline

| Step | Owner | Est. time |
|------|-------|-----------|
| ORCID registration | All authors | 10 min each |
| Push to GitHub | Guohao Lyu | 30 min |
| Regenerate figures | Guohao Lyu | 5 min |
| Format check & proofread | Guohao Lyu / Lichuan Gu | 2 hr |
| Grant number verification | Lichuan Gu / Ailian Zhou | 15 min |
| ScholarOne registration | Lichuan Gu / Ailian Zhou | 30 min |
| Final author approval | All authors | 1 hr |
| **Click Submit** | Lichuan Gu or Ailian Zhou | — |

**Estimated time to submission: 1–2 days of coordinated effort.**

---

*This checklist was generated 2026-07-09. Check https://academic.oup.com/bib/pages/author-guidelines for any BiB policy changes before submission.*
