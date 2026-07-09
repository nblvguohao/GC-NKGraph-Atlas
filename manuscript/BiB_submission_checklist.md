# Briefings in Bioinformatics — Submission Checklist

> Doubles as a readiness gate. Status: ✅ done · ⚠️ partial · ⬜ pending (server or
> writing). Verify every item against the journal's current *Instructions to
> Authors* at submission time — Oxford requirements change.
> Last updated: 2026-07-07.

## A. Article type & scope
- [x] ✅ Article type chosen: **Problem Solving Protocol** (original method +
  benchmarking). Confirm this is still the intended type vs. a standard research
  article at submission.
- [x] ✅ Scope fit: computational method of broad interest to bioinformatics
  readership; benchmarked against baselines; transferable abstraction.

## B. Manuscript structure
- [x] ✅ Title (≤ ~150 chars incl. spaces recommended) — present.
- [x] ✅ **Key Points** section (3–5 bullets) — required by BiB; present, reframed.
- [x] ✅ Structured/clear Abstract with Motivation/Results/Availability/Contact.
- [x] ✅ Introduction, Methods, Results, Discussion, Conclusion.
- [x] ✅ **Related work (§1.4)** — written (Scissor, CellChat, NicheNet, scFEA,
  CytoTRACE, CIBERSORTx/quanTIseq, MOGONET, HGT), refs [17]–[24] added.
- [x] ✅ **Methods §2.4 (scRNA)** — written to match `run_scrna_v2.py` +
  `qc_filter.py` (loading, QC thresholds, scVI, Leiden, marker annotation).
- [x] ✅ Running head present (GC-NKGraph-Atlas).

## C. Results completeness (the real gate)
- [x] ✅ Table 1 (datasets) — filled with real n.
- [x] ✅ §3.2 + Table 2 (positive-control recovery, multi-resolution) — filled.
- [x] ✅ §3.5 + Table 4 (de-circularized tumor-intrinsic candidates) — filled.
- [x] ✅ GNN internal CV numbers (§3.4) — filled.
- [x] ✅ **Table 3 baseline rows** — baselines pulled from A100, merged on identical
  folds, paired tests done; Table 3 + §3.4 filled (GNN on par with LightGBM/XGBoost,
  beats ElasticNet/SVM/MLP).
- [x] ✅ **External validation (§3.3, Table 5)** — GEO probe remapping done on A100
  (`run_geo_external_validation.py`); effector coupling replicates in both cohorts
  (r=0.42/0.62, p≪1e-13); NK markers 6/7 and 7/7.
- [x] ✅ **Figures 1–5** — publication-quality figures generated
  (`src/figures/make_figures.py`, Okabe-Ito CVD-safe palette, PDF+PNG): Fig1 Arm A
  recovery, Fig2 Arm B + external validation, Fig3 targets, Fig4 model comparison,
  Fig5 mechanism-card concept diagram. All numbers read from result tables; visually
  inspected.
- [ ] ⚠️ Ablation (§3.7) — either run a minimal ablation (graph edges on/off) or
  remove the placeholder heading.

## D. Front matter & declarations (BiB-required)
- [x] ✅ Author list + affiliations — matches `作者信息和基金.txt`.
- [x] ✅ Corresponding authors marked: L. Gu (glc@ahau.edu.cn), A. Zhou
  (zhouailian@caas.cn).
- [x] ✅ ORCID 占位字段已加（7 作者，待各自填 16 位 iD）。
- [x] ✅ **Author Contributions (CRediT)** statement — added (verify individual
  roles reflect actual contributions before submission).
- [x] ✅ **Funding** statement — present (NSFC 32472007, 62301006, 62301008; Anhui
  2308085MF217, 2308085QF202). Verify grant numbers/typos against award letters.
- [x] ✅ **Competing Interests** statement — added to manuscript.
- [x] ✅ **Data Availability** statement — present (public cohorts + code).
- [x] ✅ **Code availability**: repo structure ready for clean-clone reproduction —
  `README.md` (staged reproduction commands + synthetic mode), `requirements.txt`,
  `environment.yml` (GEOparse added), `.gitignore` fixed so the paper snapshot is
  tracked. `[repository URL]` replaced with https://github.com/nblvguohao/GC-NKGraph-Atlas
  throughout; **MIT `LICENSE` added.** **Remaining:** actually create/push the public repo.
- [x] ✅ Ethics statement — added ("Not applicable"; public de-identified data).

## E. References & formatting
- [x] ✅ Anchor + core method references present (16). 
- [x] ✅ References expanded to [24]; **[17]–[24] verified online** against primary
  sources (all volume/page/year correct as written, 2026-07-07).
- [x] ✅ References in BiB numbered Vancouver format [1]–[24]. Verify all against
  PubMed primary records.
- [x] ✅ Figures: vector PDFs, ≥300 dpi PNG, self-contained legends.
- [x] ✅ Supplementary: key tables in `submission_package/02_supplementary_tables/`.
  ⚠️ Add a Supplementary Methods note + index before submission.

## F. Submission-system items
- [x] ✅ Cover letter finalized (date July 7, 2026; 3 reviewers with verified emails).
- [x] ✅ Suggested reviewers (Satu Mustjoki, Fabian J. Theis, Xia Li) with emails.
- [x] ✅ Author biographies (~30 words each, 7 authors) added to manuscript.
- [x] ✅ ORCID iD placeholders added (7 authors, [0000-...] template).
- [x] ✅ Keywords ≤ 6 (condensed from 8).
- [x] ✅ Title page / abstract / key-points in single-block format (BiB single-blind).
- [ ] ⬜ ORCID iDs: each author fills their 16-digit iD at orcid.org.
- [ ] ⬜ Register corresponding authors in Oxford ScholarOne system.
- [ ] ⬜ Verify no dual submission; all authors approve.
- [ ] ⬜ Highlights / graphical abstract if requested by the editorial office.

## G. Honesty & reproducibility guardrails (project-specific, keep)
- [x] ✅ Every axis claim bounded by "transcriptional program permissive-of /
  associated-with" language (POSITIONING §6).
- [x] ✅ Positive-control outcome reported honestly per-hypothesis (Table 2), no
  cherry-picking.
- [x] ✅ Target list de-circularized; NK readout separated from tumor-intrinsic
  targets.
- [ ] ⬜ Pre-registration note: log the calibration change (metabolic edge now
  anchored on the effector arm, since the topology arm did not recover).

---

### The items that actually block submission (updated 2026-07-09)
1. ~~Table 3 baselines~~ — ✅ **done** (pulled from A100, merged, tested).
2. ~~External validation~~ — ✅ **done** (GEO remapped on A100; both cohorts replicate).
3. ~~Publication figures~~ — ✅ Fig 1–5 done (including mechanism-card concept diagram).
4. ~~References, ORCID placeholders, cover letter, author bios, keywords~~ — ✅ all done.
5. ~~Binary vs 4-class alignment~~ — ✅ **done** (Methods §2.4 now explains binary choice).
6. ~~Limitations reordering~~ — ✅ **done** (leads with proxy≠topology and no experimental validation).
7. **Public code repository** — the ONLY remaining blocker (push to GitHub).
8. **Each author fills ORCID iD** — individual author action.
9. **Grant number / CRediT verification** — corresponding author action.

See `manuscript/PRE_SUBMISSION_CHECKLIST.md` for the complete per-author action list.
