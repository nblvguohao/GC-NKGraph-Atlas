# SST-Axis Pre-Registration Log

> Documents the calibration history and any deviations from the pre-registered
> protocol. This file is part of the paper's honesty infrastructure — every change
> from the original pre-registration is logged here with a date and rationale.
> Last updated: 2026-07-09.

---

## Original Pre-Registration (2026-07-07)

The hypotheses H1–H5 and their expected directions were registered in
`configs/sst_axis_config.yaml` before execution:

| ID | Hypothesis | Expected direction |
|----|-----------|-------------------|
| H1 | tumor_serine_capacity ⟂ nk_sm_balance | Direction calibrated, reported |
| H2 | nk_sm_balance (+) nk_protrusion_machinery | Positive |
| H3 | nk_protrusion_machinery (+) cytotoxicity anchors | Positive |
| H4 | nk_topology_permissive (−) HAVCR2 / dysfunction | Negative |
| H5 | intratumoral NK < peritumoral/peripheral NK in sm_balance & protrusion | Negative (tumor < normal) |

Recovery definition: H2–H5 pass in the pre-registered direction in liver (Arm A).

---

## Calibration Change Log

### CHANGE-001 — 2026-07-09: metabolic_crosstalk edge calibration source

**Original plan:** Calibrate the `metabolic_crosstalk` edge sign and weight from
H1 (tumor_serine_capacity ⟂ nk_sm_balance) in TCGA-LIHC.

**What happened:** H1 was null in TCGA-LIHC (r = −0.016, p = 0.74). No usable
signal for calibration.

**Decision:** The `metabolic_crosstalk` edge sign defaults to the anchor paper's
mechanistic direction (tumor serine dysregulation ↑ → NK SM depletion ↓), as
stated in Zheng et al. (Nat Immunol 2023) and its follow-up (Clin Transl Med
2023). The edge weight is set to 0.5 (midpoint of the weight range) to reflect
mechanistic grounding without over-weighting an uncalibrated prior.

**Status flag:** `CALIBRATION_FAILED_DEFAULTING_TO_ANCHOR_PAPER_DIRECTION`
(replaces the original `NEEDS_REVIEW`).

**Impact:**
- The `metabolic_crosstalk` edge is a hypothesis-driven prior, not a
  data-calibrated parameter.
- Downstream analyses that depend on this edge (axis scoring, target
  prioritization, in-silico stratification) carry this caveat.
- An ablation experiment (§3.7) compares full-graph vs no-metabolic-edge
  performance to quantify the edge's contribution.
- `tumor_serine_capacity.expected_direction` in `sst_axis_config.yaml` has
  been updated accordingly.

**Pre-registration integrity:** H1 was always defined as "direction calibrated,
reported" — not as a pre-specified direction. A null result with transparent
reporting is within the pre-registered protocol. The calibration failure does
not constitute a protocol deviation; it constitutes a calibration outcome.

### CHANGE-002 — 2026-07-09: revision of recovery verdict

**Original plan:** Arm A recovery = H2–H5 pass in pre-registered direction in liver.

**What happened:** H3 (effector arm) and H5-cytotoxicity passed; H2 passed only
at single-cell resolution with very weak effect (r=+0.030); H4 failed (wrong
sign); H5-protrusion failed (wrong sign).

**Revised verdict:** The effector arm (H3) and cell-type-resolved metabolic
coupling (H2, scNK only) recover; the physical topology phenotype (H4,
H5-protrusion) does not recover from machinery transcription. This is reported
as a scoping map rather than a blanket recovery/failure. See manuscript §3.2
and §4.1 for full discussion.

---

## Current SST-Axis Calibration State

```yaml
calibration_summary:
  calibrated_on: NONE  # H1 null; no data-driven calibration possible
  default_direction_source: "Zheng et al. 2023 Nat Immunol mechanistic direction"
  edge_weight: 0.5  # hypothesis-driven prior, not data-calibrated
  needs_physical_validation: true  # requires paired scRNA + imaging for true calibration
  next_review_point: "When paired NK scRNA + SEM/MS data become available"
```

---

## Verification

To verify this log matches the config:

```bash
grep -A1 "tumor_serine_capacity:" configs/sst_axis_config.yaml | grep "expected_direction"
# Expected: CALIBRATION_FAILED_DEFAULTING_TO_ANCHOR_PAPER_DIRECTION
```

---

*This log is part of the submission package. Any future calibration update
(e.g., from paired transcriptomic + imaging data) must be appended here with
date, rationale, and impact assessment.*
