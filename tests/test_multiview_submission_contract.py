"""Submission-facing claim and artifact contracts for the multiview audit."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "submission_bundle_BiB" / "01_manuscript" / "main_manuscript.md"
SUPPLEMENT_INDEX = ROOT / "submission_bundle_BiB" / "03_supplementary" / "SUPPLEMENTARY_INDEX.md"
TABLES = ROOT / "submission_bundle_BiB" / "03_supplementary" / "tables"


def test_no_gain_verdict_is_reported_with_pre_registered_conservative_language():
    comparisons = pd.read_csv(TABLES / "multiview_bootstrap_comparisons.tsv", sep="\t")
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    normalized = " ".join(manuscript.replace("**", "").split())

    assert set(comparisons["verdict"]) == {"no_stable_external_gain"}
    assert "view separation improves structural auditability but does not establish a predictive advantage" in normalized
    assert "stable external predictive gain" not in manuscript


def test_methods_document_masking_external_isolation_and_seed_uncertainty():
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    normalized = " ".join(manuscript.replace("**", "").split())

    for phrase in (
        "label-masked primary analysis",
        "TCGA-LIHC was never used for standardization",
        "ten fixed random seeds",
        "2,000 paired stratified bootstrap",
        "TREE [49] and GRAFT [50]",
    ):
        assert phrase in normalized


def test_supplement_indexes_every_submission_facing_multiview_artifact():
    index = SUPPLEMENT_INDEX.read_text(encoding="utf-8")
    required = {
        "multiview_benchmark_summary.tsv",
        "multiview_bootstrap_comparisons.tsv",
        "multiview_weight_stability.tsv",
        "multiview_leave_one_out_bootstrap.tsv",
        "multiview_mechanism_randomization_summary.tsv",
        "figS1_multiview",
    }

    for name in required:
        assert name in index
    for name in required - {"figS1_multiview"}:
        assert (TABLES / name).exists()


def test_supplementary_multiview_figure_has_both_formats():
    figures = ROOT / "submission_bundle_BiB" / "02_figures"

    for suffix in ("png", "pdf"):
        path = figures / f"figS1_multiview.{suffix}"
        assert path.exists() and path.stat().st_size > 0
