import pandas as pd
from pathlib import Path

from src.interpretation.run_recoverability_atlas import (
    ROOT,
    cross_mechanism_verdict,
    eligible_direct_cards,
    write_source_manifest,
)


def test_global_pattern_needs_three_cards_two_cohorts_and_direct_evidence():
    evidence = pd.DataFrame([
        {"card_id": "a", "status": "recovered", "concordant_cohorts": 2, "direct_modality": True},
        {"card_id": "b", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
        {"card_id": "c", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
    ])
    assert cross_mechanism_verdict(evidence) == "cross_mechanism_pattern_supported"


def test_source_manifest_includes_outer_pending_and_four_available_gsm_assets(tmp_path):
    output = tmp_path / "recoverability_source_manifest.tsv"
    manifest = ROOT / "configs/recoverability_atlas/real_data_manifest.yaml"

    table = write_source_manifest(manifest, output)

    outer = table.loc[table["asset"] == "gse251950"].iloc[0]
    assert outer["status"] == "pending_download"
    per_gsm = table.loc[table["asset"].str.startswith("gse251950_gsm")]
    assert len(per_gsm) == 4
    assert set(per_gsm["status"]) == {"available"}
    assert per_gsm["sha256"].str.len().eq(64).all()
    assert set(per_gsm["analysis_scope"]) == {"exploratory_four_verified_per_gsm_subset"}


def test_exploratory_spatial_rows_do_not_count_as_gate_eligible_direct_evidence():
    direct = pd.DataFrame([
        {"card_id": "tgfb", "status": "measured", "scope": "exploratory_four_verified_per_gsm_subset"},
        {"card_id": "nkg2d", "status": "measured", "scope": "confirmatory"},
    ])

    assert eligible_direct_cards(direct) == {"nkg2d"}
