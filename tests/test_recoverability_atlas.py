import pandas as pd

from src.interpretation.run_recoverability_atlas import cross_mechanism_verdict


def test_global_pattern_needs_three_cards_two_cohorts_and_direct_evidence():
    evidence = pd.DataFrame([
        {"card_id": "a", "status": "recovered", "concordant_cohorts": 2, "direct_modality": True},
        {"card_id": "b", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
        {"card_id": "c", "status": "recovered", "concordant_cohorts": 2, "direct_modality": False},
    ])
    assert cross_mechanism_verdict(evidence) == "cross_mechanism_pattern_supported"
