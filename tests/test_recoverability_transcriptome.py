import pandas as pd

from src.interpretation.run_recoverability_transcriptome import assign_recovery_status


def test_recovered_needs_direction_fdr_coverage_and_two_cohorts():
    row = pd.Series({"direction_ok": True, "p_fdr": 0.01, "coverage": 1.0, "concordant_cohorts": 2})
    assert assign_recovery_status(row) == "recovered"
    one_cohort = row.copy()
    one_cohort["concordant_cohorts"] = 1
    assert assign_recovery_status(one_cohort) == "partially_recovered"
