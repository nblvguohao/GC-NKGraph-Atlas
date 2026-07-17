from src.interpretation.recoverability_modalities import not_measured


def test_missing_protein_matrix_is_not_measured_not_negative():
    row = not_measured("protein", "MICA abundance", "public matrix unavailable")
    assert row["status"] == "not_measured"
    assert row["not_measured_reason"] == "public matrix unavailable"
