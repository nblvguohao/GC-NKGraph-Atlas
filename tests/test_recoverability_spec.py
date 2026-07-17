import pytest

from src.interpretation.recoverability_spec import load_comparisons


def test_confirmatory_comparison_requires_explicit_sign(tmp_path):
    spec = tmp_path / "spec.yaml"
    spec.write_text(
        "comparisons: [{card_id: x, comparison_id: y, layer: z, "
        "left_module: a, right_module: b, expected_sign: NEEDS_REVIEW}]\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="explicit expected_sign"):
        load_comparisons(spec)
