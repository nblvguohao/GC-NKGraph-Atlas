import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERDICT = ROOT / "submission_bundle_BiB/03_supplementary/tables/recoverability_cross_mechanism_verdict.json"
MANUSCRIPT = ROOT / "submission_bundle_BiB/01_manuscript/main_manuscript.md"


def test_universal_claim_requires_supported_verdict():
    verdict = json.loads(VERDICT.read_text(encoding="utf-8"))["verdict"]
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    if verdict != "cross_mechanism_pattern_supported":
        assert "general rule across mechanisms" not in manuscript
        assert "universal recoverability law" not in manuscript
