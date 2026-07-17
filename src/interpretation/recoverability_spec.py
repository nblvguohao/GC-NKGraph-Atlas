"""Validated, pre-specified comparisons for the recoverability atlas."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Comparison:
    card_id: str
    comparison_id: str
    layer: str
    left_module: str
    right_module: str
    expected_sign: str
    requires_purity_control: bool = True


def load_comparisons(path: str | Path) -> list[Comparison]:
    with Path(path).open(encoding="utf-8") as handle:
        rows = (yaml.safe_load(handle) or {}).get("comparisons", [])
    comparisons = []
    for row in rows:
        if row.get("expected_sign") not in {"positive", "negative"}:
            raise ValueError("confirmatory comparison requires explicit expected_sign")
        required = ("card_id", "comparison_id", "layer", "left_module", "right_module")
        if any(not row.get(field) for field in required):
            raise ValueError("comparison has missing required field")
        comparisons.append(Comparison(**row))
    if not comparisons:
        raise ValueError("at least one comparison is required")
    return comparisons
