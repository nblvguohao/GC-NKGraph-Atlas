"""Direct-modality evidence helpers; missing inputs are never imputed."""

from __future__ import annotations


def not_measured(modality: str, endpoint: str, reason: str, accession: str = "") -> dict:
    return {"modality": modality, "direct_endpoint": endpoint, "status": "not_measured",
            "not_measured_reason": reason, "accession": accession}
