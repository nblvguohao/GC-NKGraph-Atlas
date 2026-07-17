"""Contracts for formal analyses restricted to real public data."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

import yaml


_FORBIDDEN_PATH_TOKENS = ("synthetic", "mock", "demo")
_REQUIRED_FIELDS = ("accession", "source_url", "modality", "species", "sample_count", "local_path")


@dataclass(frozen=True)
class RealDataAsset:
    accession: str
    source_url: str
    modality: str
    species: str
    sample_count: int
    sha256: str | None
    local_path: str
    status: str = "available"
    retrieved_at: str | None = None


def _validate_asset(name: str, payload: dict[str, Any]) -> RealDataAsset:
    for field in _REQUIRED_FIELDS:
        if not payload.get(field):
            raise ValueError(f"asset '{name}' requires {field}")
    if payload["species"] != "Homo sapiens":
        raise ValueError(f"asset '{name}' must be Homo sapiens")
    if int(payload["sample_count"]) < 1:
        raise ValueError(f"asset '{name}' requires a positive sample_count")
    status = payload.get("status", "available")
    digest = payload.get("sha256")
    if status != "pending_download" and (not isinstance(digest, str) or len(digest) != 64):
        raise ValueError(f"asset '{name}' requires sha256 after retrieval")
    if digest is not None and (not isinstance(digest, str) or len(digest) != 64):
        raise ValueError(f"asset '{name}' has invalid sha256")
    return RealDataAsset(
        accession=str(payload["accession"]), source_url=str(payload["source_url"]),
        modality=str(payload["modality"]), species=str(payload["species"]),
        sample_count=int(payload["sample_count"]), sha256=digest.lower() if digest else None,
        local_path=str(payload["local_path"]), status=str(status),
        retrieved_at=payload.get("retrieved_at"),
    )


def load_real_data_manifest(path: str | Path) -> dict[str, RealDataAsset]:
    """Load and validate a formal-analysis input manifest."""
    with Path(path).open(encoding="utf-8") as handle:
        document = yaml.safe_load(handle) or {}
    assets = document.get("assets")
    if not isinstance(assets, dict) or not assets:
        raise ValueError("manifest requires a non-empty assets mapping")
    return {name: _validate_asset(name, payload or {}) for name, payload in assets.items()}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_real_asset(asset: RealDataAsset, path: str | Path) -> None:
    """Reject non-real inputs and unverified files before formal analysis."""
    input_path = Path(path)
    if any(token in str(input_path).lower() for token in _FORBIDDEN_PATH_TOKENS):
        raise ValueError("non-real input path is forbidden")
    if asset.species != "Homo sapiens" or asset.status == "pending_download":
        raise ValueError("asset does not satisfy real-data contract")
    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise FileNotFoundError(f"real data file is absent or empty: {input_path}")
    if not asset.sha256 or file_sha256(input_path) != asset.sha256:
        raise ValueError(f"SHA-256 mismatch for {asset.accession}: {input_path}")
