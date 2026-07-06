"""
GC-NKGraph-Atlas I/O utilities.

Standardized file reading/writing with provenance tracking.
"""

import os
import sys
import json
import yaml
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd


def load_config(config_path: str) -> dict:
    """Load a YAML config file with path resolution."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with open(path, "r") as f:
        return yaml.safe_load(f)


def save_table(
    df: pd.DataFrame,
    path: str,
    provenance: Optional[Dict[str, str]] = None,
    index: bool = False,
):
    """Save a TSV with provenance metadata as markdown comment."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if provenance is None:
        provenance = {}

    provenance.update({
        "script": sys.argv[0],
        "timestamp": __import__("datetime").datetime.now().isoformat(),
    })

    # Write provenance header as YAML-style comments
    with open(path, "w", encoding="utf-8") as f:
        f.write("# provenance:\n")
        for k, v in provenance.items():
            f.write(f"#   {k}: {v}\n")

    # Append TSV data
    df.to_csv(path, sep="\t", index=index, mode="a")


def load_table(path: str, skip_provenance: bool = True) -> pd.DataFrame:
    """Load a TSV, optionally skipping provenance comment lines."""
    if skip_provenance:
        return pd.read_csv(path, sep="\t", comment="#")
    return pd.read_csv(path, sep="\t")


def ensure_dir(path: str) -> Path:
    """Ensure a directory exists and return the Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_json(obj: Any, path: str):
    """Save JSON with pretty printing."""
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=str)


def load_json(path: str) -> Any:
    """Load JSON."""
    with open(path) as f:
        return json.load(f)
