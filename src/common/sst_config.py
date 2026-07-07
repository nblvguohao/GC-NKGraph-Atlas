"""
GC-NKGraph-Atlas Shared SST-Axis Configuration Loader.

Single source of truth for all SST (Serine-Sphingomyelin-Topology) gene modules.
All modules that need SST gene lists MUST import from here — never hard-code gene
modules in individual scripts.

Usage:
    from src.common.sst_config import load_sst_modules, get_sst_genes, SST_MODULES

    modules = load_sst_modules()          # from configs/sst_axis_config.yaml
    all_genes = get_sst_genes(modules)    # flat set of all SST genes
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Set

import yaml

# ---------------------------------------------------------------------------
# Default config path relative to project root
# ---------------------------------------------------------------------------
_DEFAULT_CONFIG_PATH = "configs/sst_axis_config.yaml"
_MECHANISM_CARD_PATH = "configs/mechanism_cards/zheng_nk_sm_topology.yaml"

# ---------------------------------------------------------------------------
# Cached module definitions (lazy-loaded)
# ---------------------------------------------------------------------------
_cached_modules: Optional[Dict[str, Dict]] = None


def _resolve_project_root() -> Path:
    """Find the project root (where configs/ lives)."""
    # Walk up from this file's directory
    current = Path(__file__).resolve().parent
    # src/common -> src -> project root
    return current.parent.parent


def load_sst_modules(config_path: Optional[str] = None) -> Dict[str, Dict]:
    """
    Load SST-axis gene modules from the canonical YAML config.

    Returns a dict keyed by module name, each value a dict with:
        - role: str           (tumor_side | immune_side | outcome_anchor | therapeutic_hook)
        - cell_type_attribution: str
        - genes: List[str]
        - expected_direction: str
    """
    global _cached_modules

    if _cached_modules is not None and config_path is None:
        return _cached_modules

    if config_path is None:
        config_path = _DEFAULT_CONFIG_PATH

    config_file = _resolve_project_root() / config_path
    if not config_file.exists():
        raise FileNotFoundError(
            f"SST axis config not found: {config_file}. "
            f"Ensure {_DEFAULT_CONFIG_PATH} exists in the project root."
        )

    with open(config_file, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    modules: Dict[str, Dict] = {}
    for name, mod_data in raw.get("sst_axis", {}).get("modules", {}).items():
        modules[name] = {
            "role": mod_data.get("role", "unknown"),
            "cell_type_attribution": mod_data.get("cell_type_attribution", "unknown"),
            "genes": list(mod_data.get("genes", [])),
            "expected_direction": mod_data.get("expected_direction", "unknown"),
        }

    _cached_modules = modules
    return modules


def get_sst_genes(modules: Optional[Dict[str, Dict]] = None) -> Set[str]:
    """Return the flat set of all SST-axis gene symbols."""
    if modules is None:
        modules = load_sst_modules()
    all_genes: Set[str] = set()
    for mod in modules.values():
        for g in mod["genes"]:
            all_genes.add(g)
    return all_genes


def get_module_for_gene(gene: str, modules: Optional[Dict[str, Dict]] = None) -> Optional[str]:
    """Return the module name for a given gene, or None."""
    if modules is None:
        modules = load_sst_modules()
    for mod_name, mod_data in modules.items():
        if gene in mod_data["genes"]:
            return mod_name
    return None


def get_genes_by_role(role: str, modules: Optional[Dict[str, Dict]] = None) -> List[str]:
    """Return all genes for modules matching a given role."""
    if modules is None:
        modules = load_sst_modules()
    genes: List[str] = []
    for mod in modules.values():
        if mod["role"] == role:
            genes.extend(mod["genes"])
    return genes


def load_mechanism_card(card_id: str = "zheng_nk_sm_topology") -> Dict:
    """Load a full mechanism card from configs/mechanism_cards/."""
    config_file = _resolve_project_root() / "configs" / "mechanism_cards" / f"{card_id}.yaml"
    if not config_file.exists():
        raise FileNotFoundError(f"Mechanism card not found: {config_file}")
    with open(config_file, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Convenience: pre-load the default modules at import time for backward compat
# ---------------------------------------------------------------------------
SST_MODULES = load_sst_modules()
SST_GENES = get_sst_genes(SST_MODULES)
