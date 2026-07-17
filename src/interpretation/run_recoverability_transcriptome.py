"""Evaluate pre-specified mechanism-card correlations on verified real cohorts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.real_data import assert_real_asset, load_real_data_manifest
from src.interpretation.recoverability_spec import load_comparisons


ROOT = Path(__file__).resolve().parents[2]
CARD_DIR = ROOT / "configs" / "mechanism_cards"


def assign_recovery_status(row: pd.Series) -> str:
    if bool(row.direction_ok) and row.p_fdr < 0.05 and row.coverage >= 0.7:
        return "recovered" if row.concordant_cohorts >= 2 else "partially_recovered"
    return "not_recovered"


def _bh_fdr(pvalues: pd.Series) -> pd.Series:
    order = pvalues.rank(method="first").astype(int).to_numpy() - 1
    values = pvalues.to_numpy(dtype=float)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = pd.Series(adjusted[::-1]).cummin().iloc[::-1].clip(upper=1).to_numpy()
    result = pd.Series(index=pvalues.index, dtype=float)
    result.iloc[order] = adjusted
    return result


def _modules() -> dict[str, dict[str, list[str]]]:
    result = {}
    for path in CARD_DIR.glob("*.yaml"):
        if path.name in {"registry.yaml", "mechanism_card.template.yaml"}:
            continue
        with path.open(encoding="utf-8") as handle:
            card = yaml.safe_load(handle)["mechanism_card"]
        result[card["id"]] = {m["name"]: m["genes"] for m in card["transcriptional_proxy"]["modules"]}
    return result


def _score(expression: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, float]:
    found = [gene for gene in genes if gene in expression.columns]
    coverage = len(found) / len(genes)
    if not found:
        return pd.Series(float("nan"), index=expression.index), coverage
    x = expression[found].apply(pd.to_numeric, errors="coerce")
    return ((x - x.mean()) / x.std(ddof=0).replace(0, 1)).mean(axis=1), coverage


def evaluate_transcriptomic_recovery(manifest_path: str | Path) -> pd.DataFrame:
    assets = load_real_data_manifest(manifest_path)
    comparisons = load_comparisons(ROOT / "configs" / "recoverability_atlas" / "card_layer_spec.yaml")
    modules = _modules()
    cohort_assets = [("TCGA-STAD", "tcga_stad_expression"), ("TCGA-LIHC", "tcga_lihc_expression"),
                     ("GSE62254", "gse62254_expression"), ("GSE84437", "gse84437_expression")]
    rows = []
    for cohort, asset_name in cohort_assets:
        asset = assets[asset_name]
        assert_real_asset(asset, asset.local_path)
        expression = pd.read_csv(asset.local_path, sep="\t", index_col=0)
        for comp in comparisons:
            left, c_left = _score(expression, modules[comp.card_id][comp.left_module])
            right, c_right = _score(expression, modules[comp.card_id][comp.right_module])
            valid = left.notna() & right.notna()
            r, p = pearsonr(left[valid], right[valid])
            rows.append({"card_id": comp.card_id, "comparison_id": comp.comparison_id, "layer": comp.layer,
                         "cohort": cohort, "accession": asset.accession, "r": r, "p": p,
                         "coverage": min(c_left, c_right), "expected_sign": comp.expected_sign,
                         "direction_ok": r > 0 if comp.expected_sign == "positive" else r < 0})
    out = pd.DataFrame(rows)
    out["p_fdr"] = _bh_fdr(out.p)
    out["concordant_cohorts"] = out.groupby("comparison_id")["direction_ok"].transform("sum")
    out["status"] = out.apply(assign_recovery_status, axis=1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=ROOT / "configs/recoverability_atlas/real_data_manifest.yaml")
    parser.add_argument("--output", default=ROOT / "submission_bundle_BiB/03_supplementary/tables/recoverability_transcriptome_per_cohort.tsv")
    args = parser.parse_args()
    out = evaluate_transcriptomic_recovery(args.manifest)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} ({len(out)} rows)")


if __name__ == "__main__":
    main()
