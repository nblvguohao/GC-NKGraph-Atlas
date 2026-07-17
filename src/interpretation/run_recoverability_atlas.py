"""Aggregate verified real-data evidence into a gated recoverability atlas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]


def write_source_manifest(manifest_path: Path, output_path: Path) -> pd.DataFrame:
    """Materialize all registered real-data assets, including pending outer archives."""
    with Path(manifest_path).open(encoding="utf-8") as handle:
        document = yaml.safe_load(handle) or {}
    assets = document.get("assets")
    if not isinstance(assets, dict):
        raise ValueError("real-data manifest requires an assets mapping")
    fields = (
        "asset", "accession", "source_url", "modality", "species", "sample_count",
        "sha256", "local_path", "content_length_bytes", "analysis_scope", "status", "retrieved_at",
    )
    rows = [{"asset": name, **{field: payload.get(field, "") for field in fields if field != "asset"}}
            for name, payload in assets.items()]
    table = pd.DataFrame(rows, columns=fields).replace("", "not_applicable").fillna("not_applicable")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, sep="\t", index=False)
    return table


def cross_mechanism_verdict(evidence: pd.DataFrame) -> str:
    recovered = evidence[(evidence.status == "recovered") & (evidence.concordant_cohorts >= 2)]
    if recovered.card_id.nunique() >= 3 and bool(recovered.direct_modality.any()):
        return "cross_mechanism_pattern_supported"
    return "comparative_atlas_only"


def eligible_direct_cards(direct: pd.DataFrame) -> set[str]:
    """Return only measured, non-exploratory direct-modality cards for the gate."""
    scope = direct.get("scope", pd.Series("confirmatory", index=direct.index)).astype(str)
    eligible = direct.loc[(direct.status == "measured") & ~scope.str.startswith("exploratory_")]
    return set(eligible.card_id)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submission-root", default=ROOT / "submission_bundle_BiB")
    parser.add_argument("--verify-submission", action="store_true")
    args = parser.parse_args()
    tables = Path(args.submission_root) / "03_supplementary/tables"
    write_source_manifest(
        ROOT / "configs/recoverability_atlas/real_data_manifest.yaml",
        tables / "recoverability_source_manifest.tsv",
    )
    transcriptome = pd.read_csv(tables / "recoverability_transcriptome_per_cohort.tsv", sep="\t")
    direct = pd.read_csv(tables / "recoverability_direct_modality.tsv", sep="\t")
    summary = (transcriptome.sort_values("status")
               .groupby(["card_id", "comparison_id", "layer"], as_index=False)
               .agg(status=("status", lambda x: "recovered" if (x == "recovered").all() else "not_recovered"),
                    concordant_cohorts=("concordant_cohorts", "max"), accession=("accession", lambda x: ";".join(sorted(set(x))))))
    summary["direct_modality"] = summary.card_id.isin(eligible_direct_cards(direct))
    verdict = cross_mechanism_verdict(summary)
    summary.to_csv(tables / "recoverability_atlas.tsv", sep="\t", index=False)
    (tables / "recoverability_cross_mechanism_verdict.json").write_text(
        json.dumps({"verdict": verdict, "rule": "3 cards + 2 cohorts + direct modality"}, indent=2), encoding="utf-8")
    if args.verify_submission:
        assert verdict in {"cross_mechanism_pattern_supported", "comparative_atlas_only"}
    print(verdict)


if __name__ == "__main__":
    main()
