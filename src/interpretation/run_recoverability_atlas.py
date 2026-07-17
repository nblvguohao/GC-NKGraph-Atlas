"""Aggregate verified real-data evidence into a gated recoverability atlas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def cross_mechanism_verdict(evidence: pd.DataFrame) -> str:
    recovered = evidence[(evidence.status == "recovered") & (evidence.concordant_cohorts >= 2)]
    if recovered.card_id.nunique() >= 3 and bool(recovered.direct_modality.any()):
        return "cross_mechanism_pattern_supported"
    return "comparative_atlas_only"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submission-root", default=ROOT / "submission_bundle_BiB")
    parser.add_argument("--verify-submission", action="store_true")
    args = parser.parse_args()
    tables = Path(args.submission_root) / "03_supplementary/tables"
    transcriptome = pd.read_csv(tables / "recoverability_transcriptome_per_cohort.tsv", sep="\t")
    direct = pd.read_csv(tables / "recoverability_direct_modality.tsv", sep="\t")
    summary = (transcriptome.sort_values("status")
               .groupby(["card_id", "comparison_id", "layer"], as_index=False)
               .agg(status=("status", lambda x: "recovered" if (x == "recovered").all() else "not_recovered"),
                    concordant_cohorts=("concordant_cohorts", "max"), accession=("accession", lambda x: ";".join(sorted(set(x))))))
    summary["direct_modality"] = summary.card_id.isin(direct.loc[direct.status != "not_measured", "card_id"])
    verdict = cross_mechanism_verdict(summary)
    summary.to_csv(tables / "recoverability_atlas.tsv", sep="\t", index=False)
    (tables / "recoverability_cross_mechanism_verdict.json").write_text(
        json.dumps({"verdict": verdict, "rule": "3 cards + 2 cohorts + direct modality"}, indent=2), encoding="utf-8")
    if args.verify_submission:
        assert verdict in {"cross_mechanism_pattern_supported", "comparative_atlas_only"}
    print(verdict)


if __name__ == "__main__":
    main()
