"""Emit only available direct-modality evidence for the recoverability atlas."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.interpretation.recoverability_modalities import (
    not_measured,
    run_visium_spot_module_adjacency,
)


ROOT = Path(__file__).resolve().parents[2]


def write_direct_modality_table(
    output: Path,
    *,
    archive_dir: Path = ROOT / "data/external/recoverability/GSE251950/per_gsm",
    n_permutations: int = 1000,
    seed: int = 20260717,
) -> pd.DataFrame:
    """Write direct evidence, preserving real spatial rows and explicit missingness."""
    rows = [
        {"card_id": "zheng_nk_sm_topology", **not_measured("metabolomics", "serine/sphingomyelin abundance", "MTBLS3303 public archive not retrieved", "MTBLS3303")},
        {"card_id": "adenosine_nk_suppression", **not_measured("metabolomics", "adenosine abundance", "MTBLS3303 public archive not retrieved", "MTBLS3303")},
        {"card_id": "nkg2d_micab_shedding", **not_measured("protein", "MICA/B abundance", "GSE122401 public archive contains RNA matrices but no sample-level protein matrix", "GSE122401")},
    ]
    spatial = run_visium_spot_module_adjacency(
        sorted(Path(archive_dir).glob("GSM*.tar.gz")),
        n_permutations=n_permutations,
        seed=seed,
    )
    table = pd.concat([pd.DataFrame(rows), spatial], ignore_index=True, sort=False)
    # A submission TSV must expose non-applicability explicitly rather than hide
    # it behind blank trailing columns (which can be misread after import).
    table = table.replace("", "not_applicable").fillna("not_applicable")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=ROOT / "submission_bundle_BiB/03_supplementary/tables/recoverability_direct_modality.tsv")
    parser.add_argument("--spatial-archive-dir", default=ROOT / "data/external/recoverability/GSE251950/per_gsm")
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260717)
    args = parser.parse_args()
    table = write_direct_modality_table(
        Path(args.output), archive_dir=Path(args.spatial_archive_dir),
        n_permutations=args.n_permutations, seed=args.seed,
    )
    print(f"wrote {args.output} ({len(table)} direct-evidence rows)")


if __name__ == "__main__":
    main()
