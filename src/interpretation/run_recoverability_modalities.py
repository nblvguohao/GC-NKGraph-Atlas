"""Emit only available direct-modality evidence for the recoverability atlas."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.interpretation.recoverability_modalities import not_measured


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=ROOT / "submission_bundle_BiB/03_supplementary/tables/recoverability_direct_modality.tsv")
    args = parser.parse_args()
    rows = [
        {"card_id": "zheng_nk_sm_topology", **not_measured("metabolomics", "serine/sphingomyelin abundance", "MTBLS3303 public archive not retrieved", "MTBLS3303")},
        {"card_id": "adenosine_nk_suppression", **not_measured("metabolomics", "adenosine abundance", "MTBLS3303 public archive not retrieved", "MTBLS3303")},
        {"card_id": "tgfb_nk_exclusion", **not_measured("spatial_transcriptomics", "CAF--NK spot-module adjacency", "GSE251950 archive not retrieved", "GSE251950")},
        {"card_id": "nkg2d_micab_shedding", **not_measured("protein", "MICA/B abundance", "GSE122401 public archive contains RNA matrices but no sample-level protein matrix", "GSE122401")},
    ]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} ({len(rows)} direct-evidence rows)")


if __name__ == "__main__":
    main()
