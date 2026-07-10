"""
Download GEO datasets for GC-NKGraph-Atlas.

Downloads:
  - GSE62254 (gastric cancer bulk, external validation)
  - GSE84437 (gastric cancer bulk, external validation)
  - GSE246662 (gastric cancer liver metastasis scRNA)
"""

import os
import sys

# Bypass SOCKS proxy
os.environ["no_proxy"] = "*"
os.environ["all_proxy"] = ""
os.environ["ALL_PROXY"] = ""

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir
import GEOparse

GEO_DATASETS = {
    "GSE62254": {
        "type": "bulk",
        "dest": "data/raw/bulk/gse62254",
    },
    "GSE84437": {
        "type": "bulk",
        "dest": "data/raw/bulk/gse84437",
    },
    "GSE246662": {
        "type": "scrna",
        "dest": "data/raw/scrna/gse246662",
    },
}


def download_geo(geo_id: str, dest_dir: str) -> bool:
    """Download a GEO dataset using GEOparse."""
    print(f"\n{'='*60}")
    print(f"Downloading {geo_id}...")
    print(f"  Destination: {dest_dir}")

    ensure_dir(dest_dir)

    try:
        gse = GEOparse.get_GEO(geo=geo_id, destdir=dest_dir, silent=True)
        title = gse.metadata.get('title', ['Unknown'])[0]
        n_samples = len(gse.gsms)
        print(f"  Title: {title}")
        print(f"  Samples: {n_samples}")
        print(f"  Status: OK")
        return True
    except Exception as e:
        print(f"  Status: FAILED - {e}")
        return False


def main():
    results = {}
    for geo_id, info in GEO_DATASETS.items():
        ok = download_geo(geo_id, info["dest"])
        results[geo_id] = ok

    print(f"\n{'='*60}")
    print("DOWNLOAD SUMMARY:")
    for geo_id, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  {geo_id}: {status}")


if __name__ == "__main__":
    main()
