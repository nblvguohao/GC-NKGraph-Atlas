"""
GC-NKGraph-Atlas GEO data downloader.

Downloads:
  - GSE62254 / ACRG (gastric cancer bulk, external validation)
  - GSE84437 (gastric cancer bulk, external validation)
  - GSE15459, GSE26942 (optional gastric cancer bulk)
  - GSE246662 (gastric cancer liver metastasis scRNA)

Requires: GEOparse (Python) or wget + manual download.
"""

import os
import sys
import argparse
from pathlib import Path
import subprocess

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.log_utils import Logger
from src.common.io_utils import load_config, ensure_dir


GEO_DATASETS = {
    # Bulk gastric cancer (external validation)
    "GSE62254": {
        "type": "bulk",
        "cancer": "gastric_cancer",
        "role": "external_validation",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62254",
        "expected_platform": "GPL570",
        "required_files": ["GSE62254_series_matrix.txt.gz"],
    },
    "GSE84437": {
        "type": "bulk",
        "cancer": "gastric_cancer",
        "role": "external_validation",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84437",
        "expected_platform": "GPL570",
        "required_files": ["GSE84437_series_matrix.txt.gz"],
    },
    # Optional gastric cancer bulk
    "GSE15459": {
        "type": "bulk",
        "cancer": "gastric_cancer",
        "role": "optional",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15459",
    },
    "GSE26942": {
        "type": "bulk",
        "cancer": "gastric_cancer",
        "role": "optional",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26942",
    },
    # scRNA (bridge dataset)
    "GSE246662": {
        "type": "scrna",
        "cancer": "gastric_cancer_liver_metastasis",
        "role": "candidate_gastric_and_liver_bridge",
        "status": "VERIFY_REQUIRED",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246662",
        "note": "gastric cancer liver metastasis scRNA with impaired NK function. Verify NK count before relying on it.",
        "required_files": [
            "GSE246662_RAW.tar",
        ],
    },
}

DEFAULT_OUTPUT_DIR_BULK = "data/raw/bulk"
DEFAULT_OUTPUT_DIR_SCRNA = "data/raw/scrna"


def check_geo_software() -> dict:
    """Check which GEO download methods are available."""
    available = {"geoquery_r": False, "geoparse": False, "wget": True}

    # Check GEOparse (Python)
    try:
        import GEOparse
        available["geoparse"] = True
    except ImportError:
        pass

    # Check for wget
    try:
        subprocess.run(["wget", "--version"], capture_output=True, check=True)
        available["wget"] = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        available["wget"] = False

    return available


def download_geo_series(
    geo_id: str,
    output_dir: str,
    logger: Logger,
    method: str = "auto",
):
    """Download a GEO series matrix and supplementary files.

    Args:
        geo_id: GEO accession (e.g., GSE62254)
        output_dir: Directory to save files
        logger: Logger instance
        method: 'auto', 'geoparse', 'wget', 'manual'
    """
    dataset_dir = os.path.join(output_dir, geo_id.lower())
    ensure_dir(dataset_dir)

    methods = check_geo_software()

    if method == "auto":
        if methods["geoparse"]:
            method = "geoparse"
        elif methods["wget"]:
            method = "wget"
        else:
            method = "manual"

    if method == "geoparse":
        import GEOparse
        try:
            gse = GEOparse.get_GEO(geo=geo_id, destdir=dataset_dir, silent=True)
            logger.ok(
                phase="DATA_DOWNLOAD",
                message=f"Downloaded {geo_id} via GEOparse to {dataset_dir}/",
                script=__file__,
            )
            # List platforms found
            platforms = list(gse.gpls.keys())
            logger.ok(
                phase="DATA_DOWNLOAD",
                message=f"  Platforms: {platforms}",
                script=__file__,
            )
        except Exception as e:
            logger.fail(
                phase="DATA_DOWNLOAD",
                message=f"GEOparse failed for {geo_id}: {e}\n  Falling back to manual download.",
                script=__file__,
            )
            method = "manual"

    if method == "wget":
        series_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:-3]}nnn/{geo_id}/matrix/{geo_id}_series_matrix.txt.gz"
        suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:-3]}nnn/{geo_id}/suppl/{geo_id}_RAW.tar"

        for url, desc in [(series_url, "series matrix"), (suppl_url, "supplementary")]:
            out_path = os.path.join(dataset_dir, os.path.basename(url))
            if not os.path.exists(out_path):
                cmd = f"wget -q -O {out_path} {url}"
                ret = subprocess.run(cmd, shell=True).returncode
                if ret == 0:
                    logger.ok(
                        phase="DATA_DOWNLOAD",
                        message=f"Downloaded {desc} for {geo_id}",
                        script=__file__,
                    )
                else:
                    logger.needs_review(
                        phase="DATA_DOWNLOAD",
                        message=f"Could not download {desc} for {geo_id} from:\n  {url}",
                        script=__file__,
                    )

    if method == "manual":
        dataset_info = GEO_DATASETS.get(geo_id, {})
        logger.needs_review(
            phase="DATA_DOWNLOAD",
            message=f"Manual download required for {geo_id}\n"
                    f"  URL: {dataset_info.get('url', 'N/A')}\n"
                    f"  Target dir: {dataset_dir}",
            script=__file__,
        )


def main():
    parser = argparse.ArgumentParser(description="Download GEO datasets for GC-NKGraph-Atlas")
    parser.add_argument("--datasets", nargs="+", default=list(GEO_DATASETS.keys()),
                       help=f"Datasets to download (default: all)")
    parser.add_argument("--method", default="auto", choices=["auto", "geoparse", "wget", "manual"])
    parser.add_argument("--only-scrna", action="store_true", help="Download only scRNA datasets")
    parser.add_argument("--only-bulk", action="store_true", help="Download only bulk datasets")
    args = parser.parse_args()

    logger = Logger()
    bulk_dir = ensure_dir(DEFAULT_OUTPUT_DIR_BULK)
    scrna_dir = ensure_dir(DEFAULT_OUTPUT_DIR_SCRNA)

    for geo_id in args.datasets:
        info = GEO_DATASETS.get(geo_id)
        if info is None:
            logger.fail(phase="DATA_DOWNLOAD", message=f"Unknown dataset: {geo_id}")
            continue

        if args.only_scrna and info["type"] != "scrna":
            continue
        if args.only_bulk and info["type"] != "bulk":
            continue

        type_log = info.get("status", "OK")
        logger.ok(
            phase="DATA_DOWNLOAD",
            message=f"Starting download: {geo_id} ({info['cancer']}, {info['type']}) | {type_log}",
            script=__file__,
        )

        if info["type"] == "scrna":
            download_geo_series(geo_id, str(scrna_dir), logger, args.method)
        else:
            download_geo_series(geo_id, str(bulk_dir), logger, args.method)

    print(f"\nDone. See {logger.log_path} for details.")


if __name__ == "__main__":
    main()
