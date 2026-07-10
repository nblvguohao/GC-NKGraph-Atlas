"""
Download all missing datasets for GC-NKGraph-Atlas.
Handles retries, proxy bypass, and integrity checks.
"""

import os
import sys
import time
import gzip
import shutil
import urllib.request
from pathlib import Path

# Bypass SOCKS proxy
os.environ["no_proxy"] = "*"
os.environ["all_proxy"] = ""
os.environ["ALL_PROXY"] = ""

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.io_utils import ensure_dir


def download_file(url: str, dest: str, retries: int = 3, timeout: int = 300) -> bool:
    """Download a file with retry logic."""
    ensure_dir(str(Path(dest).parent))

    for attempt in range(1, retries + 1):
        try:
            print(f"    Attempt {attempt}/{retries}...", end=" ", flush=True)
            start = time.time()
            urllib.request.urlretrieve(url, dest)
            elapsed = time.time() - start
            size_mb = os.path.getsize(dest) / 1e6
            if size_mb < 0.001:  # Empty or error page
                os.remove(dest)
                raise ValueError(f"Downloaded file too small ({size_mb:.4f} MB)")
            print(f"OK ({size_mb:.1f} MB, {elapsed:.0f}s)")
            return True
        except Exception as e:
            print(f"FAILED: {e}")
            if attempt < retries:
                wait = attempt * 5
                print(f"    Waiting {wait}s before retry...")
                time.sleep(wait)
    return False


# =========================================================================
# GEO Bulk Datasets
# =========================================================================

GEO_BULK = {
    "GSE62254": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/matrix/GSE62254_series_matrix.txt.gz",
        "dest": "data/raw/bulk/gse62254/GSE62254_series_matrix.txt.gz",
    },
    "GSE84437": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84437/matrix/GSE84437_series_matrix.txt.gz",
        "dest": "data/raw/bulk/gse84437/GSE84437_series_matrix.txt.gz",
    },
}

# GEO scRNA supplementary files
GEO_SCRNA_SUPP = {
    "GSE246662_RAW": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246662/suppl/GSE246662_RAW.tar",
        "dest": "data/raw/scrna/gse246662/GSE246662_RAW.tar",
    },
}

# TCGA from UCSC Xena (RSEM TPM, log2(x+1) normalized)
# These are pan-cancer files - we'll filter by project later
TCGA_XENA = {
    "tcga_rsem_gene_tpm": {
        "url": "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz",
        "dest": "data/raw/bulk/tcga/tcga_RSEM_gene_tpm.gz",
    },
    "tcga_phenotype": {
        "url": "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
        "dest": "data/raw/bulk/tcga/TCGA_phenotype.tsv.gz",
    },
}


def main():
    print("=" * 60)
    print("GC-NKGraph-Atlas Dataset Downloader")
    print("=" * 60)

    results = {}

    # ---- GEO Bulk ----
    print("\n--- GEO Bulk Expression Datasets ---")
    for gse_id, info in GEO_BULK.items():
        dest = info["dest"]
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            size_mb = os.path.getsize(dest) / 1e6
            print(f"\n{gse_id}: Already downloaded ({size_mb:.1f} MB) - SKIP")
            results[gse_id] = True
            continue

        print(f"\n{gse_id}: {info['url']}")
        ok = download_file(info["url"], dest)
        results[gse_id] = ok

    # ---- GEO scRNA Supplementary ----
    print("\n--- GEO scRNA Supplementary Files ---")
    for name, info in GEO_SCRNA_SUPP.items():
        dest = info["dest"]
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            size_mb = os.path.getsize(dest) / 1e6
            print(f"\n{name}: Already downloaded ({size_mb:.1f} MB) - SKIP")
            results[name] = True
            continue

        print(f"\n{name}: {info['url']}")
        print(f"    Note: scRNA raw data - may be large (several GB)")
        ok = download_file(info["url"], dest, retries=2, timeout=1200)
        results[name] = ok

    # ---- TCGA Xena ----
    print("\n--- TCGA from UCSC Xena ---")
    for name, info in TCGA_XENA.items():
        dest = info["dest"]
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            size_mb = os.path.getsize(dest) / 1e6
            print(f"\n{name}: Already downloaded ({size_mb:.1f} MB) - SKIP")
            results[name] = True
            continue

        print(f"\n{name}: {info['url']}")
        ok = download_file(info["url"], dest, retries=3, timeout=600)
        results[name] = ok

    # ---- Summary ----
    print(f"\n{'=' * 60}")
    print("DOWNLOAD SUMMARY:")
    for name, ok in results.items():
        status = "OK" if ok else "FAILED"
        print(f"  {name:<25} {status}")

    all_ok = all(results.values())
    print(f"\nOverall: {'ALL OK' if all_ok else 'SOME FAILURES'}")
    return all_ok


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
