"""
GC-NKGraph-Atlas TCGA data downloader.

Downloads TCGA-STAD (gastric) and TCGA-LIHC (liver) bulk expression + clinical data.

Requires: gdc-client or TCGAbiolinks (R) or cggh/cloudman (Python).
Fallback: Manual download from GDC Data Portal instructions.
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Optional

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import load_config, ensure_dir


# --- Configuration ---
TCGA_PROJECTS = {
    "TCGA-STAD": {
        "name": "Stomach Adenocarcinoma",
        "cancer_type": "gastric_cancer",
        "role": "train_primary",
    },
    "TCGA-LIHC": {
        "name": "Liver Hepatocellular Carcinoma",
        "cancer_type": "liver_hcc",
        "role": "positive_control_liver",
    },
}

# Default output paths (relative to project root)
DEFAULT_OUTPUT_DIR = "data/raw/bulk"


def check_gdc_client() -> bool:
    """Check if gdc-client is available."""
    import shutil
    return shutil.which("gdc-client") is not None


def check_gdc_manifest(project_id: str, manifest_dir: str) -> Optional[str]:
    """Check if a manifest file exists for the given project."""
    manifest_path = os.path.join(manifest_dir, f"{project_id}_manifest.txt")
    if os.path.exists(manifest_path):
        return manifest_path
    return None


def download_tcga_rnaseq(
    project_id: str,
    output_dir: str,
    logger: Logger,
    manifest_file: Optional[str] = None,
    use_tcgabiolinks: bool = False,
):
    """Download TCGA RNA-seq expression data.

    Args:
        project_id: TCGA project ID (e.g., TCGA-STAD, TCGA-LIHC)
        output_dir: Directory to save downloaded files
        logger: Logger instance
        manifest_file: Optional GDC manifest file
        use_tcgabiolinks: Use R TCGAbiolinks if True, else gdc-client
    """
    project_dir = os.path.join(output_dir, project_id.lower().replace("-", "_"))
    ensure_dir(project_dir)

    if use_tcgabiolinks:
        logger.needs_review(
            phase="DATA_DOWNLOAD",
            message=f"TCGAbiolinks (R) download for {project_id} not yet implemented. "
                    f"To use: run R script with TCGAbiolinks::GDCdownload() and "
                    f"TCGAbiolinks::GDCprepare(). Output dir: {project_dir}",
            script=__file__,
        )
        return

    if check_gdc_client():
        client_name = "gdc-client"
        cmd = f"{client_name} download -m {manifest_file} -d {project_dir}"
        logger.ok(
            phase="DATA_DOWNLOAD",
            message=f"Run gdc-client manually:\n  {cmd}",
            script=__file__,
        )
    else:
        logger.needs_review(
            phase="DATA_DOWNLOAD",
            message=f"gdc-client not found. Download TCGA-{project_id} data from "
                    f"https://portal.gdc.cancer.gov/ and place files in {project_dir}/.\n\n"
                    f"Required files:\n"
                    f"  1. Gene expression quantifications (RNA-seq, STAR - Counts)\n"
                    f"  2. Clinical data (XML or JSON)\n"
                    f"  3. Sample annotations (e.g., tissue type, survival)",
            script=__file__,
        )


def main():
    parser = argparse.ArgumentParser(description="Download TCGA data for GC-NKGraph-Atlas")
    parser.add_argument("--config", default="configs/data_config.yaml", help="Data config path")
    parser.add_argument("--project", choices=list(TCGA_PROJECTS.keys()) + ["ALL"], default="ALL")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--use-r", action="store_true", help="Use R TCGAbiolinks instead of gdc-client")
    args = parser.parse_args()

    logger = Logger()
    output_dir = ensure_dir(args.output_dir)

    projects = [args.project] if args.project != "ALL" else list(TCGA_PROJECTS.keys())

    for pid in projects:
        logger.ok(
            phase="DATA_DOWNLOAD",
            message=f"Starting download: {pid} ({TCGA_PROJECTS[pid]['name']})\n"
                    f"  Role: {TCGA_PROJECTS[pid]['role']}\n"
                    f"  Output: {output_dir}/{pid.lower().replace('-', '_')}/",
            script=__file__,
        )
        download_tcga_rnaseq(
            project_id=pid,
            output_dir=str(output_dir),
            logger=logger,
            use_tcgabiolinks=args.use_r,
        )

    print(f"\nDone. Instructions logged to {logger.log_path}")
    print(f"After download, run: python src/preprocessing/run_bulk_preprocessing.py")


if __name__ == "__main__":
    main()
