"""
GC-NKGraph-Atlas Master Pipeline Launcher.

Orchestrates the full analysis from raw data to candidate targets.
Each phase checks for existing outputs and skips if already complete.
Supports incremental runs (--from-phase / --to-phase), synthetic data
mode (--synthetic), and dry-run preview (--dry-run).

Usage:
    # Full run with real data
    python src/pipeline.py

    # Test with synthetic data
    python src/pipeline.py --synthetic

    # Dry run (print what would execute)
    python src/pipeline.py --dry-run

    # Run specific phases
    python src/pipeline.py --from-phase preprocess --to-phase model

Phases:
    download     - Download TCGA + GEO raw data
    preprocess   - Standardize bulk expression matrices
    scrna        - scRNA-seq analysis + NK state scoring
    sst_axis     - Compute SST-axis transcriptional scores
    graph        - Build heterogeneous gene graph
    baselines    - Train & evaluate 6 baseline models
    model        - Train GNN-based NK state classifier
    prioritize   - Rank candidate targets with multi-evidence
    all          - Run everything (default)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Fix Unicode output on Windows (prevents GBK encoding errors)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

# Project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from src.common.io_utils import ensure_dir


# =========================================================================
# Phase registry
# =========================================================================

PHASES: Dict[str, Dict] = {
    "synthetic": {
        "label": "Synthetic Data Generation",
        "script": "src/common/synthetic_data.py",
        "outputs": ["data/synthetic/tcga_stad_expression_synthetic.tsv"],
        "depends_on": [],
    },
    "download": {
        "label": "Data Download",
        "script": None,  # multi-script; handled specially
        "outputs": [
            "data/raw/bulk/tcga_stad/",
            "data/raw/bulk/gse62254/",
        ],
        "depends_on": [],
    },
    "preprocess": {
        "label": "Bulk Preprocessing",
        "script": "src/preprocessing/run_bulk_preprocessing.py",
        "outputs": [
            "data/processed/bulk/tcga_stad_expression.tsv",
        ],
        "depends_on": ["download"],
    },
    "scrna": {
        "label": "scRNA-seq Analysis",
        "script": "src/scrna_analysis/run_scrna_pipeline.py",
        "outputs": [
            "data/processed/scrna/gc_integrated.h5ad",
            "results/tables/nk_state_labels.tsv",
        ],
        "depends_on": ["preprocess"],
    },
    "sst_axis": {
        "label": "SST Axis Scoring",
        "script": "src/topology/sst_axis.py",
        "outputs": [
            "results/tables/sst_axis_scores_single_cell.tsv",
            "results/figures/fig9_sst_axis_scores.pdf",
        ],
        "depends_on": ["scrna"],
    },
    "graph": {
        "label": "Graph Construction",
        "script": "src/graph_construction/build_heterograph.py",
        "outputs": [
            "data/processed/graph/nodes.tsv",
            "data/processed/graph/edges.tsv",
        ],
        "depends_on": ["sst_axis"],
    },
    "baselines": {
        "label": "Baseline Models",
        "script": "src/baselines/run_all_baselines.py",
        "outputs": [
            "results/tables/baseline_internal_results.tsv",
        ],
        "depends_on": ["preprocess"],
    },
    "model": {
        "label": "GNN Model Training",
        "script": "src/models/gc_nkgraph_atlas.py",
        "outputs": [
            "results/tables/gc_nkgraph_gnn_internal_results.tsv",
        ],
        "depends_on": ["graph", "preprocess"],
    },
    "prioritize": {
        "label": "Target Prioritization",
        "script": "src/interpretation/prioritize_targets.py",
        "outputs": [
            "results/tables/candidate_evidence_matrix.tsv",
            "results/tables/top_candidate_targets.tsv",
        ],
        "depends_on": ["sst_axis", "scrna"],
    },
}


def log(msg: str, level: str = "INFO") -> None:
    prefix = {"INFO": "  ", "PHASE": "▶ ", "DONE": "✔ ", "SKIP": "○ ", "ERR": "✗ "}
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {prefix.get(level, '  ')}{msg}", flush=True)


def check_outputs(phase_name: str, synthetic: bool = False) -> bool:
    """Check if a phase's expected outputs exist."""
    phase = PHASES.get(phase_name)
    if phase is None:
        return False

    for path_pattern in phase["outputs"]:
        full_path = PROJECT_ROOT / path_pattern
        if "*" in path_pattern or "{" in path_pattern:
            # Glob pattern — check with glob
            import glob
            matches = glob.glob(str(full_path))
            if matches:
                return True
        elif full_path.exists():
            return True

    return False


def resolve_synthetic_path(path: str) -> str:
    """Map real data paths to synthetic equivalents."""
    mapping = {
        "data/processed/bulk/tcga_stad_expression.tsv":
            "data/synthetic/tcga_stad_expression_synthetic.tsv",
        "data/processed/bulk/tcga_stad_clinical.tsv":
            "data/synthetic/tcga_stad_clinical_synthetic.tsv",
        "results/tables/nk_state_labels.tsv":
            "data/synthetic/nk_state_labels_synthetic.tsv",
        "data/processed/scrna/gc_integrated.h5ad":
            "data/synthetic/gc_integrated_synthetic.h5ad",
        "data/processed/scrna/gc_nk_subset.h5ad":
            "data/synthetic/gc_nk_subset_synthetic.h5ad",
        "data/raw/prior_networks":
            "data/synthetic/prior_networks",
    }
    return mapping.get(path, path)


def run_phase(phase_name: str, synthetic: bool = False, force: bool = False,
              extra_args: Optional[List[str]] = None) -> bool:
    """Execute a single phase.

    Returns True if successful (or skipped because already done).
    """
    phase = PHASES.get(phase_name)
    if phase is None:
        log(f"Unknown phase: {phase_name}", "ERR")
        return False

    label = phase["label"]
    script = phase["script"]

    # Check if already done (unless forced)
    if not force and check_outputs(phase_name, synthetic):
        log(f"{label} — already complete, skipping", "SKIP")
        return True

    log(f"{label}...", "PHASE")

    if script is None:
        if phase_name == "download":
            return _run_download_phase(synthetic)
        log(f"No script for phase '{phase_name}'", "ERR")
        return False

    script_path = PROJECT_ROOT / script
    if not script_path.exists():
        log(f"Script not found: {script_path}", "ERR")
        return False

    cmd = [sys.executable, str(script_path)]
    if extra_args:
        cmd.extend(extra_args)

    try:
        result = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=1800,  # 30 min max per phase
        )
        if result.returncode == 0:
            log(f"{label} — done", "DONE")
            # Print last few lines of output
            for line in result.stdout.strip().split("\n")[-3:]:
                if line.strip():
                    log(f"  {line.strip()}")
            return True
        else:
            log(f"{label} — FAILED (exit {result.returncode})", "ERR")
            for line in result.stderr.strip().split("\n")[-5:]:
                log(f"  STDERR: {line.strip()}")
            return False
    except subprocess.TimeoutExpired:
        log(f"{label} — TIMEOUT (>30min)", "ERR")
        return False
    except Exception as e:
        log(f"{label} — ERROR: {e}", "ERR")
        return False


def _run_download_phase(synthetic: bool = False) -> bool:
    """Run data download (multi-script phase)."""
    if synthetic:
        log("  Synthetic mode — skipping real data download", "SKIP")
        return True

    ok = True

    # TCGA
    tcga_script = PROJECT_ROOT / "src/data_download/download_tcga.py"
    if tcga_script.exists():
        try:
            subprocess.run(
                [sys.executable, str(tcga_script), "--project", "ALL"],
                cwd=str(PROJECT_ROOT), timeout=600,
            )
        except Exception as e:
            log(f"  TCGA download: {e}", "ERR")
            ok = False

    # GEO
    geo_script = PROJECT_ROOT / "src/data_download/download_geo.py"
    if geo_script.exists():
        try:
            subprocess.run(
                [sys.executable, str(geo_script)],
                cwd=str(PROJECT_ROOT), timeout=600,
            )
        except Exception as e:
            log(f"  GEO download: {e}", "ERR")
            ok = False

    return ok


def _generate_synthetic_data() -> bool:
    """Generate the full synthetic dataset."""
    script = PROJECT_ROOT / "src/common/synthetic_data.py"
    try:
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode == 0:
            log("Synthetic data generation — done", "DONE")
            return True
        else:
            log(f"Synthetic data generation FAILED: {result.stderr[-200:]}", "ERR")
            return False
    except Exception as e:
        log(f"Synthetic data error: {e}", "ERR")
        return False


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="GC-NKGraph-Atlas Master Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python src/pipeline.py                        # Full run (real data)
  python src/pipeline.py --synthetic             # Full run (synthetic data)
  python src/pipeline.py --dry-run               # Show what would run
  python src/pipeline.py --from-phase graph      # Start from graph construction
  python src/pipeline.py --to-phase model        # Stop after model training
        """,
    )
    parser.add_argument("--synthetic", action="store_true",
                       help="Use synthetic data for end-to-end testing")
    parser.add_argument("--dry-run", action="store_true",
                       help="Print phases without executing")
    parser.add_argument("--from-phase", default="download",
                       choices=list(PHASES.keys()),
                       help="Phase to start from (default: download)")
    parser.add_argument("--to-phase", default="prioritize",
                       choices=list(PHASES.keys()),
                       help="Phase to stop at (default: prioritize)")
    parser.add_argument("--force", action="store_true",
                       help="Re-run phases even if outputs exist")
    parser.add_argument("--list-phases", action="store_true",
                       help="List all phases and exit")
    args = parser.parse_args()

    if args.list_phases:
        print("\nPipeline phases:")
        for name, info in PHASES.items():
            deps = " → ".join(info["depends_on"]) or "none"
            print(f"  {name:<14}  depends: {deps}")
            print(f"  {'':14}  {info['label']}")
        return

    # Phase order
    ordered_phases = [
        "download", "preprocess", "scrna", "sst_axis",
        "graph", "baselines", "model", "prioritize",
    ]

    # Determine range
    try:
        start_idx = ordered_phases.index(args.from_phase)
    except ValueError:
        log(f"Unknown --from-phase: {args.from_phase}", "ERR")
        sys.exit(1)
    try:
        end_idx = ordered_phases.index(args.to_phase) + 1
    except ValueError:
        log(f"Unknown --to-phase: {args.to_phase}", "ERR")
        sys.exit(1)

    run_phases = ordered_phases[start_idx:end_idx]

    # Header
    print()
    log("=" * 60)
    log("GC-NKGraph-Atlas PIPELINE")
    log(f"  Mode: {'SYNTHETIC (test)' if args.synthetic else 'REAL DATA'}")
    log(f"  Phases: {' → '.join(run_phases)}")
    log(f"  Force: {args.force}")
    log("=" * 60)
    print()

    # Dry run
    if args.dry_run:
        log("DRY RUN — no execution", "INFO")
        for p in run_phases:
            info = PHASES[p]
            done = check_outputs(p, args.synthetic) and not args.force
            status = "(already done — would skip)" if done else "(would run)"
            log(f"  {p:<14} {info['label']:<30} {status}")
        print()
        return

    # Setup
    start_time = datetime.now()

    # Generate synthetic data first if needed
    if args.synthetic:
        if not _generate_synthetic_data():
            log("Synthetic data generation failed — aborting", "ERR")
            sys.exit(1)

    # Run each phase
    results: Dict[str, bool] = {}
    for p in run_phases:
        if args.force:
            log(f"Forcing re-run of phase: {p}", "INFO")

        ok = run_phase(p, synthetic=args.synthetic, force=args.force)
        results[p] = ok
        if not ok:
            log(f"Phase '{p}' failed — stopping pipeline", "ERR")
            break

    # Summary
    elapsed = (datetime.now() - start_time).total_seconds()
    print()
    log("=" * 60)
    log("PIPELINE COMPLETE")
    log(f"  Elapsed: {elapsed:.0f}s")
    for p, ok in results.items():
        status = "PASS" if ok else "FAIL"
        level = "DONE" if ok else "ERR"
        log(f"  {p:<14} {status}", level)
    log("=" * 60)
    print()

    # Exit code
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
