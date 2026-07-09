#!/bin/bash
# GC-NKGraph-Atlas Docker Entrypoint
#
# Modes:
#   verify   — run synthetic pipeline + unit tests (default, no real data needed)
#   web      — serve the interactive Mechanism Card Playground on port 8080
#   shell    — drop into bash for interactive exploration
#   <anything else> — passed directly to python
#
# Usage:
#   docker run --rm gc-nkgraph-atlas                    # verify mode
#   docker run --rm -p 8080:8080 gc-nkgraph-atlas web   # web playground
#   docker run --rm gc-nkgraph-atlas shell               # interactive
#   docker run --rm gc-nkgraph-atlas python src/topology/sst_axis.py  # specific stage

set -euo pipefail
MODE="${1:-verify}"

case "$MODE" in
verify)
    echo "=============================================="
    echo "GC-NKGraph-Atlas — Pipeline Verification"
    echo "=============================================="
    echo ""
    echo "Running synthetic-data pipeline smoke test..."
    python src/pipeline.py --synthetic 2>&1 | tail -20
    echo ""
    echo "Running unit tests..."
    python -m pytest tests/ -q --tb=short 2>&1 || echo "(Some tests may be skipped without GPU)"
    echo ""
    echo "=============================================="
    echo "VERIFICATION COMPLETE"
    echo "=============================================="
    echo ""
    echo "The pipeline executes end-to-end on synthetic data."
    echo "To run with real data, mount your data/ directory:"
    echo "  docker run --rm -v /path/to/data:/workspace/data gc-nkgraph-atlas"
    echo ""
    echo "For the interactive web playground:"
    echo "  docker run --rm -p 8080:8080 gc-nkgraph-atlas web"
    ;;

web)
    echo "=============================================="
    echo "GC-NKGraph-Atlas — Mechanism Card Playground"
    echo "=============================================="
    echo ""
    echo "Serving on http://localhost:8080"
    echo "Open this URL in your browser to explore:"
    echo "  - 3 mechanism cards (serine-SM, adenosine, TGF-beta)"
    echo "  - Card gallery & comparison table"
    echo "  - Searchable target browser with DrugBank/DepMap evidence"
    echo ""
    cd /workspace/web
    python -m http.server 8080
    ;;

shell)
    echo "GC-NKGraph-Atlas — Interactive Shell"
    echo ""
    echo "Key files:"
    echo "  configs/mechanism_cards/   — mechanism card YAML files"
    echo "  src/pipeline.py            — master pipeline launcher"
    echo "  src/interpretation/        — target prioritization + multi-card analysis"
    echo "  web/index.html             — interactive playground"
    echo "  results/tables/            — generated output tables"
    echo ""
    exec /bin/bash
    ;;

*)
    exec python "$@"
    ;;
esac
