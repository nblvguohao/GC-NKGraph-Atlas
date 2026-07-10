#!/bin/bash
# Master launcher for the five P1 science recompute tasks.
# Run on A100:  ssh user@100.112.165.109  then:
#   cd /data/lgh/GC-NKGraph-Atlas
#   conda activate gc-nkgraph
#   bash src/a100_recompute/run_all.sh
set -euo pipefail

echo "============================================"
echo "GC-NKGraph-Atlas P1 Science Recompute"
echo "Started: $(date)"
echo "============================================"

cd /data/lgh/GC-NKGraph-Atlas

# T5 — H1 result
echo ""
echo ">>> T5: H1 result <<<"
python src/a100_recompute/run_h1_result.py && echo "T5 OK" || echo "T5 FAILED (non-fatal)"

# T15 — Effect sizes + multiple-testing (runs before T6 since T6 reads its output)
echo ""
echo ">>> T15: Effect-size reframe + BH correction <<<"
python src/a100_recompute/run_effect_size_reframe.py && echo "T15 OK" || echo "T15 FAILED (non-fatal)"

# T14 — H3 activation control
echo ""
echo ">>> T14: H3 activation control <<<"
python src/a100_recompute/run_h3_activation_control.py && echo "T14 OK" || echo "T14 FAILED (non-fatal)"

# T6 — De-circ audit
echo ""
echo ">>> T6: De-circ audit <<<"
python src/a100_recompute/run_decirc_audit.py && echo "T6 OK" || echo "T6 FAILED (non-fatal)"

# T7 — Ablation (most expensive — last)
echo ""
echo ">>> T7: Ablation (WARNING: ~10-30 min with GNN training) <<<"
python src/a100_recompute/run_ablation.py && echo "T7 OK" || echo "T7 FAILED (non-fatal)"

echo ""
echo "============================================"
echo "ALL P1 TASKS COMPLETE"
echo "Finished: $(date)"
echo "============================================"
echo ""
echo "Output files on A100:"
echo "  results/tables/sst_axis_positive_control_recovery.tsv       (T5 — includes H1)"
echo "  results/tables/sst_axis_positive_control_recovery_v2.tsv    (T15 — CIs + effect sizes + BH)"
echo "  results/tables/h3_activation_control.tsv                    (T14)"
echo "  results/tables/tumor_intrinsic_pool_audit.tsv               (T6)"
echo "  results/tables/ablation_results.tsv                         (T7)"
echo ""
echo "After scp'ing these back to local, re-run:"
echo "  cd G:/cc/GC-NKGraph-Atlas && python src/figures/make_figures.py"
echo "  to regenerate Figures 1-4 with updated numbers"
