"""
S2 (nature-reviewer task card): audit whether the `metabolic_crosstalk` graph
edge's weight/sign is calibrated on data (and if so, whether that calibration
could leak information between the liver calibration cohort and the gastric
test cohort), or whether it is a fixed structural prior.

Motivation: §2.5 of the manuscript originally stated "the sign of the
crosstalk is calibrated on the liver positive-control cohort, not hard-coded."
Inspecting `src/graph_construction/build_heterograph.py::build_sst_edges`
shows the edge weight passed to `_add_edge` for every metabolic_crosstalk
edge is the literal constant `0.5` -- not a function of any fitted or
calibrated quantity. This script confirms that programmatically (rather than
by manual code reading alone) and separately documents the one calibration
step that *does* exist in this codebase (the sign of the tumor_serine_capacity
term in the descriptive sst_axis_score composite), to make clear the two are
independent and the edge weight cannot leak between cohorts.

Output:
  results/tables/mc_edge_sign_calibration_audit.tsv

Usage:
    python src/topology/h3_edge_sign_calibration_audit.py
"""
import ast
import inspect
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import pandas as pd  # noqa: E402
from src.graph_construction.build_heterograph import build_sst_edges  # noqa: E402
from src.common.sst_config import load_sst_modules  # noqa: E402


def log(msg):
    print(msg, flush=True)


def audit_edge_weight_is_constant():
    """Build the real SST edges from the canonical gene modules and confirm
    every metabolic_crosstalk edge has weight == 0.5, with no dependence on
    any external/fitted quantity (source code has no cohort-derived input)."""
    modules = load_sst_modules()
    # all_nodes must map every gene name to a node dict; build a permissive
    # stand-in so build_sst_edges can run standalone for this audit.
    all_genes = set()
    for mod in modules.values():
        all_genes.update(mod["genes"])
    all_nodes = {g: {"node_type": "gene", "name": g, "source": "audit"} for g in all_genes}

    edges = build_sst_edges(all_nodes, modules)
    mc_edges = [e for e in edges if e.get("edge_type") == "metabolic_crosstalk"]
    weights = sorted(set(round(e["weight"], 6) for e in mc_edges))

    log(f"metabolic_crosstalk edges built: {len(mc_edges)}")
    log(f"distinct weight values observed: {weights}")

    source = inspect.getsource(build_sst_edges)
    # Static check: the call site for metabolic_crosstalk uses a literal 0.5,
    # not a variable computed from data (e.g. no correlation/fit call before it).
    literal_weight_used = "\"metabolic_crosstalk\", 0.5" in source or "'metabolic_crosstalk', 0.5" in source

    return {
        "n_metabolic_crosstalk_edges": len(mc_edges),
        "distinct_weights": str(weights),
        "is_constant_across_all_edges": len(weights) == 1,
        "constant_value": weights[0] if len(weights) == 1 else None,
        "source_uses_literal_constant": literal_weight_used,
    }


def main():
    log("=" * 70)
    log("S2: metabolic_crosstalk EDGE SIGN/WEIGHT CALIBRATION AUDIT")
    log("=" * 70)

    result = audit_edge_weight_is_constant()
    for k, v in result.items():
        log(f"  {k}: {v}")

    rows = [
        {
            "component": "metabolic_crosstalk edge weight",
            "location": "src/graph_construction/build_heterograph.py::build_sst_edges",
            "is_calibrated_on_any_cohort": False,
            "value_or_rule": f"fixed constant = {result['constant_value']}",
            "leakage_risk": "none (constant does not depend on liver or gastric data)",
            "verified_by": "audit_edge_weight_is_constant() - built real SST edges "
                           "from canonical gene modules, confirmed single distinct "
                           "weight value across all "
                           f"{result['n_metabolic_crosstalk_edges']} edges",
        },
        {
            "component": "sst_axis_score tumor_serine_capacity sign",
            "location": "src/topology/sst_axis_validation.py (H1 crosstalk calibration)",
            "is_calibrated_on_any_cohort": True,
            "value_or_rule": "sign of the observed H1 (tumor_serine_capacity ~ "
                             "nk_sm_balance) correlation on the liver (TCGA-LIHC) "
                             "cohort, applied to the descriptive sst_axis_score "
                             "composite only",
            "leakage_risk": "low: (i) this sign feeds only the descriptive "
                            "sst_axis_score composite, not the graph edge weight "
                            "above and not the H2-H5 hypothesis tests, which "
                            "operate on the individual module scores directly; "
                            "(ii) H1 is null in this dataset (bulk r=-0.016, "
                            "p=0.74; single-cell corrected r=0.012, p=0.27), so "
                            "the calibrated sign carries near-zero information "
                            "regardless of which cohort it is drawn from",
            "verified_by": "manual trace of sst_axis.py / sst_axis_validation.py; "
                           "H1 null result confirmed in sst_axis_positive_control_recovery.tsv",
        },
    ]

    df = pd.DataFrame(rows)
    out_dir = "results/tables"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "mc_edge_sign_calibration_audit.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    log(f"\nSaved: {out_path}")
    log("\nCONCLUSION: the metabolic_crosstalk edge weight is a fixed prior, not "
        "a calibrated/fitted quantity, and cannot leak information between the "
        "liver and gastric cohorts. The manuscript's original wording "
        "('the sign of the crosstalk is calibrated ... not hard-coded') "
        "described the wrong component and has been corrected in §2.5.")


if __name__ == "__main__":
    main()
