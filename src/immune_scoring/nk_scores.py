"""
GC-NKGraph-Atlas NK immune-state scoring pipeline.

Computes NK infiltration, cytotoxicity, dysfunction, and exclusion scores
from bulk transcriptomes. Generates four immune-state labels.

Usage:
    python src/immune_scoring/nk_scores.py --config configs/data_config.yaml
"""

import os
import sys
import json
import argparse
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.logging import Logger
from src.common.io_utils import load_config, save_table, ensure_dir, load_table


# ---- Core NK gene signatures (from design doc Section 5) ----
NK_MARKERS = [
    "NCAM1", "NCR1", "KLRD1", "KLRK1", "KLRF1",
    "NKG7", "GNLY", "GZMB", "PRF1", "FCGR3A",
    "XCL1", "XCL2", "CCL5", "IFNG", "TYROBP",
]

NK_CYTOTOXICITY_GENES = [
    "NKG7", "GNLY", "GZMB", "PRF1", "IFNG", "XCL1", "XCL2", "CCL5",
]

NK_DYSFUNCTION_GENES = [
    "KLRC1", "TIGIT", "CD96", "HAVCR2", "TOX",
    "ENTPD1", "TGFB1", "HIF1A", "NT5E", "ADORA2A",
]

CAF_ECM_TGFB_GENES = [
    "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN",
    "ACTA2", "FAP", "TAGLN", "TGFBI", "TGFB1", "CXCL12",
]


def mean_zscore(df: pd.DataFrame, genes: list) -> pd.Series:
    """Compute the mean of z-score normalized values for a gene set.

    For each gene in the set (if present in the dataframe), z-score normalize
    across samples, then average.
    """
    available = [g for g in genes if g in df.columns]
    if not available:
        return pd.Series(0.0, index=df.index)

    # Z-score per gene across samples
    subset = df[available]
    z = (subset - subset.mean(0)) / subset.std(0, ddof=0)
    # Replace NaN (from zero-variance genes) with 0
    z = z.fillna(0)
    return z.mean(axis=1)


def compute_nk_scores(
    expression_df: pd.DataFrame,
    logger: Logger,
    dataset_name: str = "",
) -> pd.DataFrame:
    """Compute NK immune scores for a bulk expression matrix.

    Args:
        expression_df: Samples x Genes DataFrame.
        logger: Logger instance.
        dataset_name: Name for logging.

    Returns:
        DataFrame with NK score columns.
    """
    scores = pd.DataFrame(index=expression_df.index)

    scores["NK_infiltration_score"] = mean_zscore(expression_df, NK_MARKERS)
    scores["NK_cytotoxicity_score"] = mean_zscore(expression_df, NK_CYTOTOXICITY_GENES)
    scores["NK_dysfunction_score"] = (
        mean_zscore(expression_df, NK_DYSFUNCTION_GENES)
        - mean_zscore(expression_df, NK_CYTOTOXICITY_GENES)
    )
    scores["NK_exclusion_score"] = (
        mean_zscore(expression_df, CAF_ECM_TGFB_GENES)
        - scores["NK_infiltration_score"]
    )

    n_genes = {
        "NK_MARKERS": sum(1 for g in NK_MARKERS if g in expression_df.columns),
        "NK_CYTOTOXICITY": sum(1 for g in NK_CYTOTOXICITY_GENES if g in expression_df.columns),
        "NK_DYSFUNCTION": sum(1 for g in NK_DYSFUNCTION_GENES if g in expression_df.columns),
        "CAF_ECM_TGFB": sum(1 for g in CAF_ECM_TGFB_GENES if g in expression_df.columns),
    }

    logger.ok(
        phase="NK_SCORING",
        message=f"{dataset_name}: NK scores computed. Genes found: {n_genes}",
        script=__file__,
    )

    return scores


def assign_immune_states(
    scores_df: pd.DataFrame,
    thresholds: dict = None,
    logger: Logger = None,
) -> pd.DataFrame:
    """Assign four immune states based on NK score thresholds.

    States:
        NK-hot-cytotoxic:    high infiltration, high cytotoxicity, low dysf.
        NK-hot-dysfunctional: high infiltration, low cytotoxicity, high dysf.
        NK-cold/excluded:    low infiltration, high exclusion
        NK-intermediate:     remaining

    Args:
        scores_df: DataFrame with NK score columns.
        thresholds: dict with quantile thresholds. If None, computed from data.
        logger: Logger instance.

    Returns:
        DataFrame with added 'nk_immune_state' column.
    """
    df = scores_df.copy()

    if thresholds is None:
        thresholds = {
            "infiltration_high": df["NK_infiltration_score"].median(),
            "cytotoxicity_high": df["NK_cytotoxicity_score"].median(),
            "dysfunction_high": df["NK_dysfunction_score"].median(),
            "exclusion_high": df["NK_exclusion_score"].median(),
        }
        if logger:
            logger.ok(
                phase="NK_SCORING",
                message=f"Thresholds computed from data: {thresholds}",
                script=__file__,
            )

    # Assign states
    conditions = [
        (df["NK_infiltration_score"] >= thresholds["infiltration_high"])
        & (df["NK_cytotoxicity_score"] >= thresholds["cytotoxicity_high"])
        & (df["NK_dysfunction_score"] <= thresholds["dysfunction_high"]),
        (df["NK_infiltration_score"] >= thresholds["infiltration_high"])
        & (df["NK_cytotoxicity_score"] <= thresholds["cytotoxicity_high"])
        & (df["NK_dysfunction_score"] >= thresholds["dysfunction_high"]),
        (df["NK_infiltration_score"] < thresholds["infiltration_high"])
        & (df["NK_exclusion_score"] >= thresholds["exclusion_high"]),
    ]
    choices = ["NK-hot-cytotoxic", "NK-hot-dysfunctional", "NK-cold/excluded"]
    df["nk_immune_state"] = np.select(conditions, choices, default="NK-intermediate")

    if logger:
        counts = df["nk_immune_state"].value_counts().to_dict()
        logger.ok(
            phase="NK_SCORING",
            message=f"State distribution: {counts}",
            script=__file__,
        )

    return df, thresholds


def main():
    parser = argparse.ArgumentParser(description="Compute NK immune scores")
    parser.add_argument("--config", default="configs/data_config.yaml", help="Data config")
    parser.add_argument("--output-dir", default="results/tables", help="Output directory")
    args = parser.parse_args()

    logger = Logger()
    config = load_config(args.config)
    output_dir = ensure_dir(args.output_dir)

    # Process training and external validation datasets
    all_datasets = (
        config.get("bulk_datasets", [])
        + config.get("positive_control_bulk_datasets", [])
    )

    all_scores = []
    all_labels = []
    global_thresholds = None

    for ds in all_datasets:
        expr_path = ds.get("expression_path", "")
        name = ds["name"]
        role = ds.get("role", "unknown")

        if not os.path.exists(expr_path):
            logger.skip(
                phase="NK_SCORING",
                message=f"{name}: expression not found at {expr_path}",
                script=__file__,
            )
            continue

        # Load expression
        expr = load_table(expr_path)
        # Ensure samples are rows, genes columns
        if expr.index.name == "gene" or "gene" in str(expr.index.name):
            expr = expr.T

        # Compute scores
        scores = compute_nk_scores(expr, logger, name)
        scores["dataset"] = name
        scores["role"] = role
        all_scores.append(scores)

        # For training primary, compute thresholds
        if role == "train_primary":
            labeled, global_thresholds = assign_immune_states(scores, logger=logger)
            save_table(
                labeled,
                os.path.join(args.output_dir, "nk_state_labels.tsv"),
                provenance={"dataset": name, "threshold_source": "train_primary"},
            )
            # Save thresholds for external validation
            save_table(
                pd.DataFrame([global_thresholds]),
                os.path.join(output_dir, "nk_state_thresholds.json"),
            )
            with open(os.path.join(output_dir, "nk_state_thresholds.json"), "w") as f:
                json.dump(global_thresholds, f, indent=2)
            all_labels.append(labeled)

        elif global_thresholds is not None:
            labeled, _ = assign_immune_states(scores, thresholds=global_thresholds, logger=logger)
            labeled["dataset"] = name
            all_labels.append(labeled)

    # Merge and save scores
    if all_scores:
        merged_scores = pd.concat(all_scores)
        save_table(
            merged_scores,
            os.path.join(args.output_dir, "nk_scores_bulk.tsv"),
            provenance={"script": __file__},
        )
        logger.ok(
            phase="NK_SCORING",
            message=f"All NK scores saved. Shape: {merged_scores.shape}",
            script=__file__,
        )

    if all_labels:
        merged_labels = pd.concat(all_labels)
        save_table(
            merged_labels,
            os.path.join(args.output_dir, "nk_state_labels.tsv"),
            provenance={"script": __file__},
        )

    # Write label definition document
    label_def = """# NK Immune-State Label Definition

## Scores
- NK_infiltration_score: mean z-score of {nk_markers}
- NK_cytotoxicity_score: mean z-score of {cytotoxicity_genes}
- NK_dysfunction_score: mean z-score(dysfunction) - mean z-score(cytotoxicity)
- NK_exclusion_score: mean z-score(CAF/ECM) - NK_infiltration_score

## Immune States
- NK-hot-cytotoxic: high infiltration, high cytotoxicity, low dysfunction
- NK-hot-dysfunctional: high infiltration, low cytotoxicity, high dysfunction
- NK-cold/excluded: low infiltration, high exclusion
- NK-intermediate: remaining samples

## Thresholds
Computed on training primary dataset (TCGA-STAD) as median values.
External validation uses saved thresholds from training.

## Caveats
- Gene signatures are literature-derived starting sets.
- scRNA-derived cell-type attribution may refine scores.
- These are transcriptional proxies, not functional measurements.
""".format(
        nk_markers=", ".join(NK_MARKERS),
        cytotoxicity_genes=", ".join(NK_CYTOTOXICITY_GENES),
    )

    with open(os.path.join(args.output_dir, "label_definition.md"), "w") as f:
        f.write(label_def)

    logger.ok(
        phase="NK_SCORING",
        message="Label definition written to results/tables/label_definition.md",
        script=__file__,
    )


if __name__ == "__main__":
    main()
