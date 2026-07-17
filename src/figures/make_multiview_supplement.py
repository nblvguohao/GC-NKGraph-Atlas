"""Create the supplementary multiview performance and stability figure."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DISPLAY_VARIANTS = ["no_graph", "merged_svd", "uniform_multiview", "learned_multiview"]


def make_multiview_supplement(
    summary: pd.DataFrame,
    weights: pd.DataFrame,
    output_prefix: Path,
) -> Tuple[Path, Path]:
    """Write a two-panel PNG/PDF with seed stability shown explicitly."""

    primary = summary[
        (summary["feature_mode"] == "masked") & summary["variant"].isin(DISPLAY_VARIANTS)
    ].copy()
    primary["variant"] = pd.Categorical(primary["variant"], DISPLAY_VARIANTS, ordered=True)
    primary = primary.sort_values("variant")
    learned_weights = weights[
        (weights["feature_mode"] == "masked") & (weights["variant"] == "learned_multiview")
    ].copy()
    if primary.empty or learned_weights.empty:
        raise ValueError("Masked summary and learned_multiview weights are required")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    positions = np.arange(len(primary))
    axes[0].errorbar(
        primary["AUROC_mean"],
        positions - 0.08,
        xerr=primary["AUROC_std"],
        fmt="o",
        capsize=3,
        label="AUROC",
        color="#2563eb",
    )
    axes[0].errorbar(
        primary["AUPRC_mean"],
        positions + 0.08,
        xerr=primary["AUPRC_std"],
        fmt="s",
        capsize=3,
        label="AUPRC",
        color="#dc2626",
    )
    axes[0].set_yticks(positions, primary["variant"].astype(str))
    axes[0].invert_yaxis()
    axes[0].set_xlabel("External LIHC performance (mean +/- seed SD)")
    axes[0].set_title("A  Label-masked external benchmark")
    axes[0].grid(axis="x", alpha=0.25)
    axes[0].legend(frameon=False)

    view_order = list(dict.fromkeys(learned_weights["view"]))
    values = [learned_weights.loc[learned_weights["view"] == view, "weight"] for view in view_order]
    axes[1].boxplot(values, tick_labels=view_order, showmeans=True)
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].set_ylabel("Learned global softmax weight")
    axes[1].set_title("B  View-weight stability across seeds")
    axes[1].grid(axis="y", alpha=0.25)
    fig.suptitle("Lightweight multiview fusion: external performance and weight stability")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    paths = make_multiview_supplement(
        pd.read_csv(args.summary, sep="\t"),
        pd.read_csv(args.weights, sep="\t"),
        args.output_prefix,
    )
    print("\n".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
