"""Render the categorical real-data recoverability atlas as Figure S2."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
table = ROOT / "submission_bundle_BiB/03_supplementary/tables/recoverability_atlas.tsv"
out = ROOT / "submission_bundle_BiB/02_figures"
df = pd.read_csv(table, sep="\t")
fig, ax = plt.subplots(figsize=(9, 3.5))
colors = df.status.map({"recovered": "#287d8e", "not_recovered": "#b65b50"}).fillna("#9aa0a6")
ax.scatter(df.layer, df.card_id, s=180, c=colors, marker="s")
ax.set_xlabel("Pre-specified biological layer")
ax.set_ylabel("Mechanism card")
ax.set_title("Figure S2. Real-data recoverability atlas (bulk transcriptome)")
ax.legend(handles=[
    Line2D([0], [0], marker="s", color="w", label="recovered", markerfacecolor="#287d8e", markersize=11),
    Line2D([0], [0], marker="s", color="w", label="not recovered", markerfacecolor="#b65b50", markersize=11),
], loc="upper left", frameon=False)
fig.tight_layout()
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "figS2_recoverability_atlas.png", dpi=220)
fig.savefig(out / "figS2_recoverability_atlas.pdf")
