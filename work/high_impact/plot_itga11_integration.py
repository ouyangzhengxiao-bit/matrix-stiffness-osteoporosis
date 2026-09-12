#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/研究升级/high_impact_v4"

ctx = pd.read_csv(OUT / "study_by_cell_state_effects.csv")
mat = ctx.pivot(index="cluster", columns="study", values="mean_rank_effect")
order = ["Cl_0", "Cl_1", "Cl_2", "Cl_3", "Cl_4"]
mat = mat.loc[order]
sig = pd.read_csv(OUT / "stiffness_signature_in_itga11_bulk.csv")
labels = ["ITGA11 high vs low\n(day 0)", "ITGA11 high vs low\n(after expansion)", "ITGA11 knockdown\nvs control"]

sns.set_theme(style="white", font_scale=0.9)
fig = plt.figure(figsize=(11.2, 4.8), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1])
ax1 = fig.add_subplot(gs[0, 0])
sns.heatmap(mat, cmap="RdBu_r", center=0, vmin=-0.65, vmax=0.65, annot=True, fmt=".2f",
            linewidths=.5, cbar_kws={"label": "Mean within-study rank effect"}, ax=ax1)
ax1.set_title("A  Cell-state programs respond differently across stiffness studies", loc="left", fontweight="bold")
ax1.set_xlabel("")
ax1.set_ylabel("GSE317531 cell state")
ax1.tick_params(axis="x", rotation=45)

ax2 = fig.add_subplot(gs[0, 1])
x = np.arange(len(sig))
ax2.bar(x, sig.directional_score, color=["#4863A0", "#6F8FC6", "#8A8A8A"], width=.65)
ax2.errorbar(x, sig.directional_score, yerr=1.96*sig.null_sd, fmt="none", ecolor="black", capsize=4, lw=1)
ax2.axhline(0, color="black", lw=.8)
for i, r in sig.iterrows():
    ax2.text(i, r.directional_score + (-.018 if r.directional_score < 0 else .018),
             f"q={r.BH_q:.2f}", ha="center", va="top" if r.directional_score < 0 else "bottom", fontsize=9)
ax2.set_xticks(x, labels, rotation=20, ha="right")
ax2.set_ylabel("Stiffness directional score\n(mean positive − mean negative genes)")
ax2.set_title("B  No alignment with the ITGA11 osteogenic state", loc="left", fontweight="bold")
ax2.spines[["top", "right"]].set_visible(False)
fig.suptitle("Public-data triangulation identifies context dependence rather than a universal osteogenic stiffness program",
             fontsize=12, fontweight="bold")
for ext, dpi in [("png", 220), ("tif", 300)]:
    fig.savefig(OUT / f"Figure_v4_cell_state_context.{ext}", dpi=dpi, bbox_inches="tight")
plt.close(fig)
