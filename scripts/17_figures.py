#!/usr/bin/env python
"""Figures: baseline (Pre) responder vs non-responder comparisons and a UMAP colored by cell type.

One point per patient (never per cell). Points are colored by therapy. The p-values printed on the plots are
UNCORRECTED two-sided Mann-Whitney U tests, repeated here for orientation only; the corrected results are in
results/*.csv. Inputs come from scripts 08, 09 and 12; the TCF7 values follow script 14."""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

NON_T = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
MIN_CELLS = 20            # same minimum per patient as the TCF7 analysis (script 14)
THERAPY_COLORS = {"anti-PD1": "#1f77b4", "anti-CTLA4+PD1": "#ff7f0e", "anti-CTLA4": "#2ca02c"}
sc.settings.figdir = "docs/figures"

a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")

# ---- patient-level data at baseline ----
scores = pd.read_csv("results/patient_scores.csv")                         # script 12
scores = scores[scores["timepoint"] == "Pre"]

fr = pd.read_csv("results/biopsy_celltype_fractions.csv")                  # script 09
frp = fr.groupby(["patient", "timepoint", "response", "therapy"])[["B"]].mean().reset_index()
frp = frp[frp["timepoint"] == "Pre"]

x = a.raw[:, "TCF7"].X
x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
o = a.obs[["patient", "timepoint", "response", "therapy", "enrichment", "celltype"]].astype(str).copy()
o["tcf7"] = x
o = o[(o["timepoint"] == "Pre") & (o["enrichment"] == "CD45+") & (~o["celltype"].isin(NON_T))]
tcf7 = o.groupby(["patient", "response", "therapy"]).agg(n=("tcf7", "size"), mean=("tcf7", "mean")).reset_index()
tcf7 = tcf7[tcf7["n"] >= MIN_CELLS]

panels = [
    (scores, "score_cytotoxic", "Cytotoxic gene score (T cells)"),
    (scores, "score_memory", "Memory gene score (T cells)"),
    (tcf7, "mean", "TCF7 (mean log2(TPM+1), T cells)"),
    (frp, "B", "B-cell fraction of CD45+ cells"),
]

rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, 4, figsize=(17, 4.6))
for ax, (df, col, label) in zip(axes, panels):
    R = df[df["response"] == "Responder"]
    N = df[df["response"] == "Non-responder"]
    p = mannwhitneyu(R[col], N[col], alternative="two-sided").pvalue
    ax.boxplot([R[col], N[col]], positions=[0, 1], widths=0.5, showfliers=False)
    for pos, grp in ((0, R), (1, N)):
        jitter = rng.uniform(-0.12, 0.12, len(grp))
        ax.scatter(pos + jitter, grp[col], s=36, zorder=3, edgecolor="black", linewidth=0.4,
                   c=[THERAPY_COLORS.get(t, "grey") for t in grp["therapy"]])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Responder", "Non-responder"])
    ax.set_ylabel(label)
    ax.set_title(f"n = {len(R)} vs {len(N)} patients\nuncorrected p = {p:.4f}", fontsize=10)
    print(f"{label}: n={len(R)} vs {len(N)}, median {R[col].median():.3f} vs {N[col].median():.3f}, p={p:.4f}")
handles = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=c, markeredgecolor="black", label=t)
           for t, c in THERAPY_COLORS.items()]
fig.legend(handles=handles, title="therapy", loc="center right", bbox_to_anchor=(1.0, 0.5))
fig.suptitle("Baseline (Pre) biopsies: one point per patient", y=1.02)
plt.tight_layout(rect=[0, 0, 0.93, 1])
plt.savefig("docs/figures/baseline_response.png", dpi=150, bbox_inches="tight")
plt.close()
print("wrote docs/figures/baseline_response.png")

# ---- UMAP colored by cell-type label (no patient legend, so nothing overlaps) ----
sc.pl.umap(a, color="celltype", legend_loc="right margin", frameon=False, size=8,
           title="Cell types (working annotation)", show=False, save="_celltypes.png")
print("wrote docs/figures/umap_celltypes.png")