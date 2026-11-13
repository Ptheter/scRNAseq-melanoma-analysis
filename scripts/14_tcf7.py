#!/usr/bin/env python
"""TCF7 per patient at Pre (baseline), responders vs non-responders. Motivated by the original paper.
Unsorted CD45+ cells only; patient-level medians; Mann-Whitney U."""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

MIN_CELLS = 20
NON_T = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
assert a.raw is not None and "TCF7" in a.raw.var_names, "TCF7 not in .raw"
x = a.raw[:, "TCF7"].X
x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
o = a.obs[["patient", "timepoint", "response", "enrichment", "celltype"]].copy()
o["tcf7"] = x
o["tcf7_on"] = x > 0
print("labels present:", sorted(o["celltype"].astype(str).unique()))

def run(d, name):
    d = d[(d["timepoint"] == "Pre") & (d["enrichment"] == "CD45+")]
    g = d.groupby(["patient", "response"], observed=True).agg(n=("tcf7", "size"), mean=("tcf7", "mean"),
                                                              frac=("tcf7_on", "mean")).reset_index()
    g = g[g["n"] >= MIN_CELLS]
    for col in ("mean", "frac"):
        r = g.loc[g["response"] == "Responder", col]; n = g.loc[g["response"] == "Non-responder", col]
        if len(r) < 3 or len(n) < 3:
            print(f"{name:34s} {col:5s} too few patients (n={len(r)} v {len(n)})"); continue
        p = mannwhitneyu(r, n, alternative="two-sided").pvalue
        print(f"{name:34s} {col:5s} n={len(r)} v {len(n)}  median {r.median():5.2f} v {n.median():5.2f}  p={p:.4f}")

run(o[~o["celltype"].isin(NON_T)], "all T")
run(o[o["celltype"] == "CD8-leaning T (unresolved)"], "CD8-leaning T (unresolved)")
for ct in sorted(set(o["celltype"].astype(str)) - set(NON_T)):
    if ct.startswith("CD8") and "unresolved" not in ct:
        run(o[o["celltype"] == ct], ct)
