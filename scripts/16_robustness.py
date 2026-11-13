#!/usr/bin/env python
"""Baseline robustness: (1) TCF7 within anti-PD1 biopsies; (2) memory score recomputed without TCF7."""
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu

NON_T = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
r = a.raw.to_adata(); r.obs = a.obs.copy()
x = r[:, "TCF7"].X
x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
sc.tl.score_genes(r, ["IL7R", "CCR7", "SELL", "LEF1"], score_name="mem_noTCF7", random_state=0)
cols = ["sample", "patient", "timepoint", "response", "therapy", "enrichment", "celltype", "mem_noTCF7"]
o = r.obs[cols].copy()
for c in cols[:-1]:
    o[c] = o[c].astype(str)
o["tcf7"] = x
o["tcf7_on"] = x > 0
o = o[(o["enrichment"] == "CD45+") & (~o["celltype"].isin(NON_T))].copy()
o["mem_noTCF7_adj"] = o["mem_noTCF7"] - o.groupby("celltype")["mem_noTCF7"].transform("mean")

def utest(g, col, label):
    R = g.loc[g["response"] == "Responder", col]; N = g.loc[g["response"] == "Non-responder", col]
    if len(R) < 3 or len(N) < 3:
        print(f"{label:58s} too few patients (n={len(R)} v {len(N)})"); return
    p = mannwhitneyu(R, N, alternative="two-sided").pvalue
    print(f"{label:58s} n={len(R)} v {len(N)}  median {R.median():6.2f} v {N.median():6.2f}  p={p:.4f}")

print("--- (1) TCF7 at Pre, patient-level, >= 20 cells per patient")
pre = o[o["timepoint"] == "Pre"]
subsets = [("all T, all therapies", pre),
           ("all T, anti-PD1 only", pre[pre["therapy"] == "anti-PD1"]),
           ("CD8-leaning (unresolved), anti-PD1 only",
            pre[(pre["therapy"] == "anti-PD1") & (pre["celltype"] == "CD8-leaning T (unresolved)")])]
for lab, d in subsets:
    g = d.groupby(["patient", "response"]).agg(n=("tcf7", "size"), mean=("tcf7", "mean"),
                                               frac=("tcf7_on", "mean")).reset_index()
    g = g[g["n"] >= 20]
    for col in ("mean", "frac"):
        utest(g, col, f"TCF7 [{col}], {lab}")

print("\n--- (2) memory score without TCF7 (IL7R, CCR7, SELL, LEF1), >= 30 cells per biopsy")
KEYS = ["sample", "patient", "timepoint", "response", "therapy"]
mc = ["mem_noTCF7", "mem_noTCF7_adj"]
per = o.groupby(KEYS)[mc].mean(); per["n"] = o.groupby(KEYS).size()
per = per[per["n"] >= 30].reset_index()
pt = per.groupby(["patient", "timepoint", "response", "therapy"])[mc].mean().reset_index()
for tp in ["Pre", "Post"]:
    for lab, d in [("all therapies", pt[pt["timepoint"] == tp]),
                   ("anti-PD1 only", pt[(pt["timepoint"] == tp) & (pt["therapy"] == "anti-PD1")])]:
        for c in mc:
            utest(d, c, f"memory w/o TCF7 [{c}], {tp}, {lab}")
