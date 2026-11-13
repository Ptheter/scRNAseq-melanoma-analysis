#!/usr/bin/env python
"""Gene-set scores adjusted for cell type: residual = cell score - mean score of its cell type.
Unit = patient (per-biopsy means, >= 30 cells; then averaged per patient/timepoint/response/therapy).
Mann-Whitney within timepoint; BH across the 5 scores within each analysis x timepoint."""
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 220)

MIN_CELLS = 30
sets = {
  "memory":        ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
  "exhaustion":    ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "CXCL13"],
  "cytotoxic":     ["GZMB", "PRF1", "NKG7", "GNLY", "GZMH"],
  "proliferation": ["MKI67", "TOP2A", "STMN1", "TUBB", "RRM2"],
  "interferon":    ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "RSAD2", "IFI44L"],
}
a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
r = a.raw.to_adata(); r.obs = a.obs.copy()
for k, v in sets.items():
    present = [g for g in v if g in r.var_names]
    print(k, "used:", present, "| missing:", sorted(set(v) - set(present)))
    sc.tl.score_genes(r, present, score_name=f"score_{k}", random_state=0)

non_t = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
o = r.obs.copy()
for c in ["sample", "patient", "timepoint", "response", "therapy", "enrichment", "celltype"]:
    o[c] = o[c].astype(str)
o = o[(o["enrichment"] == "CD45+") & (~o["celltype"].isin(non_t))].copy()
raw_cols = [f"score_{k}" for k in sets]
adj_cols = [f"adj_{k}" for k in sets]
for k in sets:   # residual relative to the cell type's mean (pooled over all included cells)
    o[f"adj_{k}"] = o[f"score_{k}"] - o.groupby("celltype")[f"score_{k}"].transform("mean")
print("\nT cells used:", len(o), "| cell types:", o["celltype"].nunique())
o[raw_cols + adj_cols + ["sample", "patient", "timepoint", "response", "therapy", "celltype"]].to_csv("results/cell_scores_T.csv")

KEYS = ["sample", "patient", "timepoint", "response", "therapy"]
def agg(d, cols):
    per = d.groupby(KEYS)[cols].mean()
    per["n"] = d.groupby(KEYS).size()
    per = per[per["n"] >= MIN_CELLS].reset_index()
    return per.groupby(["patient", "timepoint", "response", "therapy"])[cols].mean().reset_index()

# sanity check: reproduce script 12's patient-level scores
try:
    old = pd.read_csv("results/patient_scores.csv")
    new = agg(o, raw_cols)
    m = new.merge(old, on=["patient", "timepoint", "response", "therapy"], suffixes=("", "_old"))
    print("patient rows: new", len(new), "| old", len(old), "| matched", len(m))
    print("max abs difference vs results/patient_scores.csv:",
          {c: round(float((m[c] - m[c + "_old"]).abs().max()), 6) for c in raw_cols})
except FileNotFoundError:
    print("results/patient_scores.csv not found; skipping the reproduction check")

def test(pt, cols, label):
    rows = []
    for tp in ["Pre", "Post"]:
        x = pt[pt["timepoint"] == tp]
        R, N = x[x["response"] == "Responder"], x[x["response"] == "Non-responder"]
        if len(R) < 3 or len(N) < 3:
            print(f"{label} {tp}: too few patients (n={len(R)} vs {len(N)})"); continue
        for c in cols:
            p = mannwhitneyu(R[c], N[c], alternative="two-sided").pvalue
            rows.append({"analysis": label, "tp": tp, "score": c.split("_", 1)[1], "n_R": len(R), "n_NR": len(N),
                         "med_R": R[c].median(), "med_NR": N[c].median(), "p": p})
    t = pd.DataFrame(rows)
    if len(t):
        t["fdr"] = np.nan
        for tp in t["tp"].unique():
            mk = t["tp"] == tp
            t.loc[mk, "fdr"] = multipletests(t.loc[mk, "p"], method="fdr_bh")[1]
    return t

cd8 = o["celltype"].str.startswith("CD8")
res = pd.concat([
    test(agg(o, raw_cols), raw_cols, "all T, raw score"),
    test(agg(o, adj_cols), adj_cols, "all T, cluster-adjusted"),
    test(agg(o[o["therapy"] == "anti-PD1"], adj_cols), adj_cols, "all T, adjusted, anti-PD1 only"),
    test(agg(o[cd8], adj_cols), adj_cols, "CD8-labelled clusters, adjusted"),
])
res.to_csv("results/cluster_adjusted_tests.csv", index=False)
print("\n", res.round(4).to_string(index=False))
