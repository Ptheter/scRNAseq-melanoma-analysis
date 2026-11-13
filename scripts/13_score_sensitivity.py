import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
r = a.raw.to_adata(); r.obs = a.obs.copy()
sets = {
  "memory":        ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
  "exhaustion":    ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "CXCL13"],
  "cytotoxic":     ["GZMB", "PRF1", "NKG7", "GNLY", "GZMH"],
  "proliferation": ["MKI67", "TOP2A", "STMN1", "TUBB", "RRM2"],
  "interferon":    ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "RSAD2", "IFI44L"],
}
for k, v in sets.items():
    sc.tl.score_genes(r, [g for g in v if g in r.var_names], score_name=f"score_{k}", random_state=0)
cols = [f"score_{k}" for k in sets]
o = r.obs.copy()
for c in ["sample", "patient", "timepoint", "response", "therapy", "enrichment", "celltype"]:
    o[c] = o[c].astype(str)
o = o[o["enrichment"] == "CD45+"]
non_t = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
T = o[~o["celltype"].isin(non_t)]
key = ["sample", "patient", "timepoint", "response", "therapy"]

def patient_table(d, min_cells):
    per = d.groupby(key)[cols].mean()
    per["n"] = d.groupby(key).size()
    per = per[per["n"] >= min_cells].reset_index()
    return per.groupby(["patient", "timepoint", "response", "therapy"])[cols].mean().reset_index()

def test(pt, label):
    for tp in ["Pre", "Post"]:
        x = pt[pt["timepoint"] == tp]
        R, N = x[x["response"] == "Responder"], x[x["response"] == "Non-responder"]
        if len(R) < 3 or len(N) < 3:
            print(f"{label} {tp}: too few (n={len(R)} vs {len(N)})"); continue
        for c in cols:
            p = mannwhitneyu(R[c], N[c], alternative="two-sided").pvalue
            print(f"{label:30s} {tp:4s} {c[6:]:14s} n={len(R)} v {len(N)}  "
                  f"median {R[c].median():6.2f} v {N[c].median():6.2f}  p={p:.4f}")

# 1. one broad cell group only, so composition cannot drive the score
test(patient_table(o[o["celltype"] == "CD8-leaning T (unresolved)"], 20), "within CD8-leaning T only")
# 2. all T cells, dropping the two patients with opposite labels at Post
pt = patient_table(T, 30)
test(pt[~((pt["timepoint"] == "Post") & pt["patient"].isin(["P1", "P5"]))], "all T, Post minus P1/P5")
