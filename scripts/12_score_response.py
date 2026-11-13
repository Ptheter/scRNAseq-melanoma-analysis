import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 220)

a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
r = a.raw.to_adata()
r.obs = a.obs.copy()                      # raw = all genes, log2(TPM+1)

sets = {
  "memory":        ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
  "exhaustion":    ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "CXCL13"],
  "cytotoxic":     ["GZMB", "PRF1", "NKG7", "GNLY", "GZMH"],
  "proliferation": ["MKI67", "TOP2A", "STMN1", "TUBB", "RRM2"],
  "interferon":    ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "RSAD2", "IFI44L"],
}
for k, v in sets.items():
    present = [g for g in v if g in r.var_names]
    print(k, "genes used:", present, "| missing:", sorted(set(v) - set(present)))
    sc.tl.score_genes(r, present, score_name=f"score_{k}", random_state=0)

non_t = ["Monocyte/DC", "Macrophage", "B", "Plasma", "pDC", "NK/CD8 cytotoxic"]
o = r.obs.copy()
for c in ["sample", "patient", "timepoint", "response", "therapy", "enrichment", "celltype"]:
    o[c] = o[c].astype(str)
o = o[(o["enrichment"] == "CD45+") & (~o["celltype"].isin(non_t))]

sc_cols = [f"score_{k}" for k in sets]
per = o.groupby(["sample", "patient", "timepoint", "response", "therapy"])[sc_cols].mean()
per["n_T"] = o.groupby(["sample", "patient", "timepoint", "response", "therapy"]).size()
per = per[per["n_T"] >= 30].reset_index()
print("\nbiopsies kept (>=30 unsorted T cells):", len(per))

# same unit as script 09: one row per patient/timepoint/response/therapy
pt = per.groupby(["patient", "timepoint", "response", "therapy"])[sc_cols].mean().reset_index()
pt.to_csv("results/patient_scores.csv", index=False)

def run(d, label):
    rows = []
    for tp in ["Pre", "Post"]:
        x = d[d["timepoint"] == tp]
        R, N = x[x["response"] == "Responder"], x[x["response"] == "Non-responder"]
        if len(R) < 3 or len(N) < 3:
            print(f"{label} {tp}: too few (n={len(R)} vs {len(N)})"); continue
        for c in sc_cols:
            p = mannwhitneyu(R[c], N[c], alternative="two-sided").pvalue
            rows.append({"set": label, "timepoint": tp, "score": c.replace("score_", ""),
                         "n_resp": len(R), "n_nonresp": len(N),
                         "median_resp": R[c].median(), "median_nonresp": N[c].median(), "p": p})
    t = pd.DataFrame(rows)
    if len(t):
        t["fdr"] = np.nan
        for tp in t["timepoint"].unique():
            m = t["timepoint"] == tp
            t.loc[m, "fdr"] = multipletests(t.loc[m, "p"], method="fdr_bh")[1]
    return t

res = pd.concat([run(pt, "all therapies"), run(pt[pt["therapy"] == "anti-PD1"], "anti-PD1 only")])
res.to_csv("results/score_response_tests.csv", index=False)
print("\n", res.round(4).to_string(index=False))
