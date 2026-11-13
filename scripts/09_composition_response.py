import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 250)

a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
o = a.obs[["sample", "patient", "timepoint", "response", "therapy", "enrichment", "celltype"]].copy()
for c in o.columns:
    o[c] = o[c].astype(str)

# unsorted CD45+ cells only, so T-enrichment does not inflate T-subset fractions
o = o[o["enrichment"] == "CD45+"]
# exclude the one-patient macrophage group from comparisons
o = o[o["celltype"] != "Macrophage"]

counts = o.groupby(["sample", "celltype"]).size().unstack(fill_value=0)
tot = counts.sum(axis=1)
keep = tot[tot >= 50].index                      # drop biopsies with too few cells for stable fractions
print(f"biopsies kept (>=50 CD45+ cells): {len(keep)} of {len(tot)}")
frac = counts.loc[keep].div(tot[keep], axis=0)

meta = o.drop_duplicates("sample").set_index("sample").loc[keep, ["patient", "timepoint", "response", "therapy"]]
df = frac.join(meta)
df.to_csv("results/biopsy_celltype_fractions.csv")

print("\nbiopsies per group:\n", df.groupby(["timepoint", "response"]).size().to_string())
print("\ntherapy by group:\n", pd.crosstab([df["timepoint"], df["response"]], df["therapy"]).to_string())

# one biopsy per patient per timepoint: average repeat biopsies (_2) within patient
pt = df.groupby(["patient", "timepoint", "response", "therapy"], observed=True)[list(frac.columns)].mean().reset_index()

rows = []
for tp in ["Pre", "Post"]:
    d = pt[pt["timepoint"] == tp]
    R, N = d[d["response"] == "Responder"], d[d["response"] == "Non-responder"]
    for ct in frac.columns:
        x, y = R[ct].to_numpy(), N[ct].to_numpy()
        if len(x) < 3 or len(y) < 3:
            continue
        p = mannwhitneyu(x, y, alternative="two-sided").pvalue
        rows.append({"timepoint": tp, "celltype": ct, "n_resp": len(x), "n_nonresp": len(y),
                     "median_resp": np.median(x), "median_nonresp": np.median(y), "p": p})
res = pd.DataFrame(rows)
res["fdr"] = np.nan
for tp in res["timepoint"].unique():
    m = res["timepoint"] == tp
    res.loc[m, "fdr"] = multipletests(res.loc[m, "p"], method="fdr_bh")[1]
res = res.sort_values("p")
res.to_csv("results/composition_response_tests.csv", index=False)
print("\nResponder vs Non-responder, per patient:\n", res.round(4).to_string(index=False))
