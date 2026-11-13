import pandas as pd, scanpy as sc
pd.set_option("display.width", 200)
a = sc.read_h5ad("data/processed/sade_feldman_clustered.h5ad")
print(pd.crosstab(a.obs["leiden"], a.obs["enrichment"]))
print(a.obs.groupby("leiden", observed=True)[["n_genes", "pct_mito"]].median().round(1))
# which biopsy dominates the patient-skewed clusters?
for c in ["13", "19", "9"]:
    print(c, a.obs[a.obs.leiden == c]["sample"].value_counts().head(4).to_dict())
