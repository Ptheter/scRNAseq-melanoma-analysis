import pandas as pd, scanpy as sc
pd.set_option("display.width", 200)
a = sc.read_h5ad("data/processed/sade_feldman_clustered.h5ad")
n_all = a.obs["sample"].value_counts()
for c in ["11", "12", "14", "15", "16", "19", "21", "22"]:
    n_in = a.obs[a.obs["leiden"] == c]["sample"].value_counts()
    frac = (n_in / n_all[n_in.index]).sort_values(ascending=False)
    print(f"cluster {c}: biopsy share of own cells {frac.head(3).round(2).to_dict()} | "
          f"top biopsy counts {n_in.head(3).to_dict()}")
print()
print(pd.crosstab(a.obs["leiden"], a.obs["enrichment"]))
