import pandas as pd, scanpy as sc
pd.set_option("display.width", 200)
a = sc.read_h5ad("data/processed/sade_feldman_clustered.h5ad")
for c in ["0", "1", "6", "9", "14", "15", "21"]:
    inc = a.obs["leiden"] == c
    n_in = a.obs[inc]["sample"].value_counts()
    n_all = a.obs["sample"].value_counts()
    frac = (n_in / n_all[n_in.index]).dropna().sort_values(ascending=False)
    print(f"cluster {c}: top biopsies by share of THAT biopsy's cells in the cluster")
    print(frac.head(4).round(2).to_dict(), "| biopsy sizes:", n_all[frac.head(4).index].to_dict(), "\n")
