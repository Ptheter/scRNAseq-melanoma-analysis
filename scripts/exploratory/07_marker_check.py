import pandas as pd, scanpy as sc
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
a = sc.read_h5ad("data/processed/sade_feldman_clustered.h5ad")
r = a.raw.to_adata()
genes = [g for g in ["PDCD1","LAG3","TIGIT","CXCL13","HAVCR2","CTLA4","FOXP3","TCF7","CCR7",
                     "GZMK","MLANA","PMEL","SOX10","COL1A1","CD3E","C1QA","MS4A1","MZB1"] if g in r.var_names]
X = r[:, genes].X
X = X.toarray() if hasattr(X, "toarray") else X
df = pd.DataFrame(X, index=r.obs_names, columns=genes)
frac = df.gt(0).groupby(a.obs["leiden"].values).mean()
print("fraction of cells expressing:\n", frac.loc[sorted(frac.index, key=int)].round(2).to_string())
print("\nPre/Post fraction by cluster:\n",
      pd.crosstab(a.obs["leiden"], a.obs["timepoint"], normalize="index").round(2).to_string())
