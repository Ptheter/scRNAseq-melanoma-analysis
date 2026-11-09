#!/usr/bin/env python
"""QC metrics for the GSE120575 matrix (log2(TPM+1)). Computes metrics and plots; removes no cells."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad("data/processed/sade_feldman_matrix.h5ad")
X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

# genes detected per cell (value > 0)
n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

# back-transform to TPM to get the mitochondrial share
T = X.copy()
T.data = np.expm1(T.data * np.log(2)).astype(np.float32)
total = np.asarray(T.sum(axis=1)).ravel()
mt_idx = np.where(adata.var_names.str.startswith("MT-"))[0]
print("MT- genes found:", len(mt_idx), list(adata.var_names[mt_idx][:15]))
pct_mito = np.asarray(T[:, mt_idx].sum(axis=1)).ravel() / total * 100

qc = adata.obs[["sample", "timepoint", "patient", "response", "enrichment"]].copy()
qc["n_genes"] = n_genes
qc["pct_mito"] = pct_mito
if "PTPRC" in adata.var_names:
    qc["PTPRC_detected"] = np.asarray((X[:, adata.var_names.get_loc("PTPRC")] > 0).todense()).ravel()
qc.to_csv("results/qc_metrics.csv")

q = [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
print("\nquantiles:\n", qc[["n_genes", "pct_mito"]].quantile(q).round(2).to_string())
print("\nmedians by sort fraction:\n",
      qc.groupby("enrichment", observed=True)[["n_genes", "pct_mito"]].median().round(2).to_string())
for t in (500, 1000, 1500):
    print(f"cells with n_genes < {t}:", int((qc["n_genes"] < t).sum()))
for t in (5, 10, 20):
    print(f"cells with pct_mito > {t}:", int((qc["pct_mito"] > t).sum()))
if "PTPRC_detected" in qc:
    print("fraction of cells with PTPRC (CD45) detected:", round(float(qc["PTPRC_detected"].mean()), 3))

per = qc.groupby("sample", observed=True).agg(cells=("n_genes", "size"),
                                              med_genes=("n_genes", "median"),
                                              med_mito=("pct_mito", "median"))
print("\nsamples with the lowest median genes/cell:\n", per.sort_values("med_genes").head(6).round(1).to_string())
print("\nsamples with the highest median mito %:\n", per.sort_values("med_mito", ascending=False).head(6).round(1).to_string())
print("\ncells per sample: min", int(per["cells"].min()), "median", int(per["cells"].median()), "max", int(per["cells"].max()))

fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].hist(qc["n_genes"], bins=60); ax[0].set_xlabel("genes detected per cell"); ax[0].set_ylabel("cells")
ax[1].hist(qc["pct_mito"], bins=60); ax[1].set_xlabel("% mitochondrial (of TPM)")
ax[2].scatter(qc["n_genes"], qc["pct_mito"], s=2, alpha=0.3)
ax[2].set_xlabel("genes detected per cell"); ax[2].set_ylabel("% mitochondrial")
plt.tight_layout(); plt.savefig("docs/figures/qc_overview.png", dpi=150); plt.close()

order = per.sort_values("med_genes").index
data = [qc.loc[qc["sample"] == s, "n_genes"].to_numpy() for s in order]
fig, ax = plt.subplots(figsize=(14, 5))
ax.boxplot(data, showfliers=False)
ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=90, fontsize=7)
ax.set_ylabel("genes detected per cell"); ax.set_title("Genes detected per cell, by biopsy")
plt.tight_layout(); plt.savefig("docs/figures/qc_by_sample.png", dpi=150); plt.close()
print("\nwrote results/qc_metrics.csv and two figures in docs/figures/")
