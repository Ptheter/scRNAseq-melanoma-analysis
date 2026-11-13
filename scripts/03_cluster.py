#!/usr/bin/env python
"""Filter, embed and cluster the GSE120575 cells. Data is log2(TPM+1): no normalize_total / log1p.
usage: python scripts/03_cluster.py [max_pct_mito=25] [leiden_resolution=0.8]"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import scanpy as sc
from scipy.sparse import csr_array

if not hasattr(csr_array, "getnnz"):
    csr_array.getnnz = csr_array._getnnz

mito_max = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
res = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8
sc.settings.figdir = "docs/figures"
sc.settings.verbosity = 1

adata = sc.read_h5ad("data/processed/sade_feldman_matrix.h5ad")
qc = pd.read_csv("results/qc_metrics.csv", index_col=0).loc[adata.obs_names]
adata.obs["n_genes"] = qc["n_genes"].to_numpy()
adata.obs["pct_mito"] = qc["pct_mito"].to_numpy()

n0 = adata.n_obs
adata = adata[adata.obs["pct_mito"] <= mito_max].copy()
print(f"cells: {n0} -> {adata.n_obs} after pct_mito <= {mito_max} (removed {n0 - adata.n_obs})")
sc.pp.filter_genes(adata, min_cells=3)
print("genes kept (detected in >= 3 cells):", adata.n_vars)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
# exclude TCR/BCR variable, joining and diversity segments from the HVG set (clonotype signal, not cell state)
vdj = adata.var_names.str.match(r"^(TR[ABGD][VDJ]|IG[HKL][VDJ])")
print("V(D)J genes flagged as HVG and dropped:", int((adata.var["highly_variable"] & vdj).sum()))
adata.var.loc[vdj, "highly_variable"] = False
noise = adata.var_names.str.match(
    r"^(MTRNR2L\d+|RN7S[KL]\d*|HSP[AB]\d|DNAJ[ABC]\d|FOS|FOSB|JUN|JUNB|JUND|"
    r"RP\d+-|RP11-|AC\d{6}\.|AL\d{6}\.|CT[ABCD]-)"
)
print("noise genes flagged as HVG and dropped:", int((adata.var["highly_variable"] & noise).sum()))
adata.var.loc[noise, "highly_variable"] = False

adata.raw = adata                      # keep all genes for marker plots and DE
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=0)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=res, key_added="leiden", random_state=0)
print("clusters:", adata.obs["leiden"].nunique())

sc.pl.umap(adata, color=["leiden", "patient", "response", "timepoint", "enrichment", "therapy",
                         "n_genes", "pct_mito"], ncols=4, wspace=0.4, show=False, save="_overview.png")

# ---- cluster sizes and how many patients/biopsies contribute (dominance by one patient = batch warning) ----
g = adata.obs.groupby("leiden", observed=True)
summ = pd.DataFrame({"cells": g.size(), "patients": g["patient"].nunique(), "biopsies": g["sample"].nunique()})
top = g["patient"].agg(lambda s: s.value_counts(normalize=True).iloc[0])
summ["top_patient_share"] = top.round(2)
print("\ncluster summary:\n", summ.to_string())

# ---- marker dotplot ----
markers = {
  "T cell": ["CD3D", "CD3E", "CD4", "CD8A", "CD8B"],
  "naive/memory": ["IL7R", "CCR7", "SELL", "TCF7"],
  "cytotoxic": ["GZMB", "GZMK", "PRF1", "NKG7", "GNLY"],
  "exhaustion": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4"],
  "Treg": ["FOXP3", "IL2RA"],
  "NK": ["KLRD1", "NCAM1"],
  "B/plasma": ["MS4A1", "CD79A", "CD19", "MZB1"],
  "myeloid": ["CD14", "LYZ", "CD68", "CD163", "C1QA", "TREM2"],
  "DC": ["CLEC9A", "CD1C", "LILRA4"],
  "proliferation": ["MKI67", "TOP2A"],
  "non-immune": ["MLANA", "PMEL", "SOX10", "COL1A1"],
}
present = {k: [m for m in v if m in adata.raw.var_names] for k, v in markers.items()}
missing = [m for v in markers.values() for m in v if m not in adata.raw.var_names]
print("\nmarkers not found:", missing)
present = {k: v for k, v in present.items() if v}
sc.pl.dotplot(adata, present, groupby="leiden", standard_scale="var", show=False, save="_markers.png")

# ---- top genes per cluster ----
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True, n_genes=25)
deg = sc.get.rank_genes_groups_df(adata, group=None)
deg.to_csv("results/cluster_markers.csv", index=False)
print("\ntop 6 genes per cluster:")
for c in sorted(adata.obs["leiden"].unique(), key=int):
    print(c, ", ".join(deg[deg["group"] == c].head(6)["names"]))

adata.write_h5ad("data/processed/sade_feldman_clustered.h5ad", compression="gzip")
print("\nwrote data/processed/sade_feldman_clustered.h5ad")
