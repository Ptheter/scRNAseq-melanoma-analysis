#!/usr/bin/env python
"""Reconcile sample labels between the matrix header and the GEO metadata table; test the TPM scale."""
import gzip
import numpy as np
import pandas as pd

MAT = "data/raw/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META = "data/raw/GSE120575_patient_ID_single_cells.txt.gz"

with gzip.open(MAT, "rt") as f:
    cells = f.readline().rstrip("\r\n").split("\t")[1:]
    labels = f.readline().rstrip("\r\n").split("\t")[1:]
mat = pd.DataFrame({"cell": cells, "label_matrix": labels})

rows = []
with gzip.open(META, "rt", encoding="latin-1") as f:
    for line in f:
        p = [x.strip() for x in line.rstrip("\r\n").split("\t")]
        if len(p) >= 7 and p[0].startswith("Sample ") and p[0][7:].strip().isdigit():
            rows.append(p[:7])
meta = pd.DataFrame(rows, columns=["sample_name", "cell", "source", "organism",
                                   "label_meta", "response", "therapy"])
m = mat.merge(meta[["cell", "label_meta", "response", "therapy"]], on="cell")

print("cells:", len(m))
print("distinct labels: matrix =", m["label_matrix"].nunique(), "| metadata =", m["label_meta"].nunique())
pat = lambda s: s.str.split("_", n=1).str[1]
print("distinct patient strings: matrix =", pat(m["label_matrix"]).nunique(), "| metadata =", pat(m["label_meta"]).nunique())

dis = m[m["label_matrix"] != m["label_meta"]]
print("\ndisagreeing cells:", len(dis))
print(dis.groupby(["label_matrix", "label_meta"]).size().sort_values(ascending=False).head(20).to_string())

print("\nmetadata labels (cells | response | therapy):")
g = m.groupby("label_meta").agg(cells=("cell", "size"),
                                response=("response", lambda s: "/".join(sorted(set(s)))),
                                therapy=("therapy", lambda s: "/".join(sorted(set(s)))))
print(g.sort_index().to_string())

import scanpy as sc
a = sc.read_h5ad("data/processed/sade_feldman_matrix.h5ad")
X = a[:500].X.copy()
X.data = np.expm1(X.data * np.log(2))  # 2**x - 1
print("\nsum of 2**x - 1 for the first 5 cells (about 1,000,000 if the data is log2(TPM+1)):")
print(np.round(np.asarray(X.sum(axis=1)).ravel()[:5]))
