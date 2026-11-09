#!/usr/bin/env python
"""Build an AnnData (cells x genes) from the GSE120575 matrix and attach cell metadata.

The GEO matrix holds log2(TPM + 1) values (2**x - 1 sums to ~1e6 per cell), so downstream
steps must NOT apply normalize_total / log1p again.

obs columns:
  sample      biopsy ID from the GEO sample table (48 biopsies), e.g. Pre_P1, Post_P23_2
  timepoint   Pre (baseline) or Post (on treatment)
  patient     patient number (Pre/Post prefix and any _2 suffix removed)
  response    Responder / Non-responder as labelled in GEO (recorded per biopsy)
  therapy     anti-PD1 / anti-CTLA4 / anti-CTLA4+PD1
  enrichment  sort fraction taken from the matrix label: CD45+ (unsorted), T_enriched, myeloid_enriched
"""
import gzip
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

MAT = "data/raw/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META = "data/raw/GSE120575_patient_ID_single_cells.txt.gz"
OUT = "data/processed/sade_feldman_matrix.h5ad"

# ---- header and layout check (fail fast, before the long read) ----
with gzip.open(MAT, "rt") as f:
    cells = f.readline().rstrip("\r\n").split("\t")[1:]
    labels = f.readline().rstrip("\r\n").split("\t")[1:]
    first = f.readline().rstrip("\r\n").split("\t")
n = len(cells)
assert len(labels) == n, f"label row has {len(labels)} fields, expected {n}"
nvals = len(first) - 1
print("cells in header:", n, "| values in first data row:", nvals, flush=True)
if nvals == n + 1 and first[-1] == "":
    print("note: data rows end with a tab; ignoring the extra empty column", flush=True)
elif nvals != n:
    raise SystemExit(f"layout mismatch: {nvals} values vs {n} cells; last fields: {first[-3:]}")
if len(set(cells)) != n:
    raise SystemExit("duplicated cell IDs in the matrix header")

# ---- read the matrix in chunks as sparse float32 ----
genes, blocks = [], []
reader = pd.read_csv(MAT, sep="\t", header=None, skiprows=2, index_col=0,
                     usecols=range(0, n + 1),
                     dtype={i: np.float32 for i in range(1, n + 1)}, chunksize=2000)
for chunk in reader:
    block = chunk.to_numpy(dtype=np.float32)
    if np.isnan(block).any():
        raise SystemExit(f"NaN values found near gene {chunk.index[0]}")
    genes.extend(chunk.index.astype(str))
    blocks.append(sp.csr_matrix(block))
    if len(genes) % 10000 == 0:
        print("genes read:", len(genes), flush=True)
print("genes read (total):", len(genes), flush=True)
X = sp.vstack(blocks, format="csr").T.tocsr()  # cells x genes
del blocks

# ---- cell metadata from the GEO sample table (latin-1: the template has a non-UTF-8 character) ----
rows = []
with gzip.open(META, "rt", encoding="latin-1") as f:
    for line in f:
        p = [x.strip() for x in line.rstrip("\r\n").split("\t")]
        if len(p) >= 7 and p[0].startswith("Sample ") and p[0][7:].strip().isdigit():
            rows.append(p[:7])
meta = pd.DataFrame(rows, columns=["sample_name", "cell", "source", "organism",
                                   "sample", "response", "therapy"])
if meta["cell"].nunique() != len(meta):
    raise SystemExit("duplicated cell IDs in the metadata table")
meta = meta.set_index("cell").reindex(cells)
if meta["response"].isna().any():
    raise SystemExit(f"{int(meta['response'].isna().sum())} cells have no metadata")

# matrix label = GEO sample label + optional sort suffix
suffix = np.array([a[len(b):] if a.startswith(b) else "MISMATCH"
                   for a, b in zip(labels, meta["sample"])])
enrichment = np.select([suffix == "", suffix == "_T_enriched", suffix == "_myeloid_enriched"],
                       ["CD45+", "T_enriched", "myeloid_enriched"], default="other")
if (enrichment == "other").any():
    raise SystemExit("unexpected label suffixes: " + str(sorted(set(suffix[enrichment == "other"]))))

obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
obs["sample"] = meta["sample"].to_numpy()
obs["timepoint"] = obs["sample"].str.split("_").str[0]
obs["patient"] = obs["sample"].str.extract(r"^(?:Pre|Post)_(P\d+)", expand=False)
obs["response"] = meta["response"].to_numpy()
obs["therapy"] = meta["therapy"].to_numpy()
obs["enrichment"] = enrichment
assert obs["timepoint"].isin(["Pre", "Post"]).all(), "unexpected timepoint labels"
assert obs["patient"].notna().all(), "could not parse a patient number"
for col in ("response", "therapy"):
    assert (obs.groupby("sample")[col].nunique() == 1).all(), f"{col} varies within a biopsy"
for col in ("timepoint", "patient", "response", "therapy", "enrichment", "sample"):
    obs[col] = obs[col].astype("category")

adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
adata.var_names_make_unique()

# ---- diagnostics ----
print("shape (cells x genes):", adata.shape)
print("nonzero fraction:", round(X.nnz / (X.shape[0] * X.shape[1]), 4), "| max value:", round(float(X.max()), 2))
s = obs.drop_duplicates("sample")
print("\nbiopsies:", len(s), "| patients:", s["patient"].nunique())
print(pd.crosstab(s["response"], s["timepoint"]))
print("\nbiopsies by therapy:\n", s["therapy"].value_counts().to_string())
print("\npatients with both Pre and Post:", int((s.groupby("patient", observed=True)["timepoint"].nunique() == 2).sum()))
inc = s.groupby("patient", observed=True)["response"].nunique()
print("patients whose biopsies carry both response labels:", sorted(inc[inc > 1].index.astype(str)))
print("\ncells by sort fraction:\n", obs["enrichment"].value_counts().to_string())

adata.write_h5ad(OUT, compression="gzip")
print("wrote", OUT)
