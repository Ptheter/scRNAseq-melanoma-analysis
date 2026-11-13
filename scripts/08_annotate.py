#!/usr/bin/env python
"""Map the 23 Leiden clusters from 03_cluster.py to 17 cell-type labels (working annotation).

The mapping below was made by hand from top marker genes and marker detection fractions of ONE specific
clustering run, so it is only valid if that exact clustering is reproduced. The check below compares the
cluster sizes with the reference run and stops if they differ."""
import pandas as pd, scanpy as sc
pd.set_option("display.width", 200)

# cluster sizes of the reference run (15,509 cells after the pct_mito <= 25 filter)
EXPECTED_SIZES = {"0": 1658, "1": 1357, "2": 1139, "3": 1112, "4": 977, "5": 947, "6": 839, "7": 828,
                  "8": 772, "9": 771, "10": 720, "11": 664, "12": 656, "13": 613, "14": 523, "15": 433,
                  "16": 368, "17": 284, "18": 216, "19": 190, "20": 186, "21": 163, "22": 93}

a = sc.read_h5ad("data/processed/sade_feldman_clustered.h5ad")
sizes = a.obs["leiden"].astype(str).value_counts().to_dict()
if sizes != EXPECTED_SIZES:
    diff = {c: (sizes.get(c), EXPECTED_SIZES.get(c)) for c in sorted(set(sizes) | set(EXPECTED_SIZES), key=int)
            if sizes.get(c) != EXPECTED_SIZES.get(c)}
    raise SystemExit("Cluster sizes differ from the reference run (cluster: found, expected): "
                     f"{diff}\nThe cluster-to-label mapping in this script would be wrong. "
                     "Check package versions against environment.yml and redo the annotation by hand if needed.")

m = {"3":"CD4 Treg","5":"CD4 activated (Tfh-like)","10":"CD8 exhausted-like","11":"CD8 activated effector",
     "18":"CD8 GZMK+ EM","1":"CD8 XCL1+","7":"Proliferating T","20":"Naive/CM T","2":"Memory T (IEG)",
     "8":"IFN-response T","13":"NK/CD8 cytotoxic","0":"CD8-leaning T (unresolved)","4":"CD8-leaning T (unresolved)","6":"CD8-leaning T (unresolved)",
     "19":"CD8-leaning T (unresolved)","9":"Monocyte/DC","14":"Macrophage","12":"B","15":"B","16":"B",
     "21":"Plasma","22":"Plasma","17":"pDC"}
a.obs["celltype"] = a.obs["leiden"].astype(str).map(m).astype("category")
assert a.obs["celltype"].notna().all()
r = a.raw.to_adata()
g = [x for x in ["CD4","CD8A","CD8B"] if x in r.var_names]
X = r[:, g].X; X = X.toarray() if hasattr(X, "toarray") else X
print(pd.DataFrame(X > 0, index=a.obs_names, columns=g).groupby(a.obs["leiden"].values).mean().round(2).to_string())
a.write_h5ad("data/processed/sade_feldman_annotated.h5ad", compression="gzip")
print(a.obs["celltype"].value_counts().to_string())