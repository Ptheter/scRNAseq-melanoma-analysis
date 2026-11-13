import pandas as pd, scanpy as sc
pd.set_option("display.width", 200)
a = sc.read_h5ad("data/processed/sade_feldman_annotated.h5ad")
o = a.obs[["sample", "patient", "timepoint", "response", "therapy"]].astype(str)
b = o.drop_duplicates("sample").sort_values(["patient", "timepoint", "sample"])
b["cells"] = b["sample"].map(o["sample"].value_counts())

print("biopsies per patient/timepoint, with labels:")
print(b.to_string(index=False))

for col in ["response", "therapy"]:
    within = b.groupby(["patient", "timepoint"])[col].nunique()
    print(f"\n{col} differs between biopsies of the SAME patient and timepoint:",
          [f"{p}/{t}" for (p, t), n in within.items() if n > 1])
    across = b.groupby("patient")[col].nunique()
    print(f"{col} differs between Pre and Post for the same P-number:", sorted(across[across > 1].index))
