import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
pd.set_option("display.width", 200)
df = pd.read_csv("results/biopsy_celltype_fractions.csv")
ct = [c for c in df.columns if c not in ("sample", "patient", "timepoint", "response", "therapy")]
pt = df.groupby(["patient", "timepoint", "response", "therapy"])[ct].mean().reset_index()

hits = {"Pre": ["Naive/CM T", "Proliferating T", "B"], "Post": ["IFN-response T", "Proliferating T"]}
for tp, cts in hits.items():
    d = pt[pt["timepoint"] == tp]
    for c in cts:
        print(f"\n{tp} | {c}: per-patient fraction")
        print(d.sort_values(["response", c])[["patient", "response", "therapy", c]].round(3).to_string(index=False))

print("\nanti-PD1 only (uncorrected):")
for tp, cts in hits.items():
    d = pt[(pt["timepoint"] == tp) & (pt["therapy"] == "anti-PD1")]
    R, N = d[d["response"] == "Responder"], d[d["response"] == "Non-responder"]
    for c in cts:
        if len(R) >= 3 and len(N) >= 3:
            p = mannwhitneyu(R[c], N[c], alternative="two-sided").pvalue
            print(f"{tp} {c}: n={len(R)} vs {len(N)}, median {R[c].median():.3f} vs {N[c].median():.3f}, p={p:.3f}")
        else:
            print(f"{tp} {c}: too few patients (n={len(R)} vs {len(N)})")
