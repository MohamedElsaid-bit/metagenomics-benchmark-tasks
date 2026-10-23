"""Task 1 exploration: does the T2D microbiome signal in MetaCardis survive metformin stratification?"""
import sys
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA = sys.argv[1]

meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0]
meta = meta[meta.study_name == "MetaCardis_2020_a"].set_index("sample_id")
ab = list(pyreadr.read_r(f"{DATA}/MetaCardis_2020_a.relative_abundance.rda").values())[0]
if "rowname" in ab.columns:
    ab = ab.set_index("rowname")

# species-level rows only
sp = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")].copy()
sp.index = sp.index.str.split("|s__", regex=False).str[-1]
print("duplicate names:", sp.index[sp.index.duplicated()].tolist()[:10])
sp = sp.groupby(level=0).sum()
sp = sp / sp.sum() * 100  # renormalize to species total
print("species x samples:", sp.shape)

meta = meta.loc[sp.columns]
tx = meta.treatment.astype(str)
meta["metformin"] = tx.str.contains("metformin")
ANTIDIAB = {"antidiab", "metformin", "insulin", "bolus", "dppiv", "glp_1", "su", "sglt2", "glitazone", "acarbose"}
meta["any_antidiab"] = tx.str.split(";").apply(lambda t: bool(ANTIDIAB & set(t)))
dis = meta.disease.astype(str)

groups = {
    "healthy": dis == "healthy",
    "T2D_all": dis.str.contains("T2D"),
}
groups["T2D_met"] = groups["T2D_all"] & meta.metformin
groups["T2D_nomet"] = groups["T2D_all"] & ~meta.metformin
groups["T2D_untreated"] = groups["T2D_all"] & ~meta.any_antidiab
groups["T2D_only_met"] = (dis == "T2D") & meta.metformin
groups["T2D_only_nomet"] = (dis == "T2D") & ~meta.metformin
for g, m in groups.items():
    sub = meta[m]
    print(f"{g:16s} n={m.sum():4d} age={sub.age.mean():.1f} BMI={sub.BMI.mean():.1f}")

print("\ntreatment strings among T2D without metformin (top 15):")
print(tx[groups["T2D_nomet"]].value_counts().head(15))


def compare(a, b, prev=0.10, q=0.05):
    x, y = sp.loc[:, groups[a]], sp.loc[:, groups[b]]
    both = pd.concat([x, y], axis=1)
    keep = (both > 0).mean(axis=1) >= prev
    rows = []
    for s in both.index[keep]:
        p = mannwhitneyu(x.loc[s], y.loc[s]).pvalue
        rows.append((s, p, x.loc[s].mean(), y.loc[s].mean(), (x.loc[s] > 0).mean(), (y.loc[s] > 0).mean()))
    r = pd.DataFrame(rows, columns=["species", "p", f"mean_{a}", f"mean_{b}", f"prev_{a}", f"prev_{b}"]).set_index("species")
    r["q"] = multipletests(r.p, method="fdr_bh")[1]
    sig = r[r.q < q]
    print(f"\n== {a} vs {b}: tested {keep.sum()}, significant (BH q<{q}) {len(sig)}")
    for s in ["Escherichia_coli", "Intestinibacter_bartlettii", "Roseburia_intestinalis", "Roseburia_hominis", "Anaerostipes_hadrus",
              "Faecalibacterium_prausnitzii", "Akkermansia_muciniphila"]:
        if s in r.index:
            print(f"   {s:30s} q={r.loc[s,'q']:.2e} {r.loc[s].iloc[1]:.3f} vs {r.loc[s].iloc[2]:.3f}")
    return r


res = {}
for a, b in [("T2D_all", "healthy"), ("T2D_met", "healthy"), ("T2D_nomet", "healthy"),
             ("T2D_untreated", "healthy"), ("T2D_met", "T2D_nomet"),
             ("T2D_only_met", "healthy"), ("T2D_only_nomet", "healthy")]:
    res[(a, b)] = compare(a, b)
