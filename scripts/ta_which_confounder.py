"""Task A: which single adjustment removes the marginal statin / sulfonylurea / insulin signal?"""
import sys
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA = sys.argv[1]
meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0]
meta = meta[meta.study_name == "MetaCardis_2020_a"].set_index("sample_id")
ab = list(pyreadr.read_r(f"{DATA}/MetaCardis_2020_a.relative_abundance.rda").values())[0]
if "rowname" in ab.columns:
    ab = ab.set_index("rowname")
sp = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")].copy()
sp.index = sp.index.str.split("|s__", regex=False).str[-1]
sp = sp.groupby(level=0).sum()
sp = sp / sp.sum() * 100
meta = meta.loc[sp.columns]
tok = meta.treatment.astype(str).str.split(";")
DRUGS = {"metformin": "metformin", "proton_pump_inhibitor": "PPI", "sulfonylurea": "su",
         "statin": "statin", "insulin": "insulin", "dpp4": "dppiv", "aspirin": "asa"}
for n, t in DRUGS.items():
    meta[n] = tok.apply(lambda x, t=t: t in x).astype(int)
m = meta[meta.disease.astype(str).str.contains("T2D")].copy()
X = sp.loc[:, m.index]
keep = X.index[(X > 0).mean(axis=1) >= 0.10]
pc = X[X > 0].min().min() / 2
lg = np.log(X + pc)
clr = lg.sub(lg.mean(axis=0), axis=1)


def nsig(p):
    return int((multipletests(p, method="fdr_bh")[1] < 0.05).sum())


print("co-prescription with metformin:")
for d in ["statin", "sulfonylurea", "insulin", "proton_pump_inhibitor"]:
    u, nu = m.index[m[d] == 1], m.index[m[d] == 0]
    print(f"  {d:22s} users on metformin {m.loc[u,'metformin'].mean():.3f}  non-users on metformin {m.loc[nu,'metformin'].mean():.3f}")

for target in ["statin", "sulfonylurea", "insulin"]:
    a, b = m.index[m[target] == 1], m.index[m[target] == 0]
    marg = nsig([mannwhitneyu(X.loc[s, a], X.loc[s, b]).pvalue for s in keep])
    print(f"\n{target}: marginal {marg}")
    for adj in ["metformin", "proton_pump_inhibitor", "dpp4", "aspirin", "age+BMI", "country"]:
        d = m[[target]].copy()
        if adj == "age+BMI":
            d = d.join(m[["age", "BMI"]]); form = f"y ~ {target} + age + BMI"
        elif adj == "country":
            d = d.join(m[["country"]]); form = f"y ~ {target} + C(country)"
        else:
            d = d.join(m[[adj]]); form = f"y ~ {target} + {adj}"
        d = d.dropna()
        p = []
        for s in keep:
            d["y"] = clr.loc[s, d.index]
            p.append(smf.ols(form, data=d).fit().pvalues[target])
        print(f"   adjusted for {adj:22s} -> {nsig(p):3d}")
    # stratified: within metformin users and non-users separately
    for mv in [1, 0]:
        sub = m.index[m.metformin == mv]
        aa = [i for i in sub if m.loc[i, target] == 1]
        bb = [i for i in sub if m.loc[i, target] == 0]
        if min(len(aa), len(bb)) >= 20:
            print(f"   within metformin={mv} (n={len(aa)}/{len(bb)}) -> {nsig([mannwhitneyu(X.loc[s, aa], X.loc[s, bb]).pvalue for s in keep]):3d}")
