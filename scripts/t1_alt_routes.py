"""Task 1 (prompt version b): do alternative analysis routes to 'T2D itself' give the same species calls?

Route A: untreated T2D vs healthy (stratification), adjusted for country/age/BMI.
Route B: all T2D vs healthy, CLR OLS with metformin as covariate (+ country/age/BMI).
Route C: route B plus every other glucose lowering drug class as covariates.
Route D: T2D without metformin (any other treatment allowed) vs healthy, adjusted.
"""
import sys
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.formula.api as smf
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
pc = sp[sp > 0].min().min() / 2
logsp = np.log(sp + pc)
clr = logsp - logsp.mean()

tok = meta.treatment.astype(str).str.split(";")
ANTIDIAB = {"antidiab", "metformin", "insulin", "bolus", "dppiv", "glp_1", "su", "sglt2", "glitazone", "acarbose"}
for d in ["metformin", "dppiv", "su", "insulin", "glp_1"]:
    meta[d] = tok.apply(lambda t, d=d: d in t).astype(int)
meta["anyd"] = tok.apply(lambda t: bool(ANTIDIAB & set(t)))
dis = meta.disease.astype(str)
t2d = dis.str.contains("T2D")
healthy = dis == "healthy"
SPECIES = ["Escherichia_coli", "Intestinibacter_bartlettii", "Faecalibacterium_prausnitzii",
           "Roseburia_hominis", "Akkermansia_muciniphila"]

routes = {
    "A untreated T2D vs healthy": (t2d & ~meta.anyd, healthy, "y ~ grp + C(country) + age + BMI"),
    "B all T2D, +metformin": (t2d, healthy, "y ~ grp + metformin + C(country) + age + BMI"),
    "C all T2D, +all antidiab": (t2d, healthy, "y ~ grp + metformin + dppiv + su + insulin + glp_1 + C(country) + age + BMI"),
    "D no-metformin T2D vs healthy": (t2d & (meta.metformin == 0), healthy, "y ~ grp + C(country) + age + BMI"),
    "B0 all T2D, +metformin, no age/BMI": (t2d, healthy, "y ~ grp + metformin + C(country)"),
}
for name, (a, b, form) in routes.items():
    d = meta.loc[a | b, ["age", "BMI", "country", "metformin", "dppiv", "su", "insulin", "glp_1"]].copy()
    d["grp"] = a[d.index].astype(int)
    if "age" in form:
        d = d.dropna()
    both = sp.loc[:, d.index]
    idx = both.index[(both > 0).mean(axis=1) >= 0.10]
    p, coef = [], []
    for s in idx:
        d["y"] = clr.loc[s, d.index]
        f = smf.ols(form, data=d).fit()
        p.append(f.pvalues["grp"]); coef.append(f.params["grp"])
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=idx)
    coef = pd.Series(coef, index=idx)
    det = "  ".join(f"{s.split('_')[0][:6]}:{q[s]:.1e}{'+' if coef[s] > 0 else '-'}" for s in SPECIES if s in q)
    print(f"{name:36s} n={len(d):4d} sig={int((q < 0.05).sum()):3d} | {det}", flush=True)
