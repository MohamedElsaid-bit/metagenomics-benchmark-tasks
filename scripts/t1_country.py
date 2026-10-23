"""Task 1: does the metformin effect survive stratification / adjustment by recruitment country?"""
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
pc = sp[sp > 0].min().min() / 2
logsp = np.log(sp + pc)
clr = logsp - logsp.mean()

tok = meta.treatment.astype(str).str.split(";")
ANTIDIAB = {"antidiab", "metformin", "insulin", "bolus", "dppiv", "glp_1", "su", "sglt2", "glitazone", "acarbose"}
met = tok.apply(lambda t: "metformin" in t)
anyd = tok.apply(lambda t: bool(ANTIDIAB & set(t)))
dis = meta.disease.astype(str)
t2d = dis.str.contains("T2D")
healthy = dis == "healthy"
SPECIES = ["Escherichia_coli", "Intestinibacter_bartlettii", "Faecalibacterium_prausnitzii",
           "Roseburia_hominis", "Akkermansia_muciniphila"]


def wilcox(a, b, label):
    both = sp.loc[:, a | b]
    idx = both.index[(both > 0).mean(axis=1) >= 0.10]
    p = pd.Series([mannwhitneyu(sp.loc[s, a], sp.loc[s, b]).pvalue for s in idx], index=idx)
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=idx)
    det = "  ".join(f"{s.split('_')[0][:6]}:{q[s]:.1e}{'+' if sp.loc[s, a].mean() > sp.loc[s, b].mean() else '-'}"
                    for s in SPECIES if s in q)
    print(f"{label:38s} n={a.sum():3d}/{b.sum():3d} sig={int((q < 0.05).sum()):3d} | {det}")


print("== Wilcoxon within country")
for c in ["FRA", "DEU", "DNK"]:
    cm = meta.country == c
    wilcox(t2d & met & cm, t2d & ~met & cm, f"{c} T2D met vs no-met")
    wilcox(t2d & ~anyd & cm, healthy & cm, f"{c} T2D naive vs healthy")

print("\n== CLR OLS adjusted for country (+ age + BMI)")
for lab, a, b in [("T2D met vs no-met", t2d & met, t2d & ~met), ("T2D naive vs healthy", t2d & ~anyd, healthy)]:
    for form in ["y ~ grp + C(country)", "y ~ grp + C(country) + age + BMI"]:
        d = meta.loc[a | b, ["age", "BMI", "country"]].copy()
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
        print(f"{lab:22s} {form:34s} n={len(d):4d} sig={int((q < 0.05).sum()):3d} | {det}")
