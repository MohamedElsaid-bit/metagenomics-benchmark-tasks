"""Task 1: is metformin the dominant drug effect among T2D patients, or do other drugs rival it?"""
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
t2d = meta.disease.astype(str).str.contains("T2D")
drugs = ["metformin", "statin", "PPI", "dppiv", "su", "insulin", "beta_blockers", "asa", "antihta", "ace_inhib", "at2_inhib", "ca2_cbl", "thiazidique"]

print("Within T2D (n=%d): users / non-users, species sig (Wilcoxon, prev>=10%%, BH q<0.05) and adjusted for metformin+age+BMI" % t2d.sum())
rows = []
for d in drugs:
    use = tok.apply(lambda t: d in t)
    a, b = t2d & use, t2d & ~use
    both = sp.loc[:, t2d]
    idx = both.index[(both > 0).mean(axis=1) >= 0.10]
    p = [mannwhitneyu(sp.loc[s, a], sp.loc[s, b]).pvalue for s in idx]
    n_sig = int((multipletests(p, method="fdr_bh")[1] < 0.05).sum())
    # adjusted: CLR ~ drug + metformin + age + BMI (metformin term dropped when d is metformin)
    dd = meta.loc[t2d, ["age", "BMI"]].copy()
    dd["drug"] = use[t2d].astype(int)
    dd["met"] = tok[t2d].apply(lambda t: "metformin" in t).astype(int)
    dd = dd.dropna()
    form = "y ~ drug + age + BMI" + ("" if d == "metformin" else " + met")
    pa = []
    for s in idx:
        dd["y"] = clr.loc[s, dd.index]
        pa.append(smf.ols(form, data=dd).fit().pvalues["drug"])
    n_adj = int((multipletests(pa, method="fdr_bh")[1] < 0.05).sum())
    rows.append((d, int(a.sum()), int(b.sum()), n_sig, n_adj))
    print(f"  {d:14s} users={a.sum():4d} non={b.sum():4d} sig={n_sig:3d} adj_sig={n_adj:3d}", flush=True)

# Akkermansia: prior says metformin raises it (mouse studies). What does this cohort say?
met = tok.apply(lambda t: "metformin" in t)
s = "Akkermansia_muciniphila"
for lab, a, b in [("T2D met vs T2D no-met", t2d & met, t2d & ~met)]:
    p = mannwhitneyu(sp.loc[s, a], sp.loc[s, b]).pvalue
    print(f"\n{s} {lab}: p={p:.3g} mean {sp.loc[s, a].mean():.3f} vs {sp.loc[s, b].mean():.3f} "
          f"prev {(sp.loc[s, a] > 0).mean():.3f} vs {(sp.loc[s, b] > 0).mean():.3f}")
