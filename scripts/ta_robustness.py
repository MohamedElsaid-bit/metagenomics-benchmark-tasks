"""Task A (mgx_drug_coprescription_confounding): is each drug's microbiome association independent
of co-prescription, under every reasonable way of controlling for the other drugs?

Routes tested per drug, within T2D participants:
  R1 unadjusted Wilcoxon (users vs non-users)
  R2 CLR OLS adjusted for the other six drugs
  R3 R2 plus age, BMI, sex
  R4 R3 plus country
  R5 stratified: among metformin users only, users vs non-users of the drug
  R6 stratified: among metformin non-users only
"""
import sys
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA = sys.argv[1]
PREV = 0.10
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
DRUGS = {"metformin": "metformin", "dpp4_inhibitor": "dppiv", "proton_pump_inhibitor": "PPI",
         "sulfonylurea": "su", "statin": "statin", "aspirin": "asa", "insulin": "insulin"}
for name, t in DRUGS.items():
    meta[name] = tok.apply(lambda x, t=t: t in x).astype(int)
t2d = meta.disease.astype(str).str.contains("T2D")
m = meta[t2d].copy()
X = sp.loc[:, m.index]
pc = X[X > 0].min().min() / 2
lg = np.log(X + pc)
clr = lg.sub(lg.mean(axis=0), axis=1)
keep = X.index[(X > 0).mean(axis=1) >= PREV]
print(f"T2D participants {len(m)}, species tested {len(keep)}")
print("users per drug:", {d: int(m[d].sum()) for d in DRUGS})

# co-prescription overlap
print("\nco-prescription with metformin (fraction of that drug's users also on metformin):")
for d in DRUGS:
    if d == "metformin":
        continue
    u = m.index[m[d] == 1]
    print(f"  {d:22s} n={len(u):3d}  on metformin {m.loc[u, 'metformin'].mean():.3f}")


def nsig(pvals):
    return int((multipletests(pvals, method="fdr_bh")[1] < 0.05).sum())


def run(drug, route):
    others = [d for d in DRUGS if d != drug]
    if route == "R1":
        a, b = m.index[m[drug] == 1], m.index[m[drug] == 0]
        return nsig([mannwhitneyu(X.loc[s, a], X.loc[s, b]).pvalue for s in keep])
    if route in ("R5", "R6"):
        sub = m.index[m.metformin == (1 if route == "R5" else 0)]
        a, b = [i for i in sub if m.loc[i, drug] == 1], [i for i in sub if m.loc[i, drug] == 0]
        if min(len(a), len(b)) < 20:
            return None
        return nsig([mannwhitneyu(X.loc[s, a], X.loc[s, b]).pvalue for s in keep])
    cov = " + ".join(others)
    form = f"y ~ {drug} + {cov}"
    d = m[list(DRUGS)].copy()
    if route in ("R3", "R4"):
        d = d.join(m[["age", "BMI", "gender"]])
        form += " + age + BMI + C(gender)"
    if route == "R4":
        d = d.join(m[["country"]])
        form += " + C(country)"
    d = d.dropna()
    p = []
    for s in keep:
        d["y"] = clr.loc[s, d.index]
        p.append(smf.ols(form, data=d).fit().pvalues[drug])
    return nsig(p)


routes = ["R1", "R2", "R3", "R4", "R5", "R6"]
print(f"\n{'drug':22s} " + "".join(f"{r:>6s}" for r in routes) + "   label")
rows = {}
for drug in DRUGS:
    vals = [run(drug, r) for r in routes]
    rows[drug] = vals
    adj = [v for v, r in zip(vals, routes) if r != "R1" and v is not None]
    label = "independent" if min(adj) >= 5 else ("confounded" if max(adj) <= 2 else "AMBIGUOUS")
    print(f"{drug:22s} " + "".join(f"{'-' if v is None else v:>6}" for v in vals) + f"   {label}")
pd.DataFrame(rows, index=routes).T.to_csv("/tmp/ta_routes.csv")
