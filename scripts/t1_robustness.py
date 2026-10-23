"""Task 1 robustness: are species-level calls (T2D-intrinsic vs metformin-driven) stable across reasonable methods?"""
import sys
import numpy as np
import pandas as pd
import pyreadr
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

DATA, OUT = sys.argv[1], sys.argv[2]

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

tx = meta.treatment.astype(str).str.split(";")
ANTIDIAB = {"antidiab", "metformin", "insulin", "bolus", "dppiv", "glp_1", "su", "sglt2", "glitazone", "acarbose"}
meta["metformin"] = tx.apply(lambda t: "metformin" in t)
meta["any_antidiab"] = tx.apply(lambda t: bool(ANTIDIAB & set(t)))
dis = meta.disease.astype(str)
t2d = dis.str.contains("T2D")
print("gender values:", meta.gender.value_counts(dropna=False).to_dict())

G = {
    "healthy": dis == "healthy",
    "T2D_naive": t2d & ~meta.any_antidiab,          # no glucose-lowering drug at all
    "T2D_met": t2d & meta.metformin,
    "T2D_nomet": t2d & ~meta.metformin,
    # non-diabetic, obese, not on any glucose-lowering drug: controls for adiposity
    "obese_nondiab": (~t2d) & (dis != "healthy") & (meta.BMI >= 30) & ~meta.any_antidiab,
}
for g, m in G.items():
    print(f"{g:14s} n={m.sum():4d} age={meta[m].age.mean():.1f} BMI={meta[m].BMI.mean():.1f}")

# CLR with pseudocount = half the smallest nonzero value
pc = sp[sp > 0].min().min() / 2
logsp = np.log(sp + pc)
clr = logsp - logsp.mean()


def tested(a, b, prev):
    both = sp.loc[:, G[a] | G[b]]
    return both.index[(both > 0).mean(axis=1) >= prev]


def run(a, b, method, prev=0.10):
    idx = tested(a, b, prev)
    ma, mb = G[a], G[b]
    p, eff = [], []
    for s in idx:
        if method == "wilcoxon":
            p.append(mannwhitneyu(sp.loc[s, ma], sp.loc[s, mb]).pvalue)
            eff.append(np.sign(sp.loc[s, ma].median() - sp.loc[s, mb].median()) or np.sign(sp.loc[s, ma].mean() - sp.loc[s, mb].mean()))
        elif method == "clr_t":
            r = ttest_ind(clr.loc[s, ma], clr.loc[s, mb], equal_var=False)
            p.append(r.pvalue)
            eff.append(np.sign(r.statistic))
        elif method == "clr_t_cc":  # unadjusted, but restricted to samples with age/BMI/sex (same n as clr_lm_adj)
            cc = meta[["age", "BMI", "gender"]].notna().all(axis=1)
            r = ttest_ind(clr.loc[s, ma & cc], clr.loc[s, mb & cc], equal_var=False)
            p.append(r.pvalue)
            eff.append(np.sign(r.statistic))
        elif method == "clr_lm_adj":
            d = meta.loc[ma | mb, ["age", "BMI", "gender"]].copy()
            d["y"] = clr.loc[s, d.index]
            d["grp"] = ma[d.index].astype(int)
            d = d.dropna()
            f = smf.ols("y ~ grp + age + BMI + C(gender)", data=d).fit()
            p.append(f.pvalues["grp"])
            eff.append(np.sign(f.params["grp"]))
    r = pd.DataFrame({"p": p, "dir": eff}, index=idx)
    r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r


comparisons = [("T2D_naive", "healthy"), ("T2D_naive", "obese_nondiab"), ("T2D_met", "T2D_nomet"), ("T2D_met", "healthy")]
methods = [("wilcoxon", 0.10), ("wilcoxon", 0.05), ("wilcoxon", 0.25), ("clr_t", 0.10), ("clr_t_cc", 0.10), ("clr_lm_adj", 0.10)]

summary, calls = [], {}
for a, b in comparisons:
    for meth, prev in methods:
        r = run(a, b, meth, prev)
        key = f"{a}_vs_{b}|{meth}|prev{prev}"
        calls[key] = r
        summary.append((f"{a} vs {b}", meth, prev, len(r), int((r.q < 0.05).sum())))
        print(f"{a:10s} vs {b:14s} {meth:11s} prev={prev:.2f} tested={len(r):3d} sig={int((r.q < 0.05).sum()):3d}", flush=True)

pd.DataFrame(summary, columns=["comparison", "method", "prev", "tested", "sig_q05"]).to_csv(f"{OUT}/t1_robustness_summary.csv", index=False)

# per-species q and direction across all method x comparison combos
long = []
for key, r in calls.items():
    comp, meth, prev = key.split("|")
    for s, row in r.iterrows():
        long.append((s, comp, meth, prev, row.q, row.dir))
long = pd.DataFrame(long, columns=["species", "comparison", "method", "prev", "q", "dir"])
long.to_csv(f"{OUT}/t1_robustness_species.csv", index=False)

# Robust category per species:
#  "metformin": sig (q<0.05) in T2D_met vs T2D_nomet under every method, and NOT sig (q>=0.20) in T2D_naive vs healthy under every method
#  "disease":   sig in T2D_naive vs healthy AND T2D_naive vs obese_nondiab under every method, same direction
piv = long.pivot_table(index="species", columns=["comparison", "method", "prev"], values="q")
dpiv = long.pivot_table(index="species", columns=["comparison", "method", "prev"], values="dir")


def allcols(comp):
    return [c for c in piv.columns if c[0] == comp]


out = []
for s in piv.index:
    met_q = piv.loc[s, allcols("T2D_met_vs_T2D_nomet")]
    nh_q = piv.loc[s, allcols("T2D_naive_vs_healthy")]
    no_q = piv.loc[s, allcols("T2D_naive_vs_obese_nondiab")]
    if met_q.notna().all() and (met_q < 0.05).all() and nh_q.notna().all() and (nh_q >= 0.20).all():
        cat = "metformin"
    elif nh_q.notna().all() and (nh_q < 0.05).all() and no_q.notna().all() and (no_q < 0.05).all() \
            and dpiv.loc[s, allcols("T2D_naive_vs_healthy") + allcols("T2D_naive_vs_obese_nondiab")].nunique() == 1:
        cat = "disease"
    else:
        cat = "ambiguous"
    out.append((s, cat, met_q.max(), nh_q.min(), no_q.min(),
                dpiv.loc[s, allcols("T2D_met_vs_T2D_nomet")].mean()))
cat = pd.DataFrame(out, columns=["species", "category", "max_q_met_vs_nomet", "min_q_naive_vs_healthy",
                                 "min_q_naive_vs_obese", "dir_met"]).set_index("species")
cat.to_csv(f"{OUT}/t1_species_categories.csv")
print("\nrobust categories:", cat.category.value_counts().to_dict())
for c in ["metformin", "disease"]:
    print(f"\n-- {c}")
    print(cat[cat.category == c].sort_values("max_q_met_vs_nomet").to_string(float_format=lambda x: f"{x:.2e}"))
print("\nmarker species:")
print(cat.loc[[s for s in ["Escherichia_coli", "Intestinibacter_bartlettii", "Roseburia_hominis",
                           "Faecalibacterium_prausnitzii", "Anaerostipes_hadrus", "Akkermansia_muciniphila"] if s in cat.index]]
      .to_string(float_format=lambda x: f"{x:.2e}"))
