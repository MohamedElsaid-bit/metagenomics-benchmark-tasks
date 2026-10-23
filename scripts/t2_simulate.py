"""Task 2 design: how many oral plaque samples, relabeled as CRC, are needed to create false CRC markers?

Stool cohort: ZellerG_2014 (CRC vs control, adenomas dropped).
Oral source: HMP_2012 supra/subgingival plaque + GhensiP_2019 + ShiB_2015 subgingival plaque.
For each injection size and random seed, compare species calls (CRC vs control) with and without the
injected samples, under several reasonable methods.
"""
import sys
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

DATA = sys.argv[1]
meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0].set_index("sample_id")


def load(study):
    name, ab = next(iter(pyreadr.read_r(f"{DATA}/{study}.relative_abundance.rda").items()))
    if "rowname" in ab.columns:
        ab = ab.set_index("rowname")
    if ab.index.str.contains(r"\|s__").any():
        ab = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")]
        ab.index = ab.index.str.split("|s__", regex=False).str[-1]
    return ab.groupby(level=0).sum()


z = load("ZellerG_2014")
zm = meta.loc[z.columns]
crc = list(zm.index[zm.study_condition == "CRC"])
ctl = list(zm.index[zm.study_condition == "control"])
oral = pd.concat([load(s) for s in ["HMP_2012", "GhensiP_2019", "ShiB_2015"]], axis=1).fillna(0)
om = meta.loc[oral.columns]
pool = list(om.index[om.body_subsite.isin(["subgingival_plaque", "supragingival_plaque"])])
print(f"stool CRC={len(crc)} control={len(ctl)}; oral plaque pool={len(pool)}")

ALL = pd.concat([z, oral[pool]], axis=1).fillna(0)
ALL = ALL / ALL.sum() * 100
TRUE = ["Fusobacterium_nucleatum", "Peptostreptococcus_stomatis", "Dialister_pneumosintes",
        "Porphyromonas_asaccharolytica", "Parvimonas_micra", "Gemella_morbillorum"]
DECOY = ["Rothia_dentocariosa", "Haemophilus_parainfluenzae", "Streptococcus_sanguinis", "Actinomyces_naeslundii",
         "Porphyromonas_gingivalis", "Corynebacterium_matruchotii", "Tannerella_forsythia", "Streptococcus_mitis",
         "Veillonella_parvula", "Prevotella_intermedia", "Treponema_denticola", "Neisseria_subflava"]


def calls(case, control, method, prev):
    X = ALL[case + control]
    keep = X.index[(X > 0).mean(axis=1) >= prev]
    if method == "wilcoxon":
        p = [mannwhitneyu(ALL.loc[s, case], ALL.loc[s, control]).pvalue for s in keep]
        d = [np.sign(ALL.loc[s, case].mean() - ALL.loc[s, control].mean()) for s in keep]
    else:  # CLR Welch t
        pc = X[X > 0].min().min() / 2
        L = np.log(X + pc); C = L - L.mean()
        r = [ttest_ind(C.loc[s, case], C.loc[s, control], equal_var=False) for s in keep]
        p = [x.pvalue for x in r]; d = [np.sign(x.statistic) for x in r]
    q = multipletests(p, method="fdr_bh")[1]
    out = pd.Series("not_tested", index=TRUE + DECOY)
    for s, qq, dd in zip(keep, q, d):
        if s in out.index:
            out[s] = ("up" if dd > 0 else "down") if qq < 0.05 else "ns"
    return out, int((q < 0.05).sum())


METHODS = [("wilcox