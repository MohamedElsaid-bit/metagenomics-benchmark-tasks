"""Task 3: are the six graded species calls stable across plausible flag sets an agent might derive?"""
import sys
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

T3 = sys.argv[1]
ab = pd.read_csv(f"{T3}/data/abundance_species.tsv", sep="\t", index_col=0)
md = pd.read_csv(f"{T3}/data/sample_metadata.tsv", sep="\t", index_col=0)
truth = json.load(open(f"{T3}/truth_PRIVATE.json"))
TRUE = sorted(truth["contaminated"].keys())
GRADED = ["Fusobacterium_nucleatum", "Peptostreptococcus_stomatis", "Parvimonas_micra",
          "Rothia_dentocariosa", "Eubacterium_rectale", "Streptococcus_salivarius"]

STRICT = {"Rothia_dentocariosa", "Rothia_aeria", "Rothia_mucilaginosa", "Corynebacterium_matruchotii",
          "Lautropia_mirabilis", "Capnocytophaga_gingivalis", "Capnocytophaga_sputigena", "Capnocytophaga_leadbetteri",
          "Capnocytophaga_granulosa", "Aggregatibacter_segnis", "Aggregatibacter_aphrophilus", "Actinomyces_naeslundii",
          "Actinomyces_oris", "Actinomyces_johnsonii", "Streptococcus_gordonii", "Streptococcus_sanguinis",
          "Streptococcus_cristatus", "Streptococcus_oralis", "Porphyromonas_gingivalis", "Tannerella_forsythia",
          "Treponema_denticola", "Neisseria_elongata", "Neisseria_sicca", "Neisseria_mucosa", "Neisseria_flavescens",
          "Kingella_oralis", "Prevotella_denticola", "Prevotella_nigrescens", "Prevotella_oris", "Selenomonas_sputigena",
          "Atopobium_rimae", "Parascardovia_denticolens", "Campylobacter_gracilis", "Eikenella_corrodens"}
GENERA = ("Rothia_", "Actinomyces_", "Neisseria_", "Capnocytophaga_", "Corynebacterium_", "Lautropia_",
          "Aggregatibacter_", "Porphyromonas_", "Tannerella_", "Treponema_", "Kingella_", "Eikenella_",
          "Selenomonas_", "Atopobium_", "Veillonella_", "Streptococcus_", "Prevotella_", "Haemophilus_")

strict_score = ab.reindex([s for s in STRICT if s in ab.index]).sum()
genus_score = ab[[any(s.startswith(g) for g in GENERA) for s in ab.index]].sum()
rng = np.random.default_rng(7)


def calls(excl, method="wilcoxon", prev=0.10):
    crc = [c for c in md.index[md.group == "CRC"] if c not in excl]
    ctl = [c for c in md.index[md.group == "control"] if c not in excl]
    cols = crc + ctl
    sub = ab[cols]
    keep = sub.index[(sub > 0).mean(axis=1) >= prev]
    if method == "wilcoxon":
        p = [mannwhitneyu(ab.loc[s, crc], ab.loc[s, ctl]).pvalue for s in keep]
        up = [ab.loc[s, crc].mean() > ab.loc[s, ctl].mean() for s in keep]
    else:
        pc = sub[sub > 0].min().min() / 2
        lg = np.log(sub + pc)
        clr = lg.sub(lg.mean(axis=0), axis=1)
        r = [ttest_ind(clr.loc[s, crc], clr.loc[s, ctl], equal_var=False) for s in keep]
        p = [x.pvalue for x in r]
        up = [x.statistic > 0 for x in r]
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
    up = pd.Series(up, index=keep)
    out = {}
    for s in GRADED:
        out[s] = ("not_significant" if s not in q.index or q[s] >= 0.05
                  else ("enriched_in_CRC" if up[s] else "depleted_in_CRC"))
    return out


scenarios = {"truth (21)": set(TRUE)}
for thr in [5, 8, 10, 12, 15, 20]:
    scenarios[f"strict score > {thr}%"] = set(strict_score.index[strict_score > thr])
for thr in [10, 15, 20, 25]:
    scenarios[f"genus score > {thr}%"] = set(genus_score.index[genus_score > thr])
# imperfect sets: miss some true ones, add some false ones
for miss, false in [(3, 0), (7, 0), (0, 5), (7, 5), (10, 3)]:
    s = set(rng.choice(TRUE, len(TRUE) - miss, replace=False))
    pool = [c for c in ab.columns if c not in TRUE]
    s |= set(rng.choice(pool, false, replace=False))
    scenarios[f"miss {miss}, false {false} (n={len(s)})"] = s
scenarios["none excluded"] = set()

ref = {s: v for s, v in zip(GRADED, ["enriched_in_CRC", "enriched_in_CRC", "enriched_in_CRC",
                                     "not_significant", "not_significant", "not_significant"])}
print(f"{'scenario':28s} {'n':>3s}  {'TP':>3s} {'FP':>3s}  wilcoxon            clr_t")
for name, excl in scenarios.items():
    tp, fp = len(excl & set(TRUE)), len(excl - set(TRUE))
    cw, cc = calls(excl, "wilcoxon"), calls(excl, "clr_t")
    okw = all(cw[s] == ref[s] for s in GRADED)
    okc = all(cc[s] == ref[s] for s in GRADED)
    dw = "" if okw else " wrong: " + ",".join(s.split("_")[0][:5] for s in GRADED if cw[s] != ref[s])
    dc = "" if okc else " wrong: " + ",".join(s.split("_")[0][:5] for s in GRADED if cc[s] != ref[s])
    print(f"{name:28s} {len(excl):3d}  {tp:3d} {fp:3d}  {'MATCH' if okw else 'DIFF ':5s}{dw:14s} {'MATCH' if okc else 'DIFF ':5s}{dc}")

print("\nflag-set size produced by each threshold rule:")
for thr in [5, 8, 10, 12, 15, 20]:
    s = set(strict_score.index[strict_score > thr])
    print(f"  strict > {thr:2d}%: n={len(s):3d} (TP {len(s & set(TRUE))}, FP {len(s - set(TRUE))})")
