"""Task 3 robustness: stability of candidate graded species, and how detectable the contamination is."""
import sys
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

T3 = sys.argv[1]
ab = pd.read_csv(f"{T3}/data/abundance_species.tsv", sep="\t", index_col=0)
md = pd.read_csv(f"{T3}/data/sample_metadata.tsv", sep="\t", index_col=0)
truth = json.load(open(f"{T3}/truth_PRIVATE.json"))
hits = sorted(truth["contaminated"].keys())
clean_ids = [c for c in ab.columns if c not in hits]

CAND = ["Fusobacterium_nucleatum", "Peptostreptococcus_stomatis", "Parvimonas_micra", "Gemella_morbillorum",
        "Rothia_dentocariosa", "Streptococcus_oralis", "Eubacterium_rectale", "Streptococcus_salivarius",
        "Roseburia_intestinalis", "Ruminococcus_bromii", "Faecalibacterium_prausnitzii"]


def da(crc, ctl, method, prev):
    cols = list(crc) + list(ctl)
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
    return q, pd.Series(up, index=keep)


crc_all = list(md.index[md.group == "CRC"])
ctl_all = list(md.index[md.group == "control"])
crc_c = [c for c in crc_all if c not in hits]
ctl_c = [c for c in ctl_all if c not in hits]

combos = [(m, p) for m in ["wilcoxon", "clr_t"] for p in [0.05, 0.10, 0.20]]
print("CLEAN analysis calls per method (prev filter 5/10/20 pct, Wilcoxon then CLR t):")
print(f"{'species':32s} " + " ".join(f"{m[:4]}{int(p*100):02d}" for m, p in combos) + "   | naive(wilc10)")
res_clean, res_naive = {}, {}
for m, p in combos:
    res_clean[(m, p)] = da(crc_c, ctl_c, m, p)
    res_naive[(m, p)] = da(crc_all, ctl_all, m, p)
for s in CAND:
    cells = []
    for m, p in combos:
        q, up = res_clean[(m, p)]
        cells.append("  --  " if s not in q.index else (f"{'+' if up[s] else '-'}{q[s]:.0e}" if q[s] < 0.05 else "  ns  "))
    qn, un = res_naive[("wilcoxon", 0.10)]
    nv = "  --  " if s not in qn.index else (f"{'+' if un[s] else '-'}{qn[s]:.0e}" if qn[s] < 0.05 else "  ns  ")
    print(f"{s:32s} " + " ".join(cells) + f"   | {nv}")

print("\n--- detectability of the contaminated samples ---")
ORAL = [s for s in ab.index if s in {
    "Rothia_dentocariosa", "Rothia_aeria", "Rothia_mucilaginosa", "Corynebacterium_matruchotii",
    "Lautropia_mirabilis", "Capnocytophaga_gingivalis", "Capnocytophaga_sputigena", "Capnocytophaga_leadbetteri",
    "Capnocytophaga_granulosa", "Aggregatibacter_segnis", "Aggregatibacter_aphrophilus", "Actinomyces_naeslundii",
    "Actinomyces_oris", "Actinomyces_johnsonii", "Streptococcus_gordonii", "Streptococcus_sanguinis",
    "Streptococcus_cristatus", "Streptococcus_oralis", "Porphyromonas_gingivalis", "Tannerella_forsythia",
    "Treponema_denticola", "Neisseria_elongata", "Neisseria_sicca", "Neisseria_mucosa", "Neisseria_flavescens",
    "Kingella_oralis", "Prevotella_denticola", "Prevotella_nigrescens", "Prevotella_oris", "Selenomonas_sputigena",
    "Atopobium_rimae", "Parascardovia_denticolens", "Campylobacter_gracilis", "Eikenella_corrodens"}]
oral = ab.loc[ORAL].sum()
print(f"oral-marker score: contaminated {oral[hits].min():.2f}-{oral[hits].max():.2f} | clean {oral[clean_ids].min():.3f}-{oral[clean_ids].max():.2f}")
print(f"gap: {'CLEAN' if oral[hits].min() > oral[clean_ids].max() else 'OVERLAP'} "
      f"(ratio {oral[hits].min() / max(oral[clean_ids].max(), 1e-9):.2f}x)")

pc = ab[ab > 0].min().min() / 2
lg = np.log(ab + pc)
clr = lg.sub(lg.mean(axis=0), axis=1)
pcs = PCA(n_components=4).fit_transform(clr.T.values)
lab = np.array([c in hits for c in ab.columns])
print("\nPCA on CLR: separation of contaminated vs clean per component")
for i in range(4):
    a, b = pcs[lab, i], pcs[~lab, i]
    overlap = not (a.max() < b.min() or b.max() < a.min())
    print(f"  PC{i+1}: contaminated [{a.min():7.2f},{a.max():7.2f}] clean [{b.min():7.2f},{b.max():7.2f}] "
          f"{'OVERLAP' if overlap else 'SEPARATES'}  (AUC-like gap {abs(np.median(a) - np.median(b)):.2f})")
print(f"\ncontaminated fractions used: {sorted(round(v['fraction'], 2) for v in truth['contaminated'].values())}")
print("groups of contaminated:", pd.Series([v["group"] for v in truth["contaminated"].values()]).value_counts().to_dict())
