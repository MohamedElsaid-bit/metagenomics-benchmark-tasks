"""Task 2 robustness: which species calls are stable across reasonable DA methods, with and without the spiked samples?"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

T2 = sys.argv[1]
ab = pd.read_csv(f"{T2}/data/abundance_species.tsv", sep="\t", index_col=0)
md = pd.read_csv(f"{T2}/data/sample_metadata.tsv", sep="\t", index_col=0)
spike = sorted(pd.read_csv(f"{T2}/id_map_PRIVATE.tsv", sep="\t", index_col=0)
               .query("source.str.startswith('SID')").index)
if not spike:  # Ghensi IDs do not start with SID; detect by source prefix difference
    m = pd.read_csv(f"{T2}/id_map_PRIVATE.tsv", sep="\t", index_col=0)
    spike = sorted(m.index[~m.source.str.startswith("CCIS")])
print("spiked:", spike)

MARKERS = ["Fusobacterium_nucleatum", "Parvimonas_micra", "Peptostreptococcus_stomatis",
           "Porphyromonas_asaccharolytica", "Solobacterium_moorei", "Gemella_morbillorum",
           "Faecalibacterium_prausnitzii", "Eubacterium_rectale", "Bacteroides_vulgatus"]


def da(crc, ctl, method, prev):
    cols = list(crc) + list(ctl)
    sub = ab[cols]
    keep = sub.index[(sub > 0).mean(axis=1) >= prev]
    if method == "wilcoxon":
        p = [mannwhitneyu(ab.loc[s, crc], ab.loc[s, ctl]).pvalue for s in keep]
        up = [ab.loc[s, crc].median() > ab.loc[s, ctl].median()
              or (ab.loc[s, crc].median() == ab.loc[s, ctl].median() and ab.loc[s, crc].mean() > ab.loc[s, ctl].mean())
              for s in keep]
    else:  # clr_t
        pc = sub[sub > 0].min().min() / 2
        lg = np.log(sub + pc)
        clr = lg.sub(lg.mean(axis=0), axis=1)
        r = [ttest_ind(clr.loc[s, crc], clr.loc[s, ctl], equal_var=False) for s in keep]
        p = [x.pvalue for x in r]
        up = [x.statistic > 0 for x in r]
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
    return q, pd.Series(up, index=keep)


crc_all = list(md.index[md.group == "CRC"])
ctl = list(md.index[md.group == "control"])
crc_clean = [c for c in crc_all if c not in spike]

rows = []
for method in ["wilcoxon", "clr_t"]:
    for prev in [0.05, 0.10, 0.20]:
        qd, upd = da(crc_all, ctl, method, prev)
        qc, upc = da(crc_clean, ctl, method, prev)
        rows.append(("dirty", method, prev, int((qd < 0.05).sum()), qd, upd))
        rows.append(("clean", method, prev, int((qc < 0.05).sum()), qc, upc))
        print(f"{method:9s} prev={prev:.2f} sig: with spike {int((qd < 0.05).sum()):3d}  clean {int((qc < 0.05).sum()):3d}")

print("\nper species call stability (clean analysis):")
print(f"{'species':32s} " + "  ".join(f"{m[:4]}{int(p*100):02d}" for m in ['wilcoxon', 'clr_t'] for p in [0.05, 0.10, 0.20]))
for s in MARKERS:
    cells = []
    for which in ["clean", "dirty"]:
        for (w, method, prev, n, q, up) in rows:
            if w != which:
                continue
            if s not in q.index:
                cells.append("  --  ")
            else:
                cells.append(f"{'+' if up[s] else '-'}{q[s]:.0e}" if q[s] < 0.05 else f" ns   ")
    print(f"{s:32s} clean[{' '.join(cells[:6])}]  dirty[{' '.join(cells[6:])}]")
