"""Task 2 exploration: CRC stool cohort vs oral plaque profiles (curatedMetagenomicData, MetaPhlAn3)."""
import sys
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA = sys.argv[1]
meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0].set_index("sample_id")


def load(study):
    d = pyreadr.read_r(f"{DATA}/{study}.relative_abundance.rda")
    name, ab = next(iter(d.items()))
    assert study in name, (study, name)
    if "rowname" in ab.columns:
        ab = ab.set_index("rowname")
    if ab.index.str.contains(r"\|s__").any():  # full lineage format
        ab = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")]
        ab.index = ab.index.str.split("|s__", regex=False).str[-1]
    ab = ab.groupby(level=0).sum()
    print(f"{study}: {ab.shape}, index example {ab.index[0]}")
    return ab / ab.sum() * 100


studies = ["ZellerG_2014", "WirbelJ_2018", "GhensiP_2019", "HMP_2012", "ShiB_2015"]
tabs = {s: load(s) for s in studies}

z = tabs["ZellerG_2014"]
zm = meta.loc[z.columns]
print(zm.study_condition.value_counts().to_dict(), zm.country.value_counts().to_dict())

# oral plaque samples
plaque_ids = {
    "GhensiP_2019": meta.index[(meta.study_name == "GhensiP_2019")],
    "HMP_2012": meta.index[(meta.study_name == "HMP_2012") & meta.body_subsite.isin(["subgingival_plaque", "supragingival_plaque"])],
    "ShiB_2015": meta.index[(meta.study_name == "ShiB_2015")],
}

MARKERS = ["Fusobacterium_nucleatum", "Parvimonas_micra", "Peptostreptococcus_stomatis", "Gemella_morbillorum",
           "Porphyromonas_asaccharolytica", "Solobacterium_moorei", "Porphyromonas_gingivalis", "Prevotella_intermedia",
           "Streptococcus_sanguinis", "Rothia_dentocariosa", "Haemophilus_parainfluenzae", "Neisseria_subflava",
           "Actinomyces_naeslundii", "Corynebacterium_matruchotii", "Treponema_denticola", "Tannerella_forsythia",
           "Bacteroides_vulgatus", "Faecalibacterium_prausnitzii"]


def stats(ab, ids, label):
    ids = [i for i in ids if i in ab.columns]
    sub = ab.reindex(index=MARKERS, columns=ids).fillna(0)
    print(f"\n{label} (n={len(ids)}): prevalence / mean %")
    for s in MARKERS:
        print(f"   {s:32s} prev={(sub.loc[s] > 0).mean():.2f} mean={sub.loc[s].mean():.3f}")


crc = zm.index[zm.study_condition == "CRC"]
ctl = zm.index[zm.study_condition == "control"]
stats(z, crc, "Zeller CRC stool")
stats(z, ctl, "Zeller control stool")
for s, ids in plaque_ids.items():
    stats(tabs[s], ids, f"{s} plaque")

# baseline CRC vs control differential abundance in Zeller (no contamination)
keep = z.index[(z[list(crc) + list(ctl)] > 0).mean(axis=1) >= 0.10]
p = pd.Series([mannwhitneyu(z.loc[s, crc], z.loc[s, ctl]).pvalue for s in keep], index=keep)
q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
up = (z.loc[keep, crc].mean(axis=1) > z.loc[keep, ctl].mean(axis=1))
print(f"\nZeller CRC vs control: tested {len(keep)}, sig q<0.05 {int((q < 0.05).sum())}")
print(pd.DataFrame({"q": q, "up_in_CRC": up}).sort_values("q").head(20).to_string())
