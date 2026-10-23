"""Task 3 design check: what oral-taxa level occurs naturally in stool, and what mixing fraction is detectable?"""
import sys
import numpy as np
import pandas as pd
import pyreadr

DATA = sys.argv[1]
meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0].set_index("sample_id")


def load(study):
    name, ab = next(iter(pyreadr.read_r(f"{DATA}/{study}.relative_abundance.rda").items()))
    assert study in name
    if "rowname" in ab.columns:
        ab = ab.set_index("rowname")
    if ab.index.str.contains(r"\|s__").any():
        ab = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")]
        ab.index = ab.index.str.split("|s__", regex=False).str[-1]
    return ab.groupby(level=0).sum()


zel, gh = load("ZellerG_2014"), load("GhensiP_2019")
zm = meta.loc[zel.columns]
stool = zel[zm.index[zm.study_condition.isin(["CRC", "control"])]]
stool = stool / stool.sum() * 100
plaque = gh / gh.sum() * 100

# Strictly oral-exclusive species: essentially never resident in the distal gut
ORAL_STRICT = ["Rothia_dentocariosa", "Rothia_aeria", "Rothia_mucilaginosa", "Corynebacterium_matruchotii",
               "Lautropia_mirabilis", "Capnocytophaga_gingivalis", "Capnocytophaga_sputigena",
               "Capnocytophaga_leadbetteri", "Capnocytophaga_granulosa", "Aggregatibacter_segnis",
               "Aggregatibacter_aphrophilus", "Actinomyces_naeslundii", "Actinomyces_oris",
               "Actinomyces_johnsonii", "Streptococcus_gordonii", "Streptococcus_sanguinis",
               "Streptococcus_cristatus", "Streptococcus_oralis", "Porphyromonas_gingivalis",
               "Tannerella_forsythia", "Treponema_denticola", "Neisseria_elongata", "Neisseria_sicca",
               "Neisseria_mucosa", "Neisseria_flavescens", "Kingella_oralis", "Prevotella_denticola",
               "Prevotella_nigrescens", "Prevotella_oris", "Selenomonas_sputigena", "Atopobium_rimae",
               "Parascardovia_denticolens", "Campylobacter_gracilis", "Eikenella_corrodens"]
present = [s for s in ORAL_STRICT if s in stool.index or s in plaque.index]
print(f"oral marker species found: {len(present)}/{len(ORAL_STRICT)}")


def score(tab):
    return tab.reindex(present).fillna(0).sum()


s_nat, s_plq = score(stool), score(plaque)
print("\nnatural oral-marker score in genuine stool (n=%d):" % stool.shape[1])
print(f"  median {s_nat.median():.3f}  p90 {s_nat.quantile(0.90):.3f}  p99 {s_nat.quantile(0.99):.3f}  max {s_nat.max():.3f}")
print(f"  top 6 values: {np.round(sorted(s_nat)[-6:], 2)}")
print(f"\nplaque donor oral-marker score (n=%d): median {s_plq.median():.1f}  min {s_plq.min():.1f}" % plaque.shape[1])

print("\nexpected oral-marker score after mixing stool with plaque at fraction f:")
for f in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    exp_lo = (1 - f) * s_nat.median() + f * s_plq.quantile(0.10)
    exp_med = (1 - f) * s_nat.median() + f * s_plq.median()
    print(f"  f={f:.2f}: pessimistic donor {exp_lo:6.2f}, median donor {exp_med:6.2f}  "
          f"(natural max {s_nat.max():.2f}) -> {'SEPARABLE' if exp_lo > 2 * s_nat.max() else 'risky'}")

# how much does the F. prausnitzii / Bacteroides signal survive dilution?
for sp in ["Faecalibacterium_prausnitzii", "Parvimonas_micra", "Gemella_morbillorum", "Fusobacterium_nucleatum"]:
    n = stool.loc[sp].mean() if sp in stool.index else 0
    p = plaque.loc[sp].mean() if sp in plaque.index else 0
    print(f"{sp:30s} stool mean {n:7.3f}  plaque mean {p:7.3f}  at f=0.25 -> {0.75 * n + 0.25 * p:7.3f}")
