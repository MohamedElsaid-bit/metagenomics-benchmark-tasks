"""Task 2: are the 12 planted plaque samples objectively separable from stool, or is the boundary fuzzy?"""
import sys
import numpy as np
import pandas as pd

T2 = sys.argv[1]
ab = pd.read_csv(f"{T2}/data/abundance_species.tsv", sep="\t", index_col=0)
md = pd.read_csv(f"{T2}/data/sample_metadata.tsv", sep="\t", index_col=0)
idmap = pd.read_csv(f"{T2}/id_map_PRIVATE.tsv", sep="\t", index_col=0)
spike = sorted(idmap.index[~idmap.source.str.startswith("CCIS")])

# Obligate gut anaerobe genera that are essentially absent from dental plaque
GUT = ["Bacteroides", "Faecalibacterium", "Roseburia", "Eubacterium", "Alistipes", "Parabacteroides",
       "Blautia", "Ruminococcus", "Bifidobacterium", "Prevotella_copri", "Akkermansia", "Dorea",
       "Coprococcus", "Anaerostipes", "Collinsella", "Odoribacter", "Barnesiella", "Clostridium"]
ORAL = ["Streptococcus", "Rothia", "Actinomyces", "Corynebacterium", "Neisseria", "Haemophilus",
        "Veillonella", "Capnocytophaga", "Leptotrichia", "Porphyromonas", "Tannerella", "Treponema",
        "Aggregatibacter", "Campylobacter", "Lautropia", "Cardiobacterium", "Kingella", "Selenomonas"]


def frac(prefixes):
    mask = ab.index.map(lambda s: any(s.startswith(p) for p in prefixes))
    return ab[mask].sum()


gut, oral = frac(GUT), frac(ORAL)
df = pd.DataFrame({"gut_pct": gut, "oral_pct": oral, "group": md.group,
                   "planted": [i in spike for i in ab.columns]})
df["bacteroides"] = ab.loc[[i for i in ab.index if i.startswith("Bacteroides")]].sum()
df["fprau"] = ab.loc["Faecalibacterium_prausnitzii"] if "Faecalibacterium_prausnitzii" in ab.index else 0.0

print("planted samples, sorted by gut_pct:")
print(df[df.planted].sort_values("gut_pct").to_string(float_format=lambda x: f"{x:8.3f}"))
print("\ngenuine stool samples with the LOWEST gut_pct:")
print(df[~df.planted].sort_values("gut_pct").head(8).to_string(float_format=lambda x: f"{x:8.3f}"))

print("\n--- separation summary ---")
print(f"planted gut_pct: min {df[df.planted].gut_pct.min():.2f}  max {df[df.planted].gut_pct.max():.2f}")
print(f"stool   gut_pct: min {df[~df.planted].gut_pct.min():.2f}  max {df[~df.planted].gut_pct.max():.2f}")
gap_lo, gap_hi = df[df.planted].gut_pct.max(), df[~df.planted].gut_pct.min()
print(f"gap between highest planted and lowest genuine stool: {gap_lo:.2f} vs {gap_hi:.2f} "
      f"({'CLEAN GAP' if gap_lo < gap_hi else 'OVERLAP'})")
for thr in [1, 5, 10, 20, 30]:
    sel = set(df.index[df.gut_pct < thr])
    print(f"  threshold gut_pct < {thr:2d}%: selects {len(sel):3d} samples, "
          f"{'exactly the planted set' if sel == set(spike) else f'{len(sel & set(spike))} planted + {len(sel - set(spike))} genuine'}")
print(f"\nBacteroides: planted max {df[df.planted].bacteroides.max():.3f}, stool min {df[~df.planted].bacteroides.min():.3f}")
print(f"F. prausnitzii: planted max {df[df.planted].fprau.max():.3f}, stool min {df[~df.planted].fprau.min():.3f}")
