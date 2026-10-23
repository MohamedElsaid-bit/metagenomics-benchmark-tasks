"""Task 3 (mgx_crc_partial_contamination): stool profiles with partial oral carry over, and the answer key.

Unlike task 2 (whole mislabelled plaque samples), each affected sample here is a genuine stool
profile mixed with 18 to 32 percent oral material, so gut taxa remain present and the samples
still look like stool. Detection requires quantifying oral-specific taxa, not outlier spotting.
"""
import sys
import json
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA, OUT = sys.argv[1], sys.argv[2]
SEED = 20261022
N_CRC_SPIKE, N_CTL_SPIKE = 20, 1
F_LO, F_HI = 0.25, 0.40

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

zel, gh = load("ZellerG_2014"), load("GhensiP_2019")
zm = meta.loc[zel.columns]
stool = zel[zm.index[zm.study_condition.isin(["CRC", "control"])]]
stool = stool / stool.sum() * 100
zm = zm.loc[stool.columns]
plaque = gh / gh.sum() * 100
# strong plaque donors only, so the mixed fraction is reliably oral
donor_score = plaque.reindex(ORAL_STRICT).fillna(0).sum()
donors = donor_score[donor_score >= 40].index.tolist()
print(f"eligible plaque donors: {len(donors)} of {plaque.shape[1]}")

rng = np.random.default_rng(SEED)
crc_ids = list(zm.index[zm.study_condition == "CRC"])
ctl_ids = list(zm.index[zm.study_condition == "control"])
# Contaminate CRC samples that do NOT carry the genuine CRC markers, so removing them
# leaves the real signal intact and only the artefactual signal disappears.
TRUE_MARKERS = ["Fusobacterium_nucleatum", "Peptostreptococcus_stomatis"]
carrier = (stool.reindex(TRUE_MARKERS).fillna(0) > 0).any(axis=0)
pool = [s for s in crc_ids if not carrier[s]]
print(f"CRC non-carriers available for contamination: {len(pool)} of {len(crc_ids)}")
hit_crc = list(rng.choice(sorted(pool), N_CRC_SPIKE, replace=False))
hit_ctl = list(rng.choice(sorted(ctl_ids), N_CTL_SPIKE, replace=False))
hits = hit_crc + hit_ctl
donor_pick = list(rng.choice(sorted(donors), len(hits), replace=False))
fracs = rng.uniform(F_LO, F_HI, len(hits))

allsp = stool.index.union(plaque.index)
mixed = stool.reindex(allsp).fillna(0.0)
plq = plaque.reindex(allsp).fillna(0.0)
truth = {}
for sid, d, f in zip(hits, donor_pick, fracs):
    mixed[sid] = (1 - f) * mixed[sid] + f * plq[d]
    truth[sid] = {"donor": d, "fraction": round(float(f), 4)}
mixed = mixed / mixed.sum() * 100
mixed = mixed[mixed.sum(axis=1) > 0]

# anonymize and shuffle
order = rng.permutation(mixed.columns)
idmap = {old: f"P{i + 1:03d}" for i, old in enumerate(order)}
mixed = mixed[order].rename(columns=idmap)
md = pd.DataFrame({
    "group": ["CRC" if s in crc_ids else "control" for s in order],
    "age": pd.to_numeric(zm.loc[order, "age"], errors="coerce").to_numpy(dtype=float),
    "gender": zm.loc[order, "gender"].astype(object).to_numpy(),
}, index=[idmap[s] for s in order])
md.index.name = "sample_id"
mixed.index.name = "species"
mixed.to_csv(f"{OUT}/abundance_species.tsv", sep="\t", float_format="%.5f")
md.to_csv(f"{OUT}/sample_metadata.tsv", sep="\t")

hit_anon = sorted(idmap[s] for s in hits)
private = {idmap[s]: {**truth[s], "group": "CRC" if s in crc_ids else "control"} for s in hits}
json.dump({"contaminated": private, "id_map": {idmap[s]: s for s in order}},
          open(f"{OUT}/../truth_PRIVATE.json", "w"), indent=2)
print(f"table {mixed.shape}; contaminated {len(hit_anon)}: {hit_anon}")

# ---- separability ----
oral = mixed.reindex([s for s in ORAL_STRICT if s in mixed.index]).fillna(0).sum()
clean_ids = [c for c in mixed.columns if c not in hit_anon]
print(f"\noral-marker score: contaminated min {oral[hit_anon].min():.2f} max {oral[hit_anon].max():.2f} | "
      f"clean max {oral[clean_ids].max():.2f} p99 {oral[clean_ids].quantile(0.99):.2f}")
print(f"  {'CLEAN GAP' if oral[hit_anon].min() > oral[clean_ids].max() else 'OVERLAP'}")
gut = mixed.reindex([s for s in mixed.index if s.startswith(("Bacteroides", "Faecalibacterium", "Roseburia", "Alistipes", "Blautia"))]).sum()
print(f"gut-core score: contaminated {gut[hit_anon].min():.1f}-{gut[hit_anon].max():.1f} | "
      f"clean {gut[clean_ids].min():.1f}-{gut[clean_ids].max():.1f}  (overlapping by design)")

# ---- differential abundance with and without the contaminated samples ----
MARK = ["Fusobacterium_nucleatum", "Peptostreptococcus_stomatis", "Parvimonas_micra",
        "Gemella_morbillorum", "Faecalibacterium_prausnitzii", "Streptococcus_sanguinis",
        "Rothia_dentocariosa", "Porphyromonas_gingivalis"]


def da(crc, ctl, label):
    cols = list(crc) + list(ctl)
    keep = mixed.index[(mixed[cols] > 0).mean(axis=1) >= 0.10]
    p = pd.Series([mannwhitneyu(mixed.loc[s, crc], mixed.loc[s, ctl]).pvalue for s in keep], index=keep)
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
    up = mixed.loc[keep, crc].mean(axis=1) > mixed.loc[keep, ctl].mean(axis=1)
    print(f"\n{label}: tested {len(keep)} sig {int((q < 0.05).sum())}")
    out = {}
    for s in MARK:
        if s in q.index:
            call = "not_significant" if q[s] >= 0.05 else ("enriched_in_CRC" if up[s] else "depleted_in_CRC")
            print(f"   {s:30s} q={q[s]:.2e} {call}")
        else:
            call = "not_significant"
            print(f"   {s:30s} below prevalence filter -> not_significant")
        out[s] = call
    return out


crc_all = list(md.index[md.group == "CRC"])
ctl_all = list(md.index[md.group == "control"])
dirty = da(crc_all, ctl_all, "NAIVE (all samples, contamination left in)")
clean = da([c for c in crc_all if c not in hit_anon], [c for c in ctl_all if c not in hit_anon],
           "CORRECT (contaminated samples removed)")
print("\nspecies whose call flips among MARK:", [s for s in MARK if dirty[s] != clean[s]])


def full(crc, ctl):
    cols = list(crc) + list(ctl)
    keep = mixed.index[(mixed[cols] > 0).mean(axis=1) >= 0.10]
    p = pd.Series([mannwhitneyu(mixed.loc[s, crc], mixed.loc[s, ctl]).pvalue for s in keep], index=keep)
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
    up = mixed.loc[keep, crc].mean(axis=1) > mixed.loc[keep, ctl].mean(axis=1)
    return q, up


qd, ud = full(crc_all, ctl_all)
qc, uc = full([c for c in crc_all if c not in hit_anon], [c for c in ctl_all if c not in hit_anon])
rows = []
for s in qd.index.union(qc.index):
    cd = "ns" if s not in qd.index or qd[s] >= 0.05 else ("up" if ud[s] else "down")
    cc = "ns" if s not in qc.index or qc[s] >= 0.05 else ("up" if uc[s] else "down")
    if cd != cc:
        rows.append((s, cd, cc, qd.get(s, np.nan), qc.get(s, np.nan),
                     float(mixed.loc[s, hit_anon].mean()) if s in mixed.index else np.nan))
flips = pd.DataFrame(rows, columns=["species", "naive", "clean", "q_naive", "q_clean", "mean_in_contaminated"])
print(f"\nALL species whose call flips: {len(flips)}")
print(flips.sort_values("q_naive").to_string(float_format=lambda x: f"{x:.3g}"))
flips.to_csv(f"{OUT}/../flips.csv", index=False)
