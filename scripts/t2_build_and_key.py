"""Task 2 (mgx_crc_sample_provenance): build the spiked cohort and compute the answer key.

Base cohort: ZellerG_2014 CRC stool (curatedMetagenomicData, release 2021-03-31).
Spike-in: subgingival plaque metagenomes from GhensiP_2019, relabelled as CRC stool samples.
The spiked samples inflate exactly the oral taxa that the CRC literature reports as markers.
"""
import sys
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA, OUT = sys.argv[1], sys.argv[2]
N_SPIKE = 12
SEED = 20261022

meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0].set_index("sample_id")


def load(study):
    name, ab = next(iter(pyreadr.read_r(f"{DATA}/{study}.relative_abundance.rda").items()))
    assert study in name
    if "rowname" in ab.columns:
        ab = ab.set_index("rowname")
    if ab.index.str.contains(r"\|s__").any():  # full lineage format
        ab = ab[ab.index.str.contains(r"\|s__") & ~ab.index.str.contains(r"\|t__")]
        ab.index = ab.index.str.split("|s__", regex=False).str[-1]
    return ab.groupby(level=0).sum()


zel, gh = load("ZellerG_2014"), load("GhensiP_2019")
zm = meta.loc[zel.columns]
base = zm.index[zm.study_condition.isin(["CRC", "control"])]
zel = zel[base]
zm = zm.loc[base]

rng = np.random.default_rng(SEED)
spike_src = list(rng.choice(sorted(gh.columns), N_SPIKE, replace=False))
spike = gh[spike_src]

ab = pd.concat([zel, spike], axis=1).fillna(0.0)
ab = ab / ab.sum() * 100

# metadata: spiked samples are presented as CRC cases, ages/sex drawn from the real CRC range
labels = pd.Series("control", index=ab.columns)
labels[zm.index[zm.study_condition == "CRC"]] = "CRC"
labels[spike_src] = "CRC"
age = pd.Series(np.nan, index=ab.columns, dtype=float)
age.loc[zm.index] = pd.to_numeric(zm.age, errors="coerce").to_numpy(dtype=float)
gender = pd.Series(pd.NA, index=ab.columns, dtype=object)
gender.loc[zm.index] = zm.gender.astype(object).to_numpy()
real_crc_age = pd.to_numeric(zm.loc[zm.study_condition == "CRC", "age"], errors="coerce").dropna().to_numpy(dtype=float)
age.loc[spike_src] = rng.choice(real_crc_age, N_SPIKE)
gender.loc[spike_src] = rng.choice(["male", "female"], N_SPIKE)

order = rng.permutation(ab.columns)
idmap = {old: f"C{i + 1:03d}" for i, old in enumerate(order)}
ab = ab[order].rename(columns=idmap)
md = pd.DataFrame({"group": labels[order].values, "age": age[order].values,
                   "gender": gender[order].values}, index=[idmap[s] for s in order])
md.index.name = "sample_id"
md["group"] = pd.Categorical(md["group"], ["control", "CRC"])

spike_anon = sorted(idmap[s] for s in spike_src)
pd.Series({idmap[s]: s for s in order}, name="source").rename_axis("anon_id").to_csv(f"{OUT}/../id_map_PRIVATE.tsv", sep="\t")
ab.index.name = "species"
ab.to_csv(f"{OUT}/abundance_species.tsv", sep="\t", float_format="%.5f")
md.to_csv(f"{OUT}/sample_metadata.tsv", sep="\t")
print("table", ab.shape, "spiked:", spike_anon)

# ---- answer key ----
MARKERS = ["Fusobacterium_nucleatum", "Parvimonas_micra", "Peptostreptococcus_stomatis",
           "Porphyromonas_asaccharolytica", "Porphyromonas_gingivalis", "Faecalibacterium_prausnitzii"]


def da(cols_crc, cols_ctl, label):
    sub = ab[list(cols_crc) + list(cols_ctl)]
    keep = sub.index[(sub > 0).mean(axis=1) >= 0.10]
    p = pd.Series([mannwhitneyu(ab.loc[s, cols_crc], ab.loc[s, cols_ctl]).pvalue for s in keep], index=keep)
    q = pd.Series(multipletests(p, method="fdr_bh")[1], index=keep)
    up = ab.loc[keep, cols_crc].mean(axis=1) > ab.loc[keep, cols_ctl].mean(axis=1)
    print(f"\n{label}: tested {len(keep)} sig {int((q < 0.05).sum())}")
    for s in MARKERS:
        if s in q:
            print(f"   {s:32s} q={q[s]:.2e} {'up' if up[s] else 'down'}")
        else:
            print(f"   {s:32s} filtered out (prevalence < 10%)")
    return q, up


crc_all = md.index[md.group == "CRC"]
ctl = md.index[md.group == "control"]
crc_clean = [c for c in crc_all if c not in spike_anon]
q_contam, up_contam = da(crc_all, ctl, "WITH spiked samples (naive analysis)")
q_clean, up_clean = da(crc_clean, ctl, "AFTER removing spiked samples (correct analysis)")

key = {
    "n_non_stool": len(spike_anon),
    "non_stool_sample_ids": spike_anon,
    "n_significant_clean": int((q_clean < 0.05).sum()),
}
for s in MARKERS:
    key[f"clean_{s}"] = ("not_significant" if s not in q_clean.index or q_clean[s] >= 0.05
                         else ("enriched_in_CRC" if up_clean[s] else "depleted_in_CRC"))
print("\nkey:", {k: v for k, v in key.items() if k != "non_stool_sample_ids"})
pd.Series(key).to_json(f"{OUT}/../key_draft.json", indent=2)
