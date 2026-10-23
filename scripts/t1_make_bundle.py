"""Build the files given to the agent for Task 1 (mgx_t2d_drug_confounding).

Source: curatedMetagenomicData 2022-10-19 release, study MetaCardis_2020_a
(ExperimentHub EH7764 relative_abundance; sampleMetadata). Identifiers that point
straight at the source publication are removed so the answer cannot be looked up.
"""
import sys
import numpy as np
import pandas as pd
import pyreadr

DATA, OUT = sys.argv[1], sys.argv[2]

meta = list(pyreadr.read_r(f"{DATA}/sampleMetadata.rda").values())[0]
meta = meta[meta.study_name == "MetaCardis_2020_a"].set_index("sample_id")
ab = list(pyreadr.read_r(f"{DATA}/MetaCardis_2020_a.relative_abundance.rda").values())[0]
if "rowname" in ab.columns:
    ab = ab.set_index("rowname")

# Replace cohort-specific sample IDs (e.g. M0x10MCx1134) with neutral ones; same random order in both files
rng = np.random.default_rng(20261022)
order = rng.permutation(ab.columns)
idmap = {old: f"S{i + 1:04d}" for i, old in enumerate(order)}
pd.Series(idmap, name="anon_id").rename_axis("original_id").to_csv(f"{OUT}/../id_map_PRIVATE.tsv", sep="\t")
ab = ab[order].rename(columns=idmap)
meta = meta.rename(index=idmap)

# Full MetaPhlAn3 table (all ranks, as the profiler emits it); agent decides how to use it
ab.index.name = "clade_name"
ab.to_csv(f"{OUT}/abundance_metaphlan3.tsv", sep="\t", float_format="%.5f")

keep = ["disease", "age", "gender", "BMI", "country", "smoker", "antibiotics_current_use",
        "treatment", "hba1c", "ldl", "triglycerides", "hscrp", "bristol_score", "number_reads"]
m = meta.loc[ab.columns, keep]
m.index.name = "sample_id"
m.to_csv(f"{OUT}/sample_metadata.tsv", sep="\t")
print(ab.shape, m.shape)
