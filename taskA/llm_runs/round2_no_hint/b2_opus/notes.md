# Drug–microbiome deconfounding in T2D participants

Script: `work/analysis.py` (plus `work/extra.py` for the per-drug single-covariate profile).
Intermediate counts: `work/intermediate.json`.

## Data and setup

- `data/sample_metadata.tsv`: 1831 participants. Subset to `disease` containing `T2D` -> **792 participants**
  (includes mixed labels such as `CAD;T2D`, `HF;CAD;T2D`). No T2D sample had an empty profile, so all 792 were kept.
- `data/abundance_metaphlan3.tsv`: multi-rank MetaPhlAn3 table. Kept rows containing `s__` and not `t__`
  -> 713 species. Columns sum to 100 per sample, so species-level rows are a complete partition.
- **Tested feature set**: species with relative abundance > 0 in >= 10% of the 792 T2D participants
  -> **195 species**. Every p-value set is BH-corrected across all 195; significance = q < 0.05.
- Abundances arcsine-square-root transformed (`asin(sqrt(x/100))`). Every result was repeated with a
  rank-inverse-normal transform as a sensitivity check; conclusions identical.
- Drug exposure = exact token match after splitting `treatment` on `;`.
  Counts in T2D: metformin 410, statin 440, PPI 150, su 134, insulin 118 (also dppiv 263, asa 313).

## Models

1. **Univariate**: species ~ drug (drug tested on its own).
2. **Full/deconfounded**: species ~ all 5 drugs + every other drug token present in >= 5% of T2D
   participants (beta_blockers, asa, dppiv, at2_inhib, ace_inhib, thiazidique, ca2_cbl, antithombo,
   bolus, diuranse, kspdiur, glp_1, xantoxinh) + age, BMI, gender, smoking, country, current
   antibiotics, log10 read depth. Umbrella class tokens (`antidiab`, `antilipid`, `antihta`) were
   dropped because they are near-deterministic unions of the specific drugs. Missing age/BMI
   (199/204, essentially the Danish subset) were median-imputed with an accompanying missingness
   indicator so that country coverage was not lost; missing categorical values got their own level.

## Results: significant species (q < 0.05 of 195)

| drug | univariate | full model | label |
|---|---|---|---|
| metformin | 94 | **49** | independent |
| PPI | 12 | **10** | independent |
| sulfonylurea | 33 | **0** | not_independent |
| statin | 24 | **0** | not_independent |
| insulin | 2 | **1** | not_independent |

Sensitivity runs (5 drugs mutually adjusted only; 5 drugs + demographics; all drug tokens >= 2% +
demographics; rank-INT outcome) all give the same picture: metformin 49–83, PPI 9–14, and
sulfonylurea / statin / insulin collapse to 0–2. Metformin is clearly the **strongest independent**
drug, with roughly five times as many retained species as PPI.

Insulin never had a substantial signal to begin with (2 species univariate) and retains at most 1,
so it is labelled `not_independent` under either reading of the criterion.

## What explains the sulfonylurea signal

Adjusting sulfonylurea for one candidate at a time (species counts at q < 0.05):

| adjusted for | n_sig |
|---|---|
| nothing | 33 |
| **metformin** | **1** |
| PPI | 32 |
| insulin | 33 |
| statin | 31 |
| dpp4 inhibitor | 22 |
| aspirin | 19 |
| country | 11 |
| age | 8 |
| BMI | 8 |

Metformin alone removes essentially the entire signal (33 -> 1); nothing else comes close. The
reason is co-prescription: **86.6%** of sulfonylurea users also take metformin, versus 44.7% of
non-users. Sulfonylurea is therefore largely a marker of metformin exposure here, and the answer to
field 3 is `metformin`.

By contrast, the statin signal has no single explanatory variable — it is knocked out by metformin,
aspirin, country, age and BMI alike (each 24 -> 0), i.e. it reflects the general confounding
structure of the cohort rather than one partner drug.

## Overlap

Among the 440 T2D statin users, 180 also take metformin -> **0.409**. (Note this is *lower* than the
overall metformin rate of 410/792 = 0.518; statin and metformin are negatively associated in this
cohort, so the statin signal is not a disguised metformin effect in the simple direction — it is
attenuated by the broader covariate set.)
