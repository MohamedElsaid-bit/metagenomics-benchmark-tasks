# Drug–microbiome deconfounding in the T2D subset

## Data and cohort

- `data/abundance_metaphlan3.tsv`: 1100 clades x 1831 samples, MetaPhlAn3 percentages, all ranks stacked.
  Species rows = clade names with 7 pipe-separated fields ending in `|s__` (713 species total, no strain rows).
- `data/sample_metadata.tsv`: 1831 participants.
- T2D subset = `disease` contains `T2D` (`T2D`, `CAD;T2D`, `HF;CAD;T2D`, `T2D;HF`, `HF;T2D`) -> **792 participants**,
  all present as columns in the abundance table.

## Testing framework

- Species level only; prevalence filter: present (>0) in >= 10% of the 792 T2D participants -> **269 species tested**.
- Response: arcsine-square-root transform of the relative abundance expressed as a proportion
  (standard variance-stabilising transform for MetaPhlAn compositional percentages).
- Test: OLS linear model, Wald t-test on the drug's coefficient.
- Multiple testing: Benjamini-Hochberg across all 269 tested species, separately for each model; significant at q < 0.05.
- Sanity check: univariate Mann-Whitney U gave the same ordering of drugs (metformin 109, su 48, statin 38, PPI 13,
  insulin 11 significant species), so the conclusions are not an artefact of the linear model.

## Models

1. **Univariate**: `species ~ drug`.
2. **Full multivariate (deconfounded)**: `species ~ drug + all other co-prescribed drugs + demographics`.
   - Drug covariates: every `treatment` token present in >= 5% of T2D participants, excluding the three
     umbrella class labels (`antidiab`, `antilipid`, `antihta`) and the `no` placeholder, since those are
     aggregates of the specific drug tokens and would be collinear with them. That gives:
     metformin, PPI, su, statin, insulin, dppiv, asa, glp_1, bolus, beta_blockers, ace_inhib, at2_inhib,
     thiazidique, ca2_cbl, diuranse, kspdiur, antithombo, xantoxinh.
   - Demographics / technical: age, BMI, gender, country (DEU/DNK/FRA dummies), current antibiotic use, smoking.
   - `age` (199 missing), `BMI` (204), `gender` (199) and `smoker` (204) were median/mode imputed with an
     explicit missingness indicator so the 792 participants are retained. A complete-case re-run (n = 588,
     no imputation) gave the same conclusions (metformin 58, PPI 6, su 0, statin 0, insulin 1).
3. **Single-covariate scan**: `species ~ drug + one candidate confounder`, to identify which single variable
   accounts for an apparent association.

## Results

Significant species (q < 0.05) out of 269 tested:

| drug | n users | univariate | + all drugs only | full multivariate | complete case | label |
|---|---|---|---|---|---|---|
| metformin | 410 | 94 | 67 | **53** | 58 | independent |
| PPI | 150 | 12 | 12 | **10** | 6 | independent |
| sulfonylurea | 134 | 33 | 1 | **0** | 0 | not_independent |
| statin | 440 | 24 | 0 | **0** | 0 | not_independent |
| insulin | 118 | 2 (MW: 11) | 1 | **1** | 1 | not_independent |

The two surviving signatures are the biologically expected ones:

- **Metformin** (53 species): *Intestinibacter bartlettii* (strongly depleted, the canonical metformin marker),
  *Clostridium clostridioforme*, *Clostridium bolteae*, *Flavonifractor plautii*, *Faecalibacterium prausnitzii*,
  *Eubacterium eligens/ventriosum*, *Gemmiger formicilis*, Escherichia/Bifidobacterium shifts. Largest independent
  set by a wide margin -> `strongest_independent_drug = metformin`.
- **PPI** (10 species): an oral-to-gut translocation signature -- *Streptococcus salivarius*, *S. parasanguinis*,
  four *Veillonella* species, *Haemophilus parainfluenzae*, *Rothia mucilaginosa*, *Bifidobacterium dentium*.
  This is the textbook PPI effect and is essentially untouched by adjustment (12 -> 10).

### Why sulfonylurea and statin collapse

Co-prescription rates among T2D participants (row drug -> fraction also on column drug):

| | metformin | PPI | su | statin | insulin | dppiv | asa |
|---|---|---|---|---|---|---|---|
| sulfonylurea (n=134) | **0.866** | 0.164 | 1.0 | 0.522 | 0.231 | 0.507 | 0.269 |
| statin (n=440) | 0.409 | 0.225 | 0.159 | 1.0 | 0.148 | 0.282 | 0.580 |
| cohort baseline | 0.518 | 0.189 | 0.169 | 0.556 | 0.149 | 0.332 | 0.395 |

**Sulfonylurea** single-covariate scan (starting from 33 significant species):

| adjusted for | significant species remaining |
|---|---|
| metformin | **1** |
| dppiv | 22 |
| aspirin (`asa`) | 19 |
| country | 11 |
| age | 8 |
| BMI | 8 |
| PPI | 32 |
| statin | 31 |
| insulin | 33 |

Only metformin collapses the sulfonylurea signal (33 -> 1, and -> 0 in the full model). 86.6% of sulfonylurea
users are also on metformin versus a 51.8% cohort baseline, so sulfonylurea is essentially a proxy for metformin
here. -> `drug_explaining_sulfonylurea_signal = metformin`.

**Statin** likewise goes 24 -> 0 as soon as metformin (or the set of all antidiabetics) is included. The statin
signal was weak to begin with (best univariate q = 1.3e-3, nothing below q = 1e-3) and statin use is *negatively*
associated with metformin (40.9% vs 51.8% baseline), so the apparent statin effect is a diluted mirror image of
the metformin effect rather than a statin effect. Adjusting for gender alone also removes it, confirming the
signal is not robust.

**Insulin** barely reaches significance on its own (2 species by linear model, 11 by Mann-Whitney, best q = 0.032)
and retains a single species (*Veillonella dispar*) after adjustment -- not a substantial independent signature.

### statin_metformin_overlap

440 T2D participants take a statin; 180 of them also take metformin -> 180/440 = 0.4091 -> **0.409**.

## Files

- `work/analyze.py` -- main analysis (univariate, full multivariate, single-covariate scan, overlap).
- `work/diag.py` -- diagnostics: q-value distributions, co-prescription matrix, complete-case sensitivity,
  stepwise adjustment for statin/su/insulin, named significant species.
- `work/raw_results.json` -- machine-readable counts from `analyze.py`.
- `work/answer.json` -- final answer.
