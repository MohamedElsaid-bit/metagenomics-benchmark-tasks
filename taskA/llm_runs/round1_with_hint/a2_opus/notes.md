# Drug–microbiome associations in the T2D subgroup: independent vs. co-prescription-confounded

## Data and cohort definition

- `data/abundance_metaphlan3.tsv`: 1831 samples, all ranks. Species rows = lines whose `clade_name`
  contains `s__` and not `t__` → **713 species**.
- `data/sample_metadata.tsv`: 1831 participants.
- T2D subgroup = `disease` contains `T2D` (`T2D`, `CAD;T2D`, `HF;CAD;T2D`, `T2D;HF`, `HF;T2D`) →
  **792 participants**.
- Drug exposure: `treatment` split on `;` and matched as **exact tokens** (important — `su` would
  otherwise substring-match `insulin`, `sglt2` etc.).

Exposure counts among the 792 T2D participants:

| drug | token | n users |
|---|---|---|
| metformin | `metformin` | 410 |
| proton pump inhibitor | `PPI` | 150 |
| sulfonylurea | `su` | 134 |
| statin | `statin` | 440 |
| insulin | `insulin` | 118 |

## Species filter and testing framework

- Prevalence filter: species present (abundance > 0) in ≥ 10% of the 792 T2D participants →
  **195 species tested**.
- Abundance transform: `log10(relative_abundance_fraction + 1e-6)`.
- BH FDR applied across the full set of 195 tested species, per drug; significance at **q < 0.05**.

Two models per drug:

1. **Marginal (unadjusted)** — the drug alone (Mann–Whitney U, and OLS with the single drug term).
2. **Mutually adjusted** — one OLS per species with **all five drug indicators entered jointly**;
   the p-value for each drug is its own partial coefficient. This is the test of whether the
   association survives once the other four medications are accounted for.

## Co-prescription structure (fraction of row-drug users also on column drug)

|            | metformin | PPI | SU | statin | insulin |
|---|---|---|---|---|---|
| metformin  | 1.000 | 0.144 | 0.283 | 0.439 | 0.202 |
| PPI        | 0.393 | 1.000 | 0.147 | 0.660 | 0.133 |
| sulfonylurea | **0.866** | 0.164 | 1.000 | 0.522 | 0.231 |
| statin     | **0.409** | 0.225 | 0.159 | 1.000 | 0.148 |
| insulin    | **0.703** | 0.169 | 0.263 | 0.551 | 1.000 |

Metformin dominates the prescription network: 87% of sulfonylurea users and 70% of insulin users
are also on metformin, and 41% of statin users are.

## Results — significant species (q < 0.05)

| drug | marginal (MWU) | marginal (OLS) | **5-drug adjusted** | +age/BMI/sex/country/abx | all-medication adjusted |
|---|---|---|---|---|---|
| metformin | 109 | 117 | **99** | 86 | 32 |
| proton pump inhibitor | 13 | 13 | **12** | 13 | 13 |
| sulfonylurea | 48 | 59 | **2** | 0 | 0 |
| statin | 38 | 39 | **7** | 0 | 0 |
| insulin | 11 | 8 | **0** | 2 | 0 |

### Robustness of the adjusted model to the abundance transform

| transform | metformin | PPI | SU | statin | insulin |
|---|---|---|---|---|---|
| log10 | 99 | 12 | 2 | 7 | 0 |
| rank (inverse-normalised) | 95 | 13 | 1 | **0** | 0 |
| arcsine-sqrt | 83 | 12 | 1 | **0** | 2 |
| CLR | 93 | 16 | 2 | 6 | 0 |

A reverse-direction check (logistic regression `drug ~ species + other four drugs`) gives the same
ordering: metformin 99, PPI 12, SU 2, statin 3, insulin 0.

### Stratified confirmation

Because metformin is the hub of the co-prescription network, each other drug was re-tested inside
the **metformin-free** T2D stratum (n = 382):

- sulfonylurea (n+ = 18): **0** significant species
- statin (n+ = 260 of 382 — ample power): **0** significant species
- insulin (n+ = 35): **0** significant species
- PPI (n+ = 91): **13** significant species — signal fully intact without metformin

And metformin re-tested within statin strata:

- metformin within statin users (n = 440): **90** significant species
- metformin within statin non-users (n = 352): **66** significant species

## Interpretation and labels

- **metformin — `independent`.** 117 marginal → 99 mutually adjusted species (85% retained), stable
  across all four transforms, and reproduced separately in statin users and statin non-users. By far
  the largest independent signal.
- **proton pump inhibitor — `independent`.** 13 → 12 species (essentially complete retention), and
  13 species in the metformin-free stratum. The surviving taxa are the canonical PPI signature of
  oral/upper-GI colonisers translocating to the gut: *Streptococcus salivarius*, *S. parasanguinis*,
  *S. mitis*, four *Veillonella* spp., *Rothia mucilaginosa*, *Haemophilus parainfluenzae*,
  *Bifidobacterium dentium*. Smaller than metformin's but genuinely its own effect.
- **sulfonylurea — `confounded_by_co_prescription`.** 59 marginal → 2 adjusted (97% lost), 0 with
  covariates, 0 in the metformin-free stratum. 86.6% of SU users are on metformin, so the marginal
  signal is a metformin shadow.
- **statin — `confounded_by_co_prescription`.** 39 marginal → 7 adjusted under log10, but 0 under
  rank and arcsine transforms, 0 after demographic adjustment, 0 after adjusting for all
  medications, and 0 among the 260 statin users in the metformin-free stratum (a well-powered test).
  The residual handful is transform-dependent noise, not a reproducible effect; statins are
  co-prescribed with metformin in 41% of users and with other cardiometabolic drugs much more.
- **insulin — `confounded_by_co_prescription`.** 8–11 marginal → 0 adjusted; 70% of insulin users
  are on metformin.

`strongest_independent_drug` = **metformin** (99 independently associated species vs 12 for PPI).

`statin_metformin_overlap` = 180 / 440 = **0.409**.

## Files

- `work/analyze.py` — main pipeline (filtering, marginal + mutually adjusted models, answer).
- `work/sensitivity.py` — transform sensitivity, all-medication adjustment, logistic reverse model,
  stratified checks, surviving-species listings.
- `work/results_detail.json` — numeric counts and co-prescription table.
- `work/answer.json` — final answer.
