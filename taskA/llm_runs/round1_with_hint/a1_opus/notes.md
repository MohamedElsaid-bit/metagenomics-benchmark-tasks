# Drug–microbiome associations in the T2D subgroup: univariate vs. co-medication-adjusted

Script: `work/analysis.py` (run with the WSL `ml-deepen` venv). Outputs: `work/answer.json`,
`work/summary_counts.csv`, `work/coprescription_overlap.csv`.

## Cohort and feature set
- Participants whose `disease` field contains `T2D`: **n = 792** (T2D, CAD;T2D, HF;CAD;T2D, T2D;HF, HF;T2D),
  all present in the abundance table.
- Species level = `clade_name` rows containing `|s__` and no `|t__` → 713 species.
- Prevalence filter: species with relative abundance > 0 in ≥ 10 % of the 792 T2D participants →
  **195 species tested**. BH FDR is applied across all 195 tested species, separately for each
  drug/model; significance = q < 0.05.
- Drug exposure = exact token match after splitting `treatment` on `;`.
  Counts in the T2D subgroup: metformin 410, statin 440, PPI 150, sulfonylurea 134, insulin 118.

## Statistical model
Abundances are zero-inflated and compositional, so each species' abundance vector was
rank-based inverse-normal transformed across the 792 participants. Then:

- **Univariate**: one drug indicator per model (OLS on the transformed abundance). Cross-checked with
  Mann–Whitney U on the raw abundances — the two give nearly identical significant-species counts.
- **Multivariate (co-medication-adjusted)**: a single OLS per species with **all five drug
  indicators entered jointly**, so each drug's coefficient is its partial association holding the
  other four medications constant.

## Co-prescription structure (P(column | row), T2D subgroup)
|            | metformin | PPI  | su   | statin | insulin |
|------------|-----------|------|------|--------|---------|
| metformin  | 1.000 | 0.144 | 0.283 | 0.439 | 0.202 |
| PPI        | 0.393 | 1.000 | 0.147 | 0.660 | 0.133 |
| su         | **0.866** | 0.164 | 1.000 | 0.522 | 0.231 |
| statin     | 0.409 | 0.225 | 0.159 | 1.000 | 0.148 |
| insulin    | **0.703** | 0.169 | 0.263 | 0.551 | 1.000 |

Sulfonylurea (86.6 %) and insulin (70.3 %) users are overwhelmingly also on metformin — the
classic setup for a borrowed association.

## Significant species (q < 0.05, BH over 195 species)
| drug | n taking | univariate (OLS) | univariate (MWU) | adjusted for other 4 drugs | retained |
|------|----------|------------------|------------------|----------------------------|----------|
| metformin | 410 | 110 | 109 | **95** | 86 % |
| PPI | 150 | 14 | 13 | **12** | 86 % |
| sulfonylurea | 134 | 49 | 48 | **1** | 2 % |
| statin | 440 | 32 | 38 | **0** | 0 % |
| insulin | 118 | 7 | 11 | **1** | 14 % |

### Sensitivity checks (same filter and FDR scheme)
- log10(abundance + 1e-5) instead of rank transform, multivariate: metformin 99, PPI 12,
  sulfonylurea 2, statin 6, insulin 0 — same qualitative split.
- Multivariate plus age, sex, BMI, country and current antibiotics: metformin 82, PPI 13,
  sulfonylurea 1, statin 0, insulin 0.

## Interpretation
- **metformin — independent.** 95/110 hits survive. Signature is the expected one: ↑ *Escherichia
  coli*, *Clostridium clostridioforme*, *C. bolteae*, *Flavonifractor plautii*; ↓ *Intestinibacter
  bartlettii*, *Faecalibacterium prausnitzii*, *Ruminococcus bicirculans*.
- **PPI — independent.** 12/14 hits survive, and they are the canonical oral-to-gut translocation
  taxa: *Streptococcus parasanguinis*, *S. salivarius*, *Veillonella parvula/atypica/dispar*,
  *Rothia mucilaginosa*, *Bifidobacterium dentium*, *Haemophilus parainfluenzae*. PPI use is only
  weakly tied to the antidiabetics, so its signal is structurally distinct.
- **sulfonylurea — confounded by co-prescription.** 49 → 1. Its univariate hits are metformin's
  hits, inherited through 86.6 % co-prescription.
- **statin — confounded by co-prescription.** 32–38 → 0. Statin is the most common drug here (440)
  and sits downstream of both metformin (41 % overlap) and PPI (66 % of PPI users are on a statin),
  so its marginal signal is fully absorbed once those are in the model.
- **insulin — confounded by co-prescription.** 7–11 → 1 (a single taxon, *Butyricimonas virosa*),
  which is not a "substantial number"; 70 % of insulin users are on metformin.

Labelling rule applied: `independent` if ≥ 5 species remain significant after adjustment
(metformin 95, PPI 12) — the gap to the others (0–1) is unambiguous, so the exact cutoff is not
load-bearing.

## Requested statistic
Among T2D statin users (n = 440), 180 also take metformin → **0.409** (verified independently with awk).
