# Task: mgx_t2d_drug_confounding

## Task prompt

You are given shotgun metagenomic profiles from a large European cardiometabolic cohort: a MetaPhlAn3 relative abundance table (all taxonomic ranks, samples in columns) and a per participant metadata table. Participants include healthy controls and patients with metabolic syndrome, impaired glucose tolerance, type 2 diabetes (T2D), coronary artery disease and heart failure. The `treatment` column lists current medications separated by semicolons.

A frequently reported gut microbiome signature of T2D is enrichment of *Escherichia coli* together with depletion of butyrate producers such as *Faecalibacterium prausnitzii* and *Roseburia hominis*. Characterize the species level differences between T2D participants and healthy controls in this cohort and determine what drives them.

Define T2D participants as everyone whose `disease` field contains `T2D`, and healthy controls as everyone whose `disease` field is exactly `healthy`.

Report your answer as JSON with exactly these fields:

1. `primary_driver`: the single factor that best explains the species level differences between T2D participants and healthy controls. One of: `type_2_diabetes`, `metformin`, `statin`, `dpp4_inhibitor`, `sulfonylurea`, `insulin`, `proton_pump_inhibitor`, `aspirin`, `BMI`, `age`, `country`.
2. `driver_fraction_T2D`: the fraction of T2D participants exposed to your `primary_driver` (0 to 1, three decimals), or `null` if it is not a yes/no exposure (`type_2_diabetes`, `BMI`, `age`, `country`).
3. `T2D_intrinsic`: for each of the five species below, the difference attributable to T2D itself, comparing T2D participants with healthy controls after accounting for any confounders you identify. One of `increased`, `decreased`, `no_robust_difference`.
4. `driver_effect_within_T2D`: for each of the five species, the difference among T2D participants exposed to your `primary_driver` compared with T2D participants not exposed. One of `increased`, `decreased`, `no_robust_difference`, or `null` for every species if `primary_driver` is not a yes/no exposure.

Species: `Escherichia_coli`, `Intestinibacter_bartlettii`, `Faecalibacterium_prausnitzii`, `Roseburia_hominis`, `Akkermansia_muciniphila`.

Example format (values are placeholders):

```json
{
  "primary_driver": "...",
  "driver_fraction_T2D": 0.000,
  "T2D_intrinsic": {"Escherichia_coli": "...", "Intestinibacter_bartlettii": "...", "Faecalibacterium_prausnitzii": "...", "Roseburia_hominis": "...", "Akkermansia_muciniphila": "..."},
  "driver_effect_within_T2D": {"Escherichia_coli": "...", "Intestinibacter_bartlettii": "...", "Faecalibacterium_prausnitzii": "...", "Roseburia_hominis": "...", "Akkermansia_muciniphila": "..."}
}
```

## Dataset

MetaCardis cohort (France, Germany, Denmark; 1,831 stool metagenomes, Ion Proton shotgun sequencing), as processed and released by curatedMetagenomicData (Bioconductor, release 2022-10-19; Pasolli et al. 2017, Nature Methods, doi:10.1038/nmeth.4468). Original cohort publication: Fromentin et al. 2022, Nature Medicine, PMID 34880489. Raw reads: ENA (accessions in curatedMetagenomicData metadata).

Files given to the agent (built by `scripts/t1_make_bundle.py`):

| File | Contents |
|---|---|
| `data/abundance_metaphlan3.tsv` | 1,100 clades (kingdom to strain) x 1,831 samples, MetaPhlAn3 relative abundance in percent, exactly as distributed by curatedMetagenomicData (ExperimentHub EH7764) |
| `data/sample_metadata.tsv` | 1,831 samples x 14 fields: disease, age, gender, BMI, country, smoker, antibiotics_current_use, treatment, hba1c, ldl, triglycerides, hscrp, bristol_score, number_reads |

Study name, PMID, curator and read accession columns were removed, and sample IDs were replaced with neutral IDs (S0001 to S1831, shuffled), so the source cohort cannot be identified from the files.

## Expected answer

```json
{
  "primary_driver": "metformin",
  "driver_fraction_T2D": 0.518,
  "T2D_intrinsic": {"Escherichia_coli": "increased OR no_robust_difference", "Intestinibacter_bartlettii": "no_robust_difference", "Faecalibacterium_prausnitzii": "no_robust_difference", "Roseburia_hominis": "no_robust_difference", "Akkermansia_muciniphila": "no_robust_difference"},
  "driver_effect_within_T2D": {"Escherichia_coli": "increased", "Intestinibacter_bartlettii": "decreased", "Faecalibacterium_prausnitzii": "decreased", "Roseburia_hominis": "decreased", "Akkermansia_muciniphila": "no_robust_difference"}
}
```

Grading is deterministic (`grader.py`, `answer_key.json`): every field must match; the fraction is accepted within plus or minus 0.01 (410 of 792 T2D participants take metformin = 0.518). Pass requires all 12 fields.

### How the answer was derived and checked for robustness

All numbers below come from `scripts/t1_explore.py`, `t1_robustness.py`, `t1_drugs.py` and `t1_country.py`. Species abundances were renormalized to the species total; tests used species present in at least 10% of the compared samples; multiple testing used Benjamini Hochberg at q < 0.05.

| Comparison | Significant species |
|---|--:|
| All T2D (792) vs healthy (234), Wilcoxon | 89 |
| Metformin T2D (410) vs healthy | 108 |
| Untreated T2D (317) vs healthy | 30 |
| Metformin vs no metformin within T2D | 109 |
| Same, CLR linear model adjusted for country, age, BMI | 93 |
| Untreated T2D vs healthy, adjusted for country, age, BMI | 1 |

Other drugs within T2D, adjusted for metformin, age and BMI: DPP4 inhibitor 26, PPI 12, sulfonylurea 1, statin 0, aspirin 0, insulin 0 (unadjusted, statins showed 38 and aspirin 66, almost entirely through co prescription with metformin). Metformin is the dominant factor under every method.

Each graded species call was required to hold under five methods (Wilcoxon at prevalence filters of 5%, 10% and 25%; CLR Welch t test; CLR linear model adjusted for age, BMI and sex) plus country adjustment. *Akkermansia muciniphila* shows no metformin effect under any method (q from 0.43 to 0.95).

`T2D_intrinsic` was checked across the different routes an agent could reasonably take to isolate T2D itself (`scripts/t1_alt_routes.py`; CLR linear models, all adjusted for country, age and BMI):

| Route | n | Significant species | Five graded species |
|---|--:|--:|---|
| Untreated T2D (no glucose lowering drug) vs healthy | 297 | 1 | all q > 0.6 |
| All T2D vs healthy, metformin as covariate | 741 | 7 | all q > 0.27 |
| All T2D vs healthy, every glucose lowering drug class as covariate | 741 | 7 | all q > 0.32 |
| T2D not on metformin vs healthy | 362 | 7 | all q > 0.53 |

*E. coli* reaches q = 0.03 in one unadjusted comparison (untreated T2D vs healthy, Wilcoxon), so both `increased` and `no_robust_difference` are accepted for it. Dropping age and BMI from the model makes *E. coli* (q = 0.005, up) and *Akkermansia* (q = 0.027, down) look T2D associated. That is not accepted: healthy controls average BMI 22.8 against 32.2 in T2D, and lower *Akkermansia* is a known correlate of adiposity, so leaving BMI out is a modeling error the task is designed to catch.

Known caveat: within Germany alone, *E. coli* trends slightly lower in metformin users (q = 0.016, n = 161 vs 52), while France, Denmark, the pooled analysis and the country adjusted model all show a strong increase. *Intestinibacter bartlettii* depletion with metformin is the most robust call, significant within each of the three countries.

## Skill tested

Recognizing that a disease signature in an observational cohort can be a treatment signature. A capable agent should parse the medication field, discover without being told that roughly half of T2D participants take metformin, stratify or adjust for it (and for age, BMI and recruitment country, which also differ between groups), see through statin and aspirin associations that exist only through co prescription with metformin, and report which species differences remain attributable to T2D itself. It should also resist two priors: that butyrate producer depletion is intrinsic to T2D, and that metformin raises *Akkermansia* (a finding from mouse studies that this human cohort does not reproduce).

## Expected LLM behavior

### Observed (3 runs, 2026-10-22)

Each run was a fresh Claude Code agent with no conversation history, given only the prompt and the two data files in an empty folder, with Python (pandas, scipy, statsmodels) available. Answers and the agents' own notes are in `llm_runs/`.

| Run | Model | Time | Result |
|---|---|--:|---|
| run1_opus | Claude Opus | ~1.7 min | PASS, 12/12 |
| run2_opus | Claude Opus | ~1.5 min | PASS, 12/12 |
| run3_sonnet | Claude Sonnet | ~1.0 min | PASS, 12/12 |

All three runs followed the same path: they parsed the `treatment` column without being told to, noticed that 410 of 792 T2D participants take metformin while almost no controls do, fit per species models on log abundance adjusted for age, BMI, sex and country, tested every drug while adjusting for the others (correctly finding statins and aspirin null once metformin is in the model), and compared T2D participants not on metformin with controls as a second check. All three correctly reported no metformin effect on *Akkermansia*. The Sonnet run hit rank deficiency warnings in its adjusted model and relied on the metformin free comparison instead; it reported *E. coli* as `increased` for `T2D_intrinsic`, which the key accepts.

One Opus run named the cohort and paper ("MetaCardis, Forslund et al.") in its notes, most likely recognizing the original sample ID prefix. Sample IDs have since been anonymized in the released bundle; the three runs above used the original IDs.

### Interpretation

Current frontier agents solve this task reliably. The reason is instructive: metformin confounding of T2D microbiome signatures is a well known result (Forslund et al. 2015), so the model's prior points toward the correct analysis rather than away from it. The task still has value as a calibration item that checks basic confounder handling (drug co prescription, BMI, country), and an agent that stops at the first T2D vs healthy comparison would fail it (conjecture; no weaker model was tested), since that comparison alone yields 89 species and the published signature. But it does not separate strong agents. A harder version would need a confounder the model has not memorized, or a setting where the literature prior is wrong; the companion task in this submission is designed around that principle.

---

**Note on this folder.** This is the first version of the confounding task; the submitted version is `taskA`. The data bundle has been removed from here because it is byte identical to `taskA/data` (same MetaCardis export, same anonymized sample IDs, built by the same script). The prompt, answer key, grader and the three test runs are kept as the record of how the task evolved.
