# Task A: mgx_drug_coprescription_confounding

**Task prompt:** Within the type 2 diabetes participants of a European cardiometabolic cohort, decide for each of five medications whether its association with gut species composition belongs to that drug or is inherited from a drug it is co-prescribed with. Report a label per drug, the strongest independent drug, the single variable that accounts for the sulfonylurea signal, and the fraction of statin users who also take metformin.

**Expected answer:** metformin and PPI `independent`; sulfonylurea, statin and insulin `not_independent`; strongest independent drug metformin; sulfonylurea signal explained by metformin; statin and metformin overlap 0.409.

**Dataset:** MetaPhlAn3 species profiles and participant metadata for the MetaCardis cohort (1,831 stool metagenomes), from curatedMetagenomicData, ExperimentHub EH7764, release 2022-10-19. Cohort paper: doi:10.1038/s41591-021-01593-2 (PMID 34880489).

**Skill explanation:** Tests whether the agent separates a drug effect from the prescribing pattern around it. Testing each drug on its own gives four of the five labels wrong, because statins, sulfonylureas and insulin are all entangled with metformin.

**LLM behavior:** Passed by every agent tested, 9 runs across 3 versions of the task, including a version with every hint about co-prescription removed. Agents parse the free text medication field unprompted and adjust the drugs against each other without being told to. Submitted as a calibration item.

---

## Full specification

## Task prompt

You are given shotgun metagenomic profiles from a European cardiometabolic cohort: a MetaPhlAn3 relative abundance table (`data/abundance_metaphlan3.tsv`, all taxonomic ranks, samples in columns, percent) and a participant metadata table (`data/sample_metadata.tsv`). The `treatment` column lists all current medications for that participant, separated by semicolons.

Restrict your analysis to participants whose `disease` field contains `T2D`. Each of the five drugs below shows an association with species level composition when it is tested on its own. Determine which of those associations are attributable to the drug itself.

Drugs and their tokens in the `treatment` column: metformin (`metformin`), proton pump inhibitor (`PPI`), sulfonylurea (`su`), statin (`statin`), insulin (`insulin`).

Use this significance criterion throughout: analyse at species level, test every species present in at least 10 percent of the T2D participants, apply Benjamini Hochberg correction across that entire set of tested species, and call a species significant at q < 0.05.

Label a drug `independent` if a substantial number of species remain significantly associated with it once everything else in the data that could account for the association is taken into account, and `not_independent` if its apparent association essentially disappears.

Report your answer as JSON with exactly these fields:

1. `drug_labels`: one label per drug, keyed by `metformin`, `proton_pump_inhibitor`, `sulfonylurea`, `statin`, `insulin`. Each value is `independent` or `not_independent`.
2. `strongest_independent_drug`: the drug with the largest number of independently associated species. One of the five drug keys above.
3. `drug_explaining_sulfonylurea_signal`: the single variable that accounts for the apparent sulfonylurea association. One of `metformin`, `proton_pump_inhibitor`, `dpp4_inhibitor`, `aspirin`, `age`, `BMI`, `country`, `none`.
4. `statin_metformin_overlap`: among T2D participants taking a statin, the fraction who are also taking metformin (0 to 1, three decimals).

Example format (values are placeholders):

```json
{
  "drug_labels": {"metformin": "...", "proton_pump_inhibitor": "...", "sulfonylurea": "...", "statin": "...", "insulin": "..."},
  "strongest_independent_drug": "...",
  "drug_explaining_sulfonylurea_signal": "...",
  "statin_metformin_overlap": 0.000
}
```

## Dataset

MetaCardis cohort (France, Germany, Denmark; 1,831 stool metagenomes, Ion Proton shotgun sequencing), as processed and released by curatedMetagenomicData (Bioconductor, release 2022-10-19; Pasolli et al. 2017, Nature Methods, doi:10.1038/nmeth.4468). Original cohort publication: Fromentin et al. 2022, Nature Medicine, PMID 34880489. Raw reads are in ENA.

Files given to the agent (built by `scripts/t1_make_bundle.py`):

| File | Contents |
|---|---|
| `data/abundance_metaphlan3.tsv` | 1,100 clades (kingdom to strain) x 1,831 samples, MetaPhlAn3 relative abundance in percent, exactly as distributed by curatedMetagenomicData (ExperimentHub EH7764) |
| `data/sample_metadata.tsv` | 1,831 samples x 14 fields: disease, age, gender, BMI, country, smoker, antibiotics_current_use, treatment, hba1c, ldl, triglycerides, hscrp, bristol_score, number_reads |

Study name, PMID, curator and read accession columns were removed, and sample IDs were replaced with neutral shuffled IDs (S0001 to S1831), so the source cohort cannot be identified from the files. The analysis population is the 792 participants whose `disease` field contains `T2D`.

## Expected answer

```json
{
  "drug_labels": {"metformin": "independent", "proton_pump_inhibitor": "independent", "sulfonylurea": "confounded_by_co_prescription", "statin": "confounded_by_co_prescription", "insulin": "confounded_by_co_prescription"},
  "strongest_independent_drug": "metformin",
  "statin_metformin_overlap": 0.409
}
```

Grading is deterministic (`grader.py`, `answer_key.json`): all seven fields must be correct; the overlap fraction is accepted within plus or minus 0.01 (180 of 440 T2D statin users also take metformin = 0.409).

### Robustness

From `scripts/ta_robustness.py`. Counts are species significant at BH q < 0.05 among the 195 species present in at least 10 percent of the 792 T2D participants. Six routes were tested, including two that use no regression at all:

| Drug | Users | R1 unadjusted | R2 adjusted for other drugs | R3 plus age, BMI, sex | R4 plus country | R5 among metformin users | R6 among metformin non users | Label |
|---|--:|--:|--:|--:|--:|--:|--:|---|
| metformin | 410 | 109 | 88 | 75 | 72 | n/a | n/a | independent |
| PPI | 150 | 13 | 14 | 12 | 11 | 8 | 13 | independent |
| sulfonylurea | 134 | 48 | 2 | 2 | 0 | 0 | n/a | confounded |
| statin | 440 | 38 | 0 | 0 | 0 | 0 | 0 | confounded |
| insulin | 118 | 11 | 0 | 0 | 0 | 0 | 0 | confounded |

Stratified routes are reported only where both strata have at least 20 participants.

Two details surfaced by the test runs are recorded here because they matter for grading:

- **Statin is the one borderline label.** It gives 0 species under CLR and rank based models and in both metformin strata, but one run found 7 residual species under a plain log10 transform, falling to 0 again once demographics or the remaining medications were added. The decisive evidence is the stratified test: among the 260 statin users within the 382 metformin free participants there are 0 significant species, which is well powered and transform free. The `confounded_by_co_prescription` label holds under every route except a single unadjusted transform choice, and that run reached the same label after checking.
- **Age and BMI are missing for every Danish participant.** All 199 Danish T2D participants lack both fields, while Germany (213) and France (380) are essentially complete. An agent that adds age and BMI to its model therefore drops a whole recruitment site by listwise deletion, taking the analysis set from 792 to 588 without warning. This is why the graded labels were verified under routes that do not use those covariates at all, including two stratified tests. One test run caught the missingness and left age and BMI out of its model, which is the right call here.
- **The token `su` is a substring of `insulin` and `sglt2`.** Medications must be matched as exact tokens after splitting on semicolons. Substring matching silently inflates the sulfonylurea group with insulin users, which blurs exactly the two labels the task is testing. One run flagged this explicitly as a necessary parsing step.

Two further drugs were tested and left out of the graded set because their label depended on the method: DPP4 inhibitors (28 species adjusted for drugs, but 2 after adding country) and aspirin (0 in the pooled adjusted models but 13 among metformin non users).

## Skill tested

Separating a drug's own effect from the prescribing pattern it sits in. The agent must parse a free text medication field into per drug indicators, recognize that the drugs are heavily co-prescribed, and adjust them against each other rather than testing each drug in isolation. Testing one drug at a time gives four of five labels wrong.

The statin result is the sharpest item. Statins are the most prescribed drug in the group (440 of 792), show 38 associated species on their own, and fall to zero under every adjusted or stratified route. The confounding also runs in a non obvious direction: statin users are *less* likely to be on metformin than non users (40.9 against 65.3 percent), so the marginal statin signal is partly an inverse shadow of the metformin signature rather than a shared one. An agent that reasons qualitatively that "statins and metformin go together, so adjust" reaches the right label by a route that does not actually describe this cohort.

## Expected LLM behavior

### Observed (6 runs over two prompt versions, 2026-10-22)

Each run was a fresh Claude Code agent with no conversation history, given only the prompt and the two data files, with Python available. Answers and notes are in `llm_runs/`.

Round 1 used an earlier prompt that stated outright that the drugs are frequently co-prescribed and that an association for one drug may reflect another, and used the label `confounded_by_co_prescription`. Round 2 is the prompt above, which removes that statement, removes the instruction to account for "the other medications", renames the label to the neutral `not_independent`, and adds the mechanism field.

| Round | Run | Model | Time | Result |
|---|---|---|--:|---|
| 1 with hint | a1_opus | Opus | ~2 min | PASS 7/7 |
| 1 with hint | a2_opus | Opus | ~3 min | PASS 7/7 |
| 1 with hint | a3_sonnet | Sonnet | ~35 s | PASS 7/7 |
| 2 no hint | b1_opus | Opus | ~3 min | PASS 8/8 |
| 2 no hint | b2_opus | Opus | ~3 min | PASS 8/8 |
| 2 no hint | b3_sonnet | Sonnet | ~70 s | PASS 8/8 |

**Removing the hint changed nothing.** Every round 2 run independently decided to adjust the drugs against each other without being told that co-prescription was the issue, and all three identified metformin as the variable behind the sulfonylurea signal by running a single covariate scan, which is the intended route. The hint was not load bearing.

All six runs built per drug indicators from the `treatment` strings by exact token matching, fitted one model per species with the drug indicators entered jointly, and recovered every label plus the overlap fraction. The stronger runs went well past what was asked: four abundance transforms, adjustment for every medication in the cohort rather than only the five named, a reverse direction logistic model, read depth as a covariate, complete case reanalysis against median imputation with missingness indicators, and stratified reanalysis inside and outside the metformin group.

Two of the data traps documented above were found by the agents rather than by the author. One run noticed that `su` substring matches `insulin` and switched to exact token matching. Another noticed that age and BMI are missing for the entire Danish subset and imputed with missingness indicators specifically so that country would not be dropped, while a third chose to leave those covariates out altogether for the same reason. Both are correct handling.

Several runs observed unprompted that statin use is *negatively* associated with metformin (40.9 against the cohort baseline of 51.8 percent), and one stated explicitly that the statin signal has no single driver, since metformin, aspirin, country, age and BMI each independently remove it. That is the correct reading, and it is why the mechanism field asks about sulfonylurea, where the answer is unique, rather than about statins. One run also noted that 66 percent of PPI users take a statin and that PPI retains the canonical oral to gut taxa (*Streptococcus salivarius*, *Veillonella*, *Rothia mucilaginosa*, *Haemophilus parainfluenzae*), which is the correct biological reading of the one non obvious `independent` label.

### Interpretation

Current frontier agents solve this reliably. The same held for an earlier T2D versus healthy design in this family, which asked which species differences survive once medication is accounted for. Nine graded runs across the three versions, no failures.

Drug confounding is a capability these models have rather than a gap. They parse free text medication fields unprompted, adjust mutually rather than stratifying naively, notice structural missingness, and do not substitute the literature prior for the data: several runs on the earlier designs correctly declined to reproduce the mouse finding that metformin raises *Akkermansia*, which this cohort does not support.

The task is therefore submitted as a calibration item rather than a discriminating one. It remains useful in a benchmark: it fails any agent that tests drugs one at a time, which is the common shortcut and gives four of five labels wrong, it carries a deterministic numeric anchor that cannot be recalled from any publication, and it contains three independent parsing and missingness traps. The companion task (`mgx_crc_partial_contamination`) is the discriminating item, and the contrast between the two is the main empirical finding of this submission.
