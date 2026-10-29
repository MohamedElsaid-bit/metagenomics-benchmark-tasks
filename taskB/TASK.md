# Task B: mgx_crc_partial_contamination

**Task prompt:** Given a colorectal cancer case control study of stool metagenomes, determine which species differ between cases and controls. Report any samples that should be excluded as unsuitable, the origin of the problem material in them, and the direction of difference for six named species.

**Expected answer:** 21 samples excluded, oral in origin; *F. nucleatum*, *P. stomatis* and *P. micra* enriched in CRC; *R. dentocariosa*, *E. rectale* and *S. salivarius* not significant.

**Dataset:** A CRC stool cohort (ZellerG_2014, ExperimentHub EH5938) in which 21 profiles were mixed with real subgingival plaque profiles (GhensiP_2019, ExperimentHub EH5542) at 26 to 40 percent, 20 of the 21 in the case arm. Both from curatedMetagenomicData, release 2021-03-31. Source papers: doi:10.15252/msb.20145645 and doi:10.1038/s41522-020-00155-7.

**Skill explanation:** Tests whether the agent checks sample provenance before interpreting a case control result, in a setting where the contamination reproduces what the literature predicts. The oral contamination inflates the canonical CRC markers and simultaneously dilutes gut commensals, manufacturing a false butyrate producer depletion. If the contamination is left in, three of the six graded species flip under a Wilcoxon test and two under a CLR Welch t test.

**LLM behavior:** Opus passed 2 of 2, Sonnet failed. All three got the statistics right; the split was on detection. Sonnet scored whole genera, which pulls in gut residents such as *Prevotella copri*, found a fuzzy boundary and left two contaminated samples in. PCA does not separate the contaminated samples, so ordination shortcuts fail.

---

## Full specification

## Task prompt

You are given a shotgun metagenomic case control study of colorectal cancer (CRC): a species level relative abundance table (`data/abundance_species.tsv`, species in rows, samples in columns, percent) and a sample table (`data/sample_metadata.tsv`) with group (`CRC` or `control`), age and gender. All samples are recorded as stool.

Determine which species differ between CRC cases and controls.

Use this significance criterion: test every species present in at least 10 percent of the samples you include, apply Benjamini Hochberg correction across that entire set of tested species, and call a species significant at q < 0.05.

Report your answer as JSON with exactly these fields:

1. `n_flagged_samples`: the number of samples you would exclude from the case control comparison as unsuitable, after checking the data. Use `0` if you would exclude none.
2. `flagged_sample_ids`: the list of those sample IDs (empty list if none).
3. `flagged_sample_origin`: if you excluded samples, the most likely biological origin of the problem material in them. One of `oral`, `skin`, `vaginal`, `nasal`, `soil`, `reagent_blank`, `not_applicable`.
4. `crc_vs_control`: for each of the six species below, the difference between CRC cases and controls in the comparison you consider correct. One of `enriched_in_CRC`, `depleted_in_CRC`, `not_significant`.

Species: `Fusobacterium_nucleatum`, `Peptostreptococcus_stomatis`, `Parvimonas_micra`, `Rothia_dentocariosa`, `Eubacterium_rectale`, `Streptococcus_salivarius`.

Example format (values are placeholders):

```json
{
  "n_flagged_samples": 0,
  "flagged_sample_ids": [],
  "flagged_sample_origin": "not_applicable",
  "crc_vs_control": {"Fusobacterium_nucleatum": "...", "Peptostreptococcus_stomatis": "...", "Parvimonas_micra": "...", "Rothia_dentocariosa": "...", "Eubacterium_rectale": "...", "Streptococcus_salivarius": "..."}
}
```

## Dataset

Built by `scripts/t3_build_and_key.py` from two curatedMetagenomicData studies (Bioconductor, release 2021-03-31; Pasolli et al. 2017, Nature Methods, doi:10.1038/nmeth.4468):

| Source | Role |
|---|---|
| ZellerG_2014 (French CRC cohort, stool; Zeller et al. 2014, Mol Syst Biol, doi:10.15252/msb.20145645), ExperimentHub EH5938 | all 114 samples (53 CRC, 61 control) |
| GhensiP_2019 (subgingival plaque; Ghensi et al. 2020, NPJ Biofilms Microbiomes, doi:10.1038/s41522-020-00155-7), ExperimentHub EH5542 | oral material mixed into 21 of those stool profiles |

Each affected sample is a genuine stool profile mixed with a genuine plaque profile at a fraction drawn uniformly from 26 to 40 percent (`new = (1 - f) * stool + f * plaque`, renormalized to 100 percent). This models partial carry over, for example saliva contamination during collection or index hopping, rather than the whole sample swap used in the companion task `mgx_crc_sample_provenance`. 20 of the 21 affected samples are CRC cases and 1 is a control, which is what makes the contamination produce a false disease signal.

Plaque donors were restricted to profiles whose oral marker content is at least 40 percent, so the mixed in material is reliably oral. The contaminated CRC samples were drawn from cases carrying neither *Fusobacterium nucleatum* nor *Peptostreptococcus stomatis*, so that removing them leaves the genuine CRC signal intact. Sample IDs are neutral and shuffled (P001 to P114). Truth is kept in `truth_PRIVATE.json` and is not given to the agent.

## Expected answer

```json
{
  "n_flagged_samples": 21,
  "flagged_sample_ids": ["P005","P009","P013","P016","P022","P029","P030","P037","P040","P043","P046","P047","P049","P056","P062","P094","P099","P107","P110","P111","P114"],
  "flagged_sample_origin": "oral",
  "crc_vs_control": {"Fusobacterium_nucleatum": "enriched_in_CRC", "Peptostreptococcus_stomatis": "enriched_in_CRC", "Parvimonas_micra": "enriched_in_CRC", "Rothia_dentocariosa": "not_significant", "Eubacterium_rectale": "not_significant", "Streptococcus_salivarius": "not_significant"}
}
```

Grading is deterministic (`grader.py`, `answer_key.json`). All ten checks must pass:

- every one of the 21 contaminated samples must be flagged (recall 1.0)
- at most 18 additional samples may be flagged, a guard set well above any realistic answer
- `n_flagged_samples` must equal the length of the submitted ID list
- `flagged_sample_origin` must be `oral`
- all six species calls must be correct

**Why the flag set is not graded as an exact match.** The two errors are not symmetric. Leaving a contaminated sample in biases the case control comparison, because the contamination is concentrated in the case arm. Dropping a clean sample only costs statistical power. There is also a defensible range of marker definitions: P031 carries *Prevotella intermedia* at 11.3 percent and P044 carries *Treponema berlinense* at 7.3 percent, both periodontal organisms, in genuinely uncontaminated stool, so an agent using a broader oral panel will flag them on reasonable grounds.

**Why the guard on extras is set where it is.** An earlier version capped extras at 3. That was wrong, because the cap itself decided an outcome: one test run flagged exactly 3 extra samples and another flagged 2, so moving the threshold by one would have flipped a pass to a fail. A grading threshold that determines a result is a threshold chosen badly.

The guard is now 18, which is 20 percent of the 93 clean samples. It exists only to reject degenerate answers, such as flagging half the cohort to guarantee recall, and sits far above anything a real analysis produces: the three test runs flagged 0, 2 and 3 extras. Over flagging within the guard is judged where it should be, in the species calls, since an exclusion set that distorts the comparison will produce wrong calls and fail on those fields. This follows the endpoint grading philosophy of the scBench-Long benchmark, which does not require intermediate analyses to follow a prescribed workflow and lets an agent reach the target surface by any valid path. `scripts/t3_flagset_sensitivity.py` confirms the species calls are unchanged across 16 plausible flag sets, including sets with 5 false positives.

### Why this is harder than whole sample contamination

In the companion task the contaminated samples are pure plaque, so they contain no *Bacteroides* and no *Faecalibacterium* at all and fall out of any ordination immediately. Here the affected samples are still 60 to 74 percent genuine stool, so those shortcuts fail:

| Detection route | Companion task (whole sample) | This task (partial) |
|---|---|---|
| Gut core taxa absent | 0 percent in all 12, unmissable | 13 to 58 percent, fully overlapping genuine stool |
| PCA on CLR, first 4 components | separates cleanly | **overlaps on every component** |
| Oral specific species score | 0.0 to 0.3 percent in stool vs 58 to 99 percent | 0.0 to 3.0 percent in stool vs 12.8 to 33.1 percent |

Detection requires going looking for oral specific taxa. The remaining gap is wide: `scripts/t3_flagset_sensitivity.py` shows that **any threshold from 3 to 12 percent on a summed oral species score returns exactly the 21 affected samples**, with no false positives and no misses. The answer is therefore well defined, but only for an agent that scores oral *species*. Scoring whole *genera* blurs the boundary, because *Prevotella copri*, *Streptococcus thermophilus* and *Veillonella* are ordinary gut residents. With the genus panel in `scripts/t3_flagset_sensitivity.py`, a cut at 10 percent flags all 21 plus 20 clean samples and a cut at 15 percent flags all 21 plus 12; there is no threshold that separates cleanly. (An earlier draft of this file quoted a genus result that did not reproduce from that script; corrected 2026-10-29.)

### Robustness of the species calls

All six calls were computed under Wilcoxon and CLR Welch t test at prevalence filters of 5, 10 and 20 percent, and across 16 plausible flag sets (threshold rules, and sets perturbed to miss or over call samples):

| Species | Naive (nothing removed) | Correct (21 removed) | Stability |
|---|---|---|---|
| *Fusobacterium nucleatum* | enriched | enriched | all methods |
| *Peptostreptococcus stomatis* | enriched | enriched | all methods |
| *Parvimonas micra* | enriched | enriched | all methods |
| *Rothia dentocariosa* | enriched (q = 0.02) | not significant (q = 0.50) | all methods |
| *Eubacterium rectale* | depleted under Wilcoxon (q = 0.01); not significant under CLR t | not significant (q = 0.06 to 0.11) | cleaned call stable; naive call depends on the test |
| *Streptococcus salivarius* | depleted (q = 0.01) | not significant (q = 0.14) | all methods |

The three true markers survive cleaning by construction, and the artefacts disappear. The cleaned calls are identical under both tests, so the answer key does not depend on the method. Sensitivity to incomplete cleaning does depend on it: under Wilcoxon, flag sets that recover 18 or more of the 21 affected samples reproduce all six calls, while under CLR t missing even 3 flips *S. salivarius* back to depleted. Sets that recover only 14 flip *E. rectale* or *S. salivarius* under both. This is the expected behaviour of an under cleaned analysis rather than a grading defect, and it is why the grader requires full recall.

## Skill tested

Detecting partial cross contamination that is correlated with case status, and distinguishing genuine disease biology from the artefact it manufactures. The agent must notice that nothing in the data obviously fails QC, construct a targeted oral taxa score using species that genuinely do not reside in the gut, find the threshold from the distribution rather than by guesswork, and then recompute the comparison.

The design works against the agent's prior in two directions at once. The contamination inflates oral taxa in the case group, which mimics the well known oral to gut translocation signature of CRC. It simultaneously dilutes ordinary gut commensals in those same samples, which manufactures an apparent depletion of butyrate producers (*E. rectale*, *S. salivarius*), another textbook CRC claim. An agent that treats every oral looking signal as contamination overcorrects and loses the three genuine markers; an agent that treats the result as biology reports two classic CRC findings that are pure artefact.

## Expected LLM behavior

### Observed (3 runs, 2026-10-22)

Each run was a fresh Claude Code agent with no conversation history, given only the prompt and the two data files, with Python available. Answers and the agents' own notes are in `llm_runs/`.

| Run | Model | Contaminated samples found | Extra flagged | Species calls | Result |
|---|---|---|--:|---|---|
| r1_opus | Opus | 21 of 21 | 3 | 6 of 6 | PASS 10/10 |
| r2_opus | Opus | 21 of 21 | 2 | 6 of 6 | PASS 10/10 |
| r3_sonnet | Sonnet | 19 of 21 | 0 | 6 of 6 | FAIL 9/10 |

**All three got every species call right, and the task still separated them.** The discriminating step was detection, not statistics. Both Opus runs built a panel of obligate oral taxa, found the bimodal split, and recovered every contaminated sample. One of them explicitly excluded *F. nucleatum*, *P. micra* and *P. stomatis* from its own detection panel on the grounds that those genuinely colonize the CRC gut and using them would be circular, which is the correct reasoning and is what keeps the three true markers intact.

The Sonnet run scored whole genera rather than oral specific species. That sweeps in *Prevotella copri*, *Streptococcus thermophilus* and *Veillonella*, which are ordinary gut residents, so its score was noisy at the boundary. It described the 7 to 9 percent range as "ambiguous", kept samples there, and missed P016 and P029 (true oral scores 13.8 and 18.7 percent). It still reached the right species calls, because 19 of 21 removals is enough for the artefacts to disappear, but it reported a comparison that still contained two contaminated samples.

Both Opus runs independently noted that the contamination is concentrated in the case arm (they counted 22 or 23 of their flagged samples as CRC) and stated that this is what manufactures the false signal. Both also observed that the three genuine markers become *more* significant after cleaning, and used that as evidence that the exclusion removes noise rather than signal.

**Difficulty compared with the whole sample version.** The companion task is passed by both models once the significance criterion is specified. This one is passed by Opus and failed by Sonnet, with the failure falling on the QC step rather than the statistics. Run times roughly doubled (6 to 8 minutes against 1 to 4). The contamination here cannot be found by ordination, so the task rewards knowing which organisms are oral specific, which is domain knowledge rather than method.

### A defect found during testing

The first version of `answer_key.json` listed the contaminated samples from an earlier build of the dataset. The build script had been re run with different parameters, which regenerated both the contamination assignment and the anonymized IDs, and the key was not regenerated with it. The first grading pass therefore scored all three runs against the wrong ID list and reported three failures that were not real. It was caught because a model flagged a sample that the stale key called clean while the data showed *Rothia dentocariosa* at 14.7 percent in it, which could not both be true.

The key is now generated from `truth_PRIVATE.json`, which the build script writes in the same run as the data, so the two cannot diverge. The general lesson for task authoring is that an answer key derived by hand from a console log is a latent defect; it has to be emitted by the same code path that produces the data.
