# Task: mgx_crc_sample_provenance

## Task prompt

You are given a shotgun metagenomic case control study of colorectal cancer (CRC): a species level relative abundance table (`data/abundance_species.tsv`, species in rows, samples in columns, percent) and a sample table (`data/sample_metadata.tsv`) with group (`CRC` or `control`), age and gender. All samples are recorded as stool.

Determine which species differ between CRC cases and controls.

Use this significance criterion: test every species present in at least 10 percent of the samples you include, apply Benjamini Hochberg correction across that entire set of tested species, and call a species significant at q < 0.05.

Report your answer as JSON with exactly these fields:

1. `n_flagged_samples`: the number of samples you would exclude from the case control comparison as unsuitable, after checking the data. Use `0` if you would exclude none.
2. `flagged_sample_ids`: the list of those sample IDs (empty list if none).
3. `flagged_sample_origin`: if you excluded samples, the most likely biological origin of the material in them. One of `oral`, `skin`, `vaginal`, `nasal`, `soil`, `reagent_blank`, `not_applicable`.
4. `crc_vs_control`: for each of the five species below, the difference between CRC cases and controls in the comparison you consider correct. One of `enriched_in_CRC`, `depleted_in_CRC`, `not_significant`.

Species: `Fusobacterium_nucleatum`, `Peptostreptococcus_stomatis`, `Parvimonas_micra`, `Gemella_morbillorum`, `Faecalibacterium_prausnitzii`.

Example format (values are placeholders):

```json
{
  "n_flagged_samples": 0,
  "flagged_sample_ids": [],
  "flagged_sample_origin": "not_applicable",
  "crc_vs_control": {"Fusobacterium_nucleatum": "...", "Peptostreptococcus_stomatis": "...", "Parvimonas_micra": "...", "Gemella_morbillorum": "...", "Faecalibacterium_prausnitzii": "..."}
}
```

## Dataset

Built by `scripts/t2_build_and_key.py` from two curatedMetagenomicData studies (Bioconductor, release 2021-03-31; Pasolli et al. 2017, Nature Methods, doi:10.1038/nmeth.4468), both MetaPhlAn3 profiles of real human shotgun metagenomes:

| Source | Role | n |
|---|---|--:|
| ZellerG_2014 (French CRC cohort, stool; Zeller et al. 2014, Mol Syst Biol, doi:10.15252/msb.20145645) | 53 CRC cases, 61 controls, unchanged | 114 |
| GhensiP_2019 (subgingival plaque, peri implant study; Ghensi et al. 2020, NPJ Biofilms Microbiomes, doi:10.1038/s41522-020-00155-7) | 12 plaque profiles relabelled as CRC stool cases | 12 |

Total 126 samples, 798 species. Sample IDs were replaced with neutral IDs (C001 to C126, shuffled); the 12 plaque samples were given ages drawn from the real CRC age distribution and random sex, so metadata alone does not identify them. Abundances were renormalized to 100 percent per sample after merging. The adenoma samples of the Zeller cohort were dropped so the comparison is a clean two group design. The plaque profiles are unmodified real data, not simulated: this mimics a sample swap or mislabelling event rather than low level cross contamination.

Files given to the agent: `data/abundance_species.tsv`, `data/sample_metadata.tsv`. The mapping back to source IDs is kept in `id_map_PRIVATE.tsv` and is not given to the agent.

## Expected answer

```json
{
  "n_flagged_samples": 12,
  "flagged_sample_ids": ["C016", "C020", "C023", "C036", "C041", "C051", "C057", "C060", "C075", "C083", "C092", "C100"],
  "flagged_sample_origin": "oral",
  "crc_vs_control": {"Fusobacterium_nucleatum": "enriched_in_CRC", "Peptostreptococcus_stomatis": "enriched_in_CRC", "Parvimonas_micra": "not_significant", "Gemella_morbillorum": "not_significant", "Faecalibacterium_prausnitzii": "not_significant"}
}
```

Grading is deterministic (`grader.py`, `answer_key.json`): the flagged ID set must match exactly, and all eight fields must be correct.

### Why these species, and robustness

From `scripts/t2_robustness.py`. Every species call was computed six ways (Wilcoxon and CLR Welch t test, each at prevalence filters of 5, 10 and 20 percent; Benjamini Hochberg q < 0.05), with and without the 12 plaque samples.

| Species | With plaque samples | After removing them | Graded as |
|---|---|---|---|
| *Fusobacterium nucleatum* | enriched, 6/6 | enriched, 4/4 tested | `enriched_in_CRC` |
| *Peptostreptococcus stomatis* | enriched, 6/6 | enriched, 6/6 | `enriched_in_CRC` |
| *Parvimonas micra* | enriched, 6/6 | not significant, 6/6 | `not_significant` |
| *Gemella morbillorum* | enriched, 6/6 | not significant, 6/6 | `not_significant` |
| *Faecalibacterium prausnitzii* | depleted, 6/6 | not significant, 6/6 | `not_significant` |

*F. nucleatum* falls below the 20 percent prevalence filter in the cleaned data and is not tested there; it is significant in all four runs that do test it. Species whose calls were method dependent (*Porphyromonas asaccharolytica*, *Solobacterium moorei*, *Eubacterium rectale*, *Bacteroides vulgatus*) were deliberately left out of the graded set.

The prompt fixes the significance criterion (prevalence filter 10 percent, Benjamini Hochberg across all tested species, q < 0.05) because the first round of test runs showed this is genuinely ambiguous otherwise: correcting across only the five named species makes *Parvimonas micra* (q = 0.010) and *Gemella morbillorum* (q = 0.018) significant, while correcting across all 267 tested species leaves them at q = 0.085 and 0.141. Both are defensible analyses, so the task must say which one it wants.

The total number of significant species is not graded, because it is strongly method dependent: 5 to 12 in the cleaned comparison against 32 to 81 with the plaque samples included. The direction of that change (far fewer findings after cleaning) is stable, but the count is not.

### The flagged set is objectively well defined

Exact set grading is only fair if the 12 samples are unambiguously separable, which `scripts/t2_separability.py` confirms:

| Measure | 12 planted samples | 114 genuine stool samples |
|---|---|---|
| Obligate gut anaerobe content (Bacteroides, Faecalibacterium, Roseburia, Eubacterium, Alistipes, Blautia, Ruminococcus and similar) | 0.00 to 11.08 percent | 26.09 to 97.44 percent |
| *Bacteroides* total | 0.000 to 0.022 percent | 0.059 percent and above |
| *Faecalibacterium prausnitzii* | exactly 0 in all 12 | nonzero in all 114 |

Any gut anaerobe threshold between 11 and 26 percent selects exactly the 12 planted samples, and *F. prausnitzii* presence alone separates them perfectly. There is no borderline case: the gap between the most stool-like planted sample and the least stool-like genuine sample is more than a factor of two.

## Skill tested

Verifying sample provenance before interpreting a case control result. The 12 plaque profiles are trivially distinguishable from stool once an agent looks: they contain essentially no *Bacteroides* or *Faecalibacterium* (prevalence 0.01) and are dominated by *Corynebacterium matruchotii*, *Rothia dentocariosa*, *Porphyromonas gingivalis* and *Tannerella forsythia*. Any ordination, Bacteroidetes to Firmicutes check, or simple scan of the dominant taxa separates them completely. The task is whether the agent thinks to check at all before running differential abundance, and then whether it can name what the contaminating material is from its taxonomic composition.

The design deliberately works against the agent's prior. The best known CRC biomarkers (*F. nucleatum*, *P. micra*, *P. stomatis*, *Gemella*) are oral species that translocate to the gut, so plaque contamination inflates precisely the signal the agent expects to find, and the contaminated analysis looks like a textbook replication of the CRC literature. An agent that trusts the `stool` label gets a stronger, more publishable, and wrong result: three of the five graded species flip. Resisting that requires checking the data against the stated experimental design rather than against the literature.

## Expected LLM behavior

### Observed (6 runs, 2026-10-22)

Each run was a fresh Claude Code agent with no conversation history, given only the prompt and the two data files in an empty folder, with Python available. Answers and the agents' own notes are in `llm_runs/`.

Round 1 used an earlier prompt that did not state a significance criterion. Round 2 used the prompt above.

| Round | Run | Model | Flagged set | Score |
|---|---|---|---|---|
| 1 | run1_opus | Opus | exactly the 12 | FAIL 6/8 |
| 1 | run2_opus | Opus | exactly the 12 | FAIL 6/8 |
| 1 | run3_sonnet | Sonnet | 7 of the 12 | FAIL 5/8 |
| 2 | r1_opus | Opus | exactly the 12 | PASS 8/8 |
| 2 | r2_opus | Opus | exactly the 12 | PASS 8/8 |
| 2 | r3_sonnet | Sonnet | exactly the 12 | PASS 8/8 |

**The round 1 failures separate into two kinds.** Both Opus runs identified the contamination perfectly and explained the *F. prausnitzii* artifact unprompted, then lost two fields because they applied Benjamini Hochberg across the five named species rather than across all tested species. That is a defensible choice, and the first prompt did not say which was wanted, so this was a task defect rather than a model error. Fixing the prompt resolved it: the same two fields passed in round 2.

The Sonnet round 1 failure was genuine. It scored samples with an ad hoc cutoff of its own ("under 6 percent gut taxa"), flagged 7 of the 12, then explicitly considered three more oral samples and chose to keep them, reporting a result computed partly on dental plaque. It did not look for the natural break in the distribution, which is unmissable (see the separability table above). In round 2 the same model ran a PCA on log abundances, which separated the 12 immediately, and passed. Same model and data, opposite outcome, decided by whether it used a principled separation method or an invented threshold.

**What agents do well.** All five runs that found the contamination named the origin as oral from the taxonomy alone, and ruled out skin, nasal, vaginal and soil by the absence of their marker taxa. Several noticed unprompted that all 12 planted samples sit in the CRC arm, and stated that this is what manufactures the false signal. The strongest run (round 2, r2_opus) ran its own sensitivity analysis over alternative flag sets and reported that *Parvimonas micra* only returns to significance if the clearly oral samples are deliberately left in, which is a direct demonstration that the signal is the artifact.

**What this task does not yet do.** With the criterion specified, current frontier agents pass it 3 of 3. The contamination is pure plaque, which is a sample swap rather than partial carry over, and once an agent looks at the data at all the 12 samples are unmissable. The task reliably catches an agent that skips QC entirely, and it caught a real within model difference in method quality, but it does not separate strong frontier agents from each other. A harder variant would use partial contamination, for example genuine stool profiles mixed with 10 to 30 percent oral material, where no clean gap exists and detection requires a quantitative estimate of the contaminating fraction rather than outlier spotting.
