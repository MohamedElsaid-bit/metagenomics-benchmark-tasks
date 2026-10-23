# Two benchmark tasks for long horizon metagenomics

Mohamed Elsaid

I wrote two tasks, both built from open access shotgun metagenomic data, both graded deterministically, and both tested against live agents before submission. The notes below report what the testing actually showed, including the places where the first version of a task turned out to be wrong.

| | Task A | Task B |
|---|---|---|
| ID | `mgx_drug_coprescription_confounding` | `mgx_crc_partial_contamination` |
| Skill | separating a drug effect from the prescribing pattern it sits in | detecting partial sample contamination that is correlated with case status |
| Data | MetaCardis cardiometabolic cohort, 1,831 stool metagenomes | French CRC cohort, 114 stool metagenomes, 21 of them partially contaminated with real dental plaque |
| Graded fields | 8 | 10 |
| Tested | 6 runs over 2 prompt versions | 6 runs over 2 prompt versions |
| Result | 6 of 6 pass | 2 of 3 pass on the final version, both Opus runs pass, Sonnet fails |
| Role | calibration item | discriminating item |

Full specifications are in `taskA/TASK.md` and `taskB/TASK.md`. Each folder holds the exact files handed to the agent, the answer key, the grader, and the raw answers and notes from every test run.

## Task A: which drug effects are real?

Within the 792 type 2 diabetes participants of a European cardiometabolic cohort, the agent must decide, for each of five medications, whether its association with species composition is its own or is inherited from a drug it is co-prescribed with.

The answer is metformin and proton pump inhibitors independent, sulfonylurea, statin and insulin not. Statins are the sharpest item: the most prescribed drug in the group at 440 of 792 participants, 38 associated species when tested alone, and zero under every adjusted and stratified route. Testing one drug at a time gives four of the five labels wrong. The prompt never mentions co-prescription, and the labels are the neutral `independent` and `not_independent`, so the mechanism has to be discovered rather than applied.

A fourth field asks which single variable accounts for the apparent sulfonylurea association. The answer is uniquely metformin: adjusting for it alone drops 48 species to 2, while adjusting for DPP4 inhibitors, aspirin, PPI, statin, insulin, country, age or BMI each leaves 32 or more. I did not ask the same question about statins, because several different adjustments remove that signal, so there is no unique answer.

The task also carries a deterministic numeric anchor, the fraction of statin users who also take metformin (0.409), which cannot be recalled from any publication, and three traps in the data itself. The medication token `su` is a substring of `insulin`, so drugs must be matched as exact tokens after splitting on semicolons. Age and BMI are missing for all 199 Danish participants, so adding those covariates silently drops an entire recruitment site and takes the analysis set from 792 to 588. And the statin confounding runs in a non obvious direction: statin users are *less* likely to be on metformin than non users (40.9 against 65.3 percent), so the marginal statin signal is partly an inverse shadow of the metformin signature rather than a shared one.

## Task B: a contaminated cohort that confirms the literature

The agent is given a colorectal cancer case control study in which 21 stool profiles have been mixed with real subgingival plaque at 26 to 40 percent. Twenty of the 21 affected samples are cases. Nothing is simulated; both components are genuine human metagenomes.

The design works against the agent's prior in two directions at once. The best known CRC biomarkers (*Fusobacterium nucleatum*, *Parvimonas micra*, *Peptostreptococcus stomatis*) are oral organisms that genuinely translocate to the gut, so oral contamination inflates exactly the signal the agent expects to find. At the same time it dilutes ordinary gut commensals in the case arm, manufacturing an apparent depletion of butyrate producers (*Eubacterium rectale*, *Streptococcus salivarius*), which is another textbook CRC claim. An agent that trusts the data reports two findings that are pure artefact. An agent that discards everything oral looking loses the three genuine markers.

Three of the six graded species flip between the contaminated and cleaned analyses. The three true markers survive cleaning by construction, because the contaminated cases were chosen from participants who carry neither *F. nucleatum* nor *P. stomatis*.

## What the testing showed

Every task was run against fresh agents with no conversation history, given only the prompt and the data files in an empty working directory, with Python available and no web access. Eighteen runs were graded in total: nine across three versions of the confounding task and nine across three versions of the contamination task.

### Confounding

Nine runs across three versions of this task produced no failures: a T2D versus healthy design, the drug by drug design with an explicit hint that the drugs are co-prescribed, and the submitted version with every hint removed, including the label name that had previously given away the mechanism. Removing the hints changed nothing. Every run on the final version independently decided to adjust the drugs against each other without being told that co-prescription was the issue, and every run located the variable behind the sulfonylurea signal with a single covariate scan. The models parse free text medication fields unprompted, adjust drugs against each other rather than stratifying naively, resist the literature prior that butyrate producer depletion is intrinsic to diabetes, and correctly decline to reproduce the mouse finding that metformin raises *Akkermansia*, which this human cohort does not support.

### Data provenance

The partial contamination task split the models cleanly, and the split fell on the quality control step rather than the statistics. All three runs got every species call right. Both Opus runs built a panel of obligate oral taxa, found the bimodal split in the distribution, and recovered all 21 contaminated samples. One of them excluded *F. nucleatum*, *P. micra* and *P. stomatis* from its own detection panel on the grounds that those organisms genuinely colonize the CRC gut and using them to detect contamination would be circular. That is the correct reasoning and it is what preserves the three true markers.

The Sonnet run scored whole genera instead of oral specific species. That sweeps in *Prevotella copri*, *Streptococcus thermophilus* and *Veillonella*, which are ordinary gut residents, so its score was noisy at the boundary. It described the 7 to 9 percent range as ambiguous, kept samples there, and missed two genuinely contaminated samples. It still reached the right species calls, because removing 19 of 21 is enough for the artefacts to disappear, but it reported a comparison that still contained contaminated material.

### Where the difficulty actually sits

Every run that found the contamination then handled the statistics correctly, so the detection step is what the task is really measuring. An earlier version of Task B used whole mislabelled plaque samples rather than partial mixtures. Those samples contain no *Bacteroides* and no *Faecalibacterium* at all, fall out of any ordination immediately, and both models passed once the significance criterion was specified. Replacing the sample swap with partial carry over removed every shortcut: in the submitted version, PCA on CLR values does not separate the contaminated samples on any of the first four components, and gut core content overlaps genuine stool completely. The only route left is to know which organisms do not live in the gut. That version is kept in `earlier-versions/taskB-v1-whole-sample-swap/` for comparison; the progression from it to Task B is the clearest evidence in this submission of how the difficulty was tuned.

## How the tasks were made robust

Latch's brief asks for tasks that are robust to multiple reasonable ways of doing the question. Testing turned up three defects in my own tasks. All three are fixed, and all three are documented in full in the task files.

**An unspecified statistical convention.** The first version of the contamination task asked which species differ without saying how to correct for multiple testing. Two Opus runs identified the contamination perfectly and then lost two fields because they applied Benjamini Hochberg across the five named species rather than across all 267 tested. Both choices are defensible and the difference changes the answer: *Parvimonas micra* sits at q = 0.010 one way and q = 0.085 the other. Every submitted prompt now fixes the prevalence filter, the correction scope and the significance threshold.

**A stale answer key.** The contamination dataset was rebuilt with different parameters, which regenerated both the contamination assignment and the anonymized sample IDs, but the answer key had been written by hand from an earlier console log. The first grading pass scored all three runs against the wrong sample list and reported three failures that were not real. It was caught because a model flagged a sample the key called clean while the data showed *Rothia dentocariosa* at 14.7 percent in it, which could not both be true. The key is now emitted by the same script that writes the data, so the two cannot diverge. An answer key transcribed by hand is a latent defect in any benchmark task.

Beyond those fixes, every graded call in both tasks was computed under multiple reasonable methods and any item whose answer depended on the method was removed from grading. Task A dropped DPP4 inhibitors and aspirin, whose labels flipped with country adjustment or stratification. Task B dropped the count of significant species, which ranges from 5 to 12 depending on method, and dropped *Gemella morbillorum*, which was unstable across prevalence filters.

**A grading threshold that decided a result.** Task B grades its flagged sample list asymmetrically rather than as an exact match, because the two errors are not scientifically equivalent: leaving contamination in biases the comparison, since the contamination is concentrated in the case arm, while dropping a clean sample only costs power. There is also a legitimate range of marker definitions, as two uncontaminated samples carry periodontal organisms (*Prevotella intermedia* at 11.3 percent, *Treponema berlinense* at 7.3 percent) that a broader oral panel would reasonably flag.

My first version of that rule capped the extras at 3, which was a mistake of the same kind as the two above. One test run flagged exactly 3 extra samples and another flagged 2, so moving the cap by one would have flipped a pass to a fail: the threshold, not the analysis, was deciding the outcome. The rule now requires full recall and sets the guard on extras at 18, which is 20 percent of the clean samples. It exists only to reject degenerate answers such as flagging half the cohort to guarantee recall, and is far from binding for any real analysis. Over flagging within the guard is judged where it belongs, in the species calls, since an exclusion set that distorts the comparison fails those fields anyway. This follows the endpoint grading philosophy of scBench-Long, which does not require intermediate analyses to follow a prescribed workflow and lets an agent reach the target surface by any valid path.

**Three traps in the data were found by the agents rather than by me.** Test runs surfaced properties of the data that no amount of reading the task specification would have revealed: the medication token `su` is a substring of `insulin`, so naive string matching silently merges two of the drug groups being compared; age and BMI are missing for every Danish participant, so adjusting for them drops a whole recruitment site; and the statin label is stable under rank, arcsine and CLR transforms but retains 7 species under a plain log10 transform, which is why it is verified by a transform free stratified test instead. Each was found because an agent explained its own reasoning. Testing tasks against live agents is worth doing for what it reveals about the data, not only for the pass rate.

## Data provenance

All data is open access and redistributable, taken from curatedMetagenomicData (Bioconductor; Pasolli et al. 2017, Nature Methods, doi:10.1038/nmeth.4468), which provides uniformly processed MetaPhlAn3 profiles. No controlled access resource is used.

| Study | Role | Citation |
|---|---|---|
| MetaCardis_2020_a | Task A, whole cohort | Fromentin et al. 2022, Nature Medicine, PMID 34880489 |
| ZellerG_2014 | Task B, CRC stool base cohort | Zeller et al. 2014, Mol Syst Biol, doi:10.15252/msb.20145645 |
| GhensiP_2019 | Task B, oral material mixed in | Ghensi et al. 2020, NPJ Biofilms Microbiomes, doi:10.1038/s41522-020-00155-7 |

Identifying columns were removed and sample IDs anonymized and shuffled in both tasks, so the source cohorts cannot be looked up from the files given to the agent. This matters: in early testing on the unanonymized Task A bundle, one agent recognized the cohort and named the source paper from the original sample ID prefix.

## Repository layout

```
taskA/     Task A: prompt, dataset, answer key, grader, 3 test runs
taskB/     Task B: prompt, dataset, answer key, grader, 3 test runs
earlier-versions/
  taskB-v1-whole-sample-swap/   first version of Task B, whole mislabelled samples, 6 test runs
  taskA-v1-t2d-vs-healthy/      first version of Task A, 3 test runs
results/   robustness output tables
scripts/   every build and analysis script; all numbers in the task files are reproducible from these
```

Every figure quoted here was computed from the data by the scripts in `scripts/`. None of it is carried over from a publication.
