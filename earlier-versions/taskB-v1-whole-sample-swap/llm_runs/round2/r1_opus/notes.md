# CRC metagenomics case-control analysis — notes

## Data
- `data/abundance_species.tsv`: 798 species x 126 samples, percent relative abundance (every column sums to 100.000 +/- 2e-5, so no read-depth information is available).
- `data/sample_metadata.tsv`: 126 samples, 65 CRC / 61 control, complete age and gender, all recorded as stool.
- Sample IDs match exactly between the two files; no missing values.

## QC: 12 samples are not stool
I scored every sample for (a) obligate-gut taxa and (b) oral/upper-airway taxa.

Strict core-gut score = sum of unambiguous gut species (`Bacteroides`/`Phocaeicola` spp.,
*Faecalibacterium prausnitzii*, *Eubacterium rectale/hallii*, *Roseburia* spp., *Blautia wexlerae*,
*Anaerostipes hadrus*, *Alistipes* spp., *Akkermansia muciniphila*, *Collinsella aerofaciens*,
*Prevotella copri*, *Escherichia coli*, etc.).

There is a clean, unambiguous gap in this score:

| sample | group | core-gut % |
|---|---|---|
| C023, C020, C051, C041, C092, C100, C075 | CRC | 0.000 |
| C057 | CRC | 0.008 |
| C036 | CRC | 0.011 |
| C060 | CRC | 0.022 |
| C016 | CRC | 0.095 |
| C083 | CRC | 1.243 |
| **gap** | | |
| C042 (next sample) | CRC | 21.4 |
| cohort median | | 66.4 |

Those same 12 samples carry 58-99 % obligate oral taxa, and their dominant species are
classic supragingival/subgingival plaque and saliva organisms, e.g.:
- C020: *Streptococcus gordonii* 22 %, *S. oralis* 17 %, *Neisseria sicca* 16 %, *Rothia dentocariosa* 11 %, *Kingella denitrificans* 9 %, *Lautropia mirabilis*
- C100: *Rothia dentocariosa* 21 %, *R. aeria* 9 %, *Actinomyces massiliensis/naeslundii*, *Corynebacterium matruchotii* 8 %
- C023: *Streptococcus sanguinis* 49 %, *Neisseria* oral taxon 014 38 %, four *Capnocytophaga* spp.
- C075 / C060 / C051: *Porphyromonas gingivalis*, *Tannerella forsythia*, *Fretibacterium fastidiosum*, *Eubacterium brachy/yurii* (periodontal consortium)
- C016 / C036 / C041: oral *Prevotella* (*denticola*, *nigrescens*, *oris*, *oralis*), *Atopobium rimae*, *Olsenella uli*, *Parascardovia denticolens*

This is the oral cavity, not skin (no *Cutibacterium/Staphylococcus* dominance), not vaginal
(no *Lactobacillus*), not nasal (no *Staphylococcus/Dolosigranulum/Moraxella*), not soil
(no environmental Proteobacteria/Actinomycetota) and not a reagent blank (rich, coherent,
host-associated communities with 26-95 species each). => `flagged_sample_origin = "oral"`.

Critically, **all 12 are labelled CRC**, so leaving them in inflates every oral-derived taxon
in the case group — exactly the confound that would create spurious "CRC-enriched" oral species.

Other QC checks that came up clean: no duplicated/replicate columns (max pairwise Pearson r on
log abundance = 0.78, the one 0.96 raw-scale pair C074/C095 is just shared *P. copri* dominance),
no all-zero samples, no metadata inconsistencies. C042 has an atypical but genuinely intestinal
community (*Acidaminococcus intestini*, *Hungatella hathewayi*, *Flavonifractor plautii*,
*Clostridioides difficile*) and was kept.

## Differential abundance (the comparison I consider correct)
- Cohort: 114 samples (53 CRC, 61 control) after dropping the 12 oral samples.
- Prevalence filter: species non-zero in >= 10 % of the 114 included samples -> 267 species tested.
- Test: two-sided Mann-Whitney U per species (relative abundance is zero-inflated and non-normal).
- Multiple testing: Benjamini-Hochberg across all 267 tested species, significance at q < 0.05.
- 12 species reach q < 0.05. Full tables in `work/results_CLEAN_excl_oral.csv` (and
  `work/results_ALL_samples.csv` for the un-cleaned comparison).

Target species (clean cohort):

| species | mean % CRC | mean % control | prev CRC | prev ctrl | p | q | call |
|---|---|---|---|---|---|---|---|
| *Fusobacterium nucleatum* | 0.0228 | 0.000 | 0.34 | 0.00 | 9.1e-07 | 8.1e-05 | enriched_in_CRC |
| *Peptostreptococcus stomatis* | 0.162 | 0.00019 | 0.47 | 0.08 | 5.3e-07 | 8.1e-05 | enriched_in_CRC |
| *Parvimonas micra* | 0.721 | 0.0123 | 0.53 | 0.43 | 0.0061 | 0.085 | not_significant |
| *Gemella morbillorum* | 0.0875 | 0.0016 | 0.34 | 0.18 | 0.0143 | 0.141 | not_significant |
| *Faecalibacterium prausnitzii* | 4.41 | 4.95 | 1.00 | 1.00 | 0.151 | 0.500 | not_significant |

Other notable hits in the clean cohort: *Dialister pneumosintes* and
*Porphyromonas asaccharolytica* up in CRC; *Eubacterium hallii*, *E. eligens*, *E. rectale*,
*Anaerostipes hadrus*, *Roseburia intestinalis*, *Blautia wexlerae*, *Ruminococcus bicirculans*,
*Bifidobacterium catenulatum*, *Streptococcus salivarius* down in CRC — i.e. the expected loss of
butyrate producers, but *F. prausnitzii* itself is not individually significant here.

## Why the cleaning matters
With all 126 samples, 32 species pass q < 0.05 and all five target species come out "significant"
(*Parvimonas* q=0.021, *Gemella* q=0.031, *F. prausnitzii* q=0.025 depleted). Those three calls are
artefacts of the 12 CRC-labelled oral samples: *P. micra* and *G. morbillorum* are oral commensals
that are abundant in plaque, and *F. prausnitzii* is 0 % in all 12 oral samples, which manufactures
an apparent depletion in CRC. *F. nucleatum* and *P. stomatis* survive the cleaning and are the
genuine signal.

## Robustness
Conclusions are stable across the analysis choices tested (see `work/sens.py`):

| flag set | test | F. nucleatum q | P. stomatis q | P. micra q | G. morbillorum q | F. prausnitzii q |
|---|---|---|---|---|---|---|
| excl 12 | Mann-Whitney | 0.0001 | 0.0001 | 0.085 | 0.141 | 0.500 |
| excl 12 | CLR + t-test | <0.0001 | <0.0001 | 0.074 | 0.087 | 0.476 |
| excl 12 | logistic, age+gender adj. | n/a* | 0.0032 | 0.159 | 0.263 | 0.583 |
| excl 11 (keep C083) | Mann-Whitney | 0.0001 | 0.0001 | 0.113 | 0.159 | 0.450 |
| excl 13 (also C042) | Mann-Whitney | 0.0001 | 0.0001 | 0.067 | 0.129 | 0.531 |

\* the logistic model fails to converge for *F. nucleatum* because it is perfectly separating
(0 % prevalence in controls); the Mann-Whitney result stands.

Cases are older than controls in the clean cohort (66.8 vs 60.6 y, t-test p = 0.004; gender
balanced, p = 0.45), but the age+gender-adjusted model gives the same five calls, so age
confounding does not change the answer.

## Answer
```json
{
  "n_flagged_samples": 12,
  "flagged_sample_ids": ["C016","C020","C023","C036","C041","C051","C057","C060","C075","C083","C092","C100"],
  "flagged_sample_origin": "oral",
  "crc_vs_control": {
    "Fusobacterium_nucleatum": "enriched_in_CRC",
    "Peptostreptococcus_stomatis": "enriched_in_CRC",
    "Parvimonas_micra": "not_significant",
    "Gemella_morbillorum": "not_significant",
    "Faecalibacterium_prausnitzii": "not_significant"
  }
}
```

## Files
- `work/explore.py`, `work/explore2.py`, `work/explore3.py`, `work/explore4.py` — QC / contamination screen
- `work/da_analysis.py` — main differential abundance analysis, writes `work/answer.json`
- `work/sens.py` — sensitivity analysis
- `work/results_CLEAN_excl_oral.csv`, `work/results_ALL_samples.csv` — per-species statistics
- `work/sample_profile.csv` — per-sample oral/gut scores
