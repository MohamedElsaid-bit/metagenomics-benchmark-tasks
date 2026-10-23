# CRC metagenomics case-control analysis — notes

## Data
- `data/abundance_species.tsv`: 798 species x 126 samples, percent relative abundance; every column sums to 100.000 (max deviation 5e-5), no negatives, no missing samples.
- `data/sample_metadata.tsv`: 126 samples, 65 CRC / 61 control, with age and gender. All recorded as stool.

## Step 1 — QC: 12 samples are not stool
I scored every sample on (a) summed abundance of a panel of obligate gut-core species
(*Faecalibacterium prausnitzii*, *Bacteroides uniformis/vulgatus/dorei/stercoris*,
*Eubacterium rectale/hallii*, *Roseburia intestinalis*, *Blautia wexlerae*,
*Alistipes putredinis*, *Ruminococcus bromii*, *Prevotella copri*,
*Parabacteroides distasonis*, *Bifidobacterium adolescentis*, *Akkermansia muciniphila*, …)
and (b) summed abundance of oral/upper-aerodigestive genera.

The distribution is bimodal with a hard gap: **12 samples contain exactly 0.000% of
every gut-core species**, while the next lowest sample has 0.35% and the rest 13–90%.
Those same 12 samples carry 68–99.7% oral taxa and 0% *Bacteroides*. Richness is also
low (26–95 species vs ~150–250 typical for stool), as expected for a
low-diversity oral/plaque community.

Flagged: **C016, C020, C023, C036, C041, C051, C057, C060, C075, C083, C092, C100**
— all 12 are labelled CRC, which is exactly the pattern that would create spurious
case-control signal.

### Origin = oral (supragingival/subgingival plaque and saliva)
Dominant taxa in these samples are unambiguous oral-cavity organisms, not gut,
skin, vaginal, nasal, soil or reagent-blank taxa:
- Plaque architects / early colonisers: *Rothia dentocariosa*, *Rothia aeria*,
  *Corynebacterium matruchotii*, *Lautropia mirabilis*, *Actinomyces naeslundii*,
  *A. oris*, *A. odontolyticus*, *A. sp. oral taxon 180/448*, *Parascardovia denticolens*,
  *Bifidobacterium dentium*.
- Oral streptococci: *S. gordonii*, *S. sanguinis*, *S. oralis*, *S. cristatus*,
  *S. parasanguinis*, *S. mutans*.
- Oral *Neisseria*/*Kingella*/*Capnocytophaga*: *Neisseria sp. oral taxon 014*,
  *N. elongata*, *N. sicca*, *Kingella denitrificans/oralis*, *Capnocytophaga
  gingivalis/sputigena/ochracea/leadbetteri/granulosa*, *Cardiobacterium hominis*.
- Periodontal pathogens: *Porphyromonas gingivalis*, *P. endodontalis*,
  *Tannerella forsythia*, *Fretibacterium fastidiosum*, *Treponema*-adjacent
  *Bacteroidetes oral taxon 272/274*, *Selenomonas sputigena*, *Eubacterium brachy/yurii*.
- Oral *Prevotella* (*denticola*, *nigrescens*, *oris*, *oralis*, *baroniae*,
  oral taxon 376) rather than gut *P. copri*.
- Marker panels confirm it: skin 0.15% mean, vaginal 0%, nasal 0%, soil 0%,
  gut-core 0% in flagged samples vs 31% gut-core in the rest.

Borderline samples I deliberately **kept**:
- **C074** (control, 63% "oral-genera" score) — this is *Prevotella copri*-dominated
  (62%) with 76% gut-core; *P. copri* is a gut organism, so the genus-level score is
  misleading. Genuine stool.
- **C042** (CRC, gut-core 0.35%) — dysbiotic but unmistakably enteric: *Acidaminococcus
  intestini*, *Hungatella hathewayi*, *Flavonifractor plautii*, *Ruminococcus gnavus*,
  *Clostridium clostridioforme/innocuum*, *Clostridioides difficile*, *Hafnia alvei*.
  A real, abnormal (possibly post-antibiotic) stool sample, not a different body site.

## Step 2 — Differential abundance
Primary test: two-sided Mann–Whitney U on relative abundances, restricted to species
with >=10% overall prevalence (267 species after cleaning), Benjamini–Hochberg FDR.
Cleaned cohort: 53 CRC vs 61 controls.

| species | mean CRC | mean ctrl | prev CRC/ctrl | MWU p | BH q (5 targets) | call |
|---|---|---|---|---|---|---|
| *Fusobacterium nucleatum* | 0.023 | 0.000 | 0.34 / 0.00 | 9.1e-07 | 1.5e-06 | enriched_in_CRC |
| *Peptostreptococcus stomatis* | 0.162 | 0.0002 | 0.47 / 0.08 | 5.3e-07 | 1.3e-06 | enriched_in_CRC |
| *Parvimonas micra* | 0.721 | 0.012 | 0.53 / 0.43 | 0.0061 | 0.010 | enriched_in_CRC |
| *Gemella morbillorum* | 0.087 | 0.0016 | 0.34 / 0.18 | 0.014 | 0.018 | enriched_in_CRC |
| *Faecalibacterium prausnitzii* | 4.41 | 4.95 | 1.00 / 1.00 | 0.151 | 0.151 | not_significant |

Confounder check: age differs between arms (CRC 66.8 vs control 60.6, t-test
p = 0.0036); gender does not (Fisher p = 0.45). Logistic regression of case status on
log10(abundance) adjusted for age and gender gives the same conclusions:
*P. stomatis* p = 2.7e-05 (+), *P. micra* p = 0.0058 (+), *G. morbillorum* p = 0.016 (+),
*F. prausnitzii* p = 0.17 (ns). *F. nucleatum* separates perfectly (0% prevalence in
controls) so the logistic fit is degenerate, but the direction is unambiguous.

Robustness: results are unchanged whether the flag set is 11 (dropping C057), 12, or
13 samples (adding C042) — all four oral-type pathogens stay significantly enriched and
*F. prausnitzii* stays non-significant (p = 0.11–0.20).

## Why the naive analysis is wrong
Including the 12 oral samples makes *F. prausnitzii* look significantly depleted in CRC
(p = 0.0025, q = 0.025, median 2.40 vs 3.94). That is entirely artefactual: all 12
contaminated samples are CRC-labelled and contain exactly 0% *F. prausnitzii*, which
drops its apparent CRC prevalence from 100% to 82%. After exclusion, *F. prausnitzii*
is present in every single sample at comparable levels and shows no significant
difference. The real butyrate-producer depletion signal in this cohort sits in
*Eubacterium hallii/eligens/rectale*, *Anaerostipes hadrus*, *Blautia wexlerae*,
*Roseburia intestinalis* and *Ruminococcus bicirculans* (all q < 0.05 after cleaning),
not in *F. prausnitzii*.

Conversely, the four oral-pathogen markers are genuine: they remain strongly enriched
in CRC after removing the samples that were literally oral communities, so their signal
does not depend on the contamination.

## Files
- `work/explore.py` — table QC, oral/gut genus scoring
- `work/profile.py` — per-sample top taxa and body-site marker panels
- `work/analyze.py` — naive vs cleaned differential abundance, confounder-adjusted models
- `work/check.py` — borderline-sample inspection and flag-set sensitivity analysis
- `work/answer.json` — final answer
