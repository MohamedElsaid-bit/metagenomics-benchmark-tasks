# CRC stool metagenome case-control analysis

Data: `data/abundance_species.tsv` — 798 species x 126 samples, percent relative abundance
(every column sums to 100.000 +/- 2e-5, so the table is already closed/complete).
`data/sample_metadata.tsv` — 126 samples, 65 CRC / 61 control, no missing age or gender.
Sample IDs match the abundance columns exactly.

## 1. QC: 12 samples are not stool

All samples are *recorded* as stool, but a per-sample body-site check shows a discrete,
non-overlapping group of 12 that cannot be faecal material.

Scoring: summed relative abundance of species that are essentially oral-cavity-exclusive
(genera *Rothia, Corynebacterium (matruchotii), Lautropia, Capnocytophaga, Kingella,
Neisseria, Actinomyces, Actinobaculum, Cardiobacterium, Tannerella, Treponema,
Fretibacterium, Anaeroglobus, Selenomonas, Olsenella, Porphyromonas, Atopobium,
Parascardovia, Slackia, Pseudopropionibacterium*, every `*_oral_taxon_*` species,
plus *Eubacterium brachy/yurii*, *Prevotella denticola/nigrescens/oris/oralis*,
*Streptococcus sanguinis/gordonii/oralis/cristatus*, *Bifidobacterium dentium*)
versus the summed abundance of obligate distal-gut taxa.

| sample | oral-only % | gut-only % | *Bacteroides*+*Phocaeicola* | *F. prausnitzii* | *E. coli* | richness |
|---|---|---|---|---|---|---|
| C023 | 99.0 | 0.2 | 0 | 0 | 0 | 66 |
| C083 | 95.7 | 1.3 | 0 | 0 | 0 | 84 |
| C100 | 94.5 | 0.0 | 0 | 0 | 0 | 29 |
| C075 | 93.0 | 1.4 | 0 | 0 | 0 | 87 |
| C060 | 91.2 | 1.0 | 0.02 | 0 | 0 | 54 |
| C092 | 90.1 | 0.1 | 0 | 0 | 0 | 51 |
| C020 | 85.1 | 0.4 | 0 | 0 | 0 | 26 |
| C016 | 84.5 | 1.6 | 0 | 0 | 0 | 62 |
| C036 | 80.5 | 0.6 | 0 | 0 | 0 | 93 |
| C041 | 70.0 | 1.8 | 0 | 0 | 0 | 27 |
| C051 | 61.5 | 0.7 | 0 | 0 | 0 | 84 |
| C057 | 57.9 | 9.4 | 0.008 | 0 | 0 | 95 |
| **next-highest sample (C034)** | **13.7** | **76.1** | 8.9 | 11.1 | 0.08 | 179 |

The gap is unambiguous: 57.9% -> 13.7%. Decision rule used: oral-only > 50% **and**
gut-only < 10%.

Why **oral** and not another site:
- Dominant taxa are dental-plaque / subgingival organisms: *Rothia dentocariosa* and
  *R. aeria*, *Corynebacterium matruchotii*, *Lautropia mirabilis*, *Actinomyces
  naeslundii/oris/odontolyticus/massiliensis* and `Actinomyces_sp_oral_taxon_180/448`,
  *Capnocytophaga gingivalis/ochracea/sputigena*, *Kingella oralis/denitrificans*,
  *Neisseria sicca/elongata* and `Neisseria_sp_oral_taxon_014`, *Streptococcus
  sanguinis/gordonii/oralis*, *Porphyromonas gingivalis/endodontalis*,
  *Tannerella forsythia*, *Prevotella denticola/nigrescens/oris*, *Atopobium rimae*,
  *Olsenella uli*, *Parascardovia denticolens*, *Veillonella parvula*.
- Simultaneously **zero** for every obligate gut marker: *Bacteroides*/*Phocaeicola*,
  *Faecalibacterium prausnitzii*, *Escherichia coli*, *Blautia*, *Roseburia*,
  *Alistipes*, *Bifidobacterium longum/adolescentis*, *Methanobrevibacter smithii*.
- Not skin (no *Cutibacterium acnes* / *Staphylococcus epidermidis* dominance — *C. acnes*
  appears only as a trace in C020), not vaginal (no *Lactobacillus crispatus/iners*,
  no *Gardnerella*), not nasal (no *Corynebacterium/Staphylococcus/Dolosigranulum*
  nasal pattern; the *Corynebacterium* present is *C. matruchotii*, a plaque species),
  not soil (no environmental Proteobacteria/Actinobacteria), not a reagent blank
  (rich 26-95 species communities, not *Burkholderia/Ralstonia/Pseudomonas* kit taxa).

**Critically, all 12 are in the CRC arm** (12/65 cases, 0/61 controls). Leaving them in
confounds group with body site and spuriously inflates every oral taxon in CRC while
deflating every gut commensal.

Borderline samples I deliberately **kept**:
- **C042** (CRC) — low in my narrow gut panel but a genuinely dysbiotic gut community
  (*Acidaminococcus intestini*, *Hungatella hathewayi*, *Flavonifractor plautii*,
  *Ruminococcus gnavus*, *Clostridioides difficile*, *Hafnia alvei*); no oral taxa at all.
- **C074 / C095 / C071** — high in a naive "oral" score only because *Prevotella copri*
  (a gut species) dominates them; all have normal *Bacteroides* and *F. prausnitzii*.
- **C051** needed care in the other direction: it scored ~11% "gut" only because
  *Eubacterium brachy* and *E. yurii* are in the *Eubacterium* genus — both are
  periodontal species. Its profile (*F. nucleatum* 18.7%, `Bacteroidetes_oral_taxon_274`
  16.7%, *Parvimonas micra* 9.6%, *Tannerella forsythia*, *Solobacterium moorei*, zero
  *Bacteroides*) is subgingival plaque, so it is flagged.

## 2. Differential abundance

CRC vs control after exclusion: n = 53 CRC, 61 control. Species filtered to those present
in >= 10% of retained samples (267 species). Age differs between arms
(66.8 vs 60.6 y, p = 0.004); gender does not (p = 0.45), so age-adjusted models were run
as a sensitivity check.

Three concordant analyses on the clean set (BH-FDR across the five pre-specified species):

| species | mean CRC | mean ctl | prev CRC | prev ctl | q (Wilcoxon, RA) | q (Wilcoxon, CLR) | q (logistic, CLR + age + sex) | call |
|---|---|---|---|---|---|---|---|---|
| *Fusobacterium nucleatum* | 0.0228 | 0.0000 | 18/53 | 0/61 | 2e-6 | 3.4e-4 | 0.077* | enriched_in_CRC |
| *Peptostreptococcus stomatis* | 0.1623 | 0.0002 | 25/53 | 5/61 | 2e-6 | 3.8e-5 | 2.0e-4 | enriched_in_CRC |
| *Parvimonas micra* | 0.7206 | 0.0123 | 28/53 | 26/61 | 0.010 | 0.0026 | 0.0054 | enriched_in_CRC |
| *Gemella morbillorum* | 0.0875 | 0.0016 | 18/53 | 11/61 | 0.018 | 0.0075 | 0.018 | enriched_in_CRC |
| *Faecalibacterium prausnitzii* | 4.406 | 4.955 | 53/53 | 61/61 | 0.151 | 0.096 | 0.156 | **not_significant** |

\* the regression coefficient for *F. nucleatum* is unstable because of perfect separation
(0/61 controls carry it); Fisher's exact test on prevalence gives p = 1.6e-7, so the
enrichment is the strongest in the panel, not a null.

*P. micra* and *G. morbillorum* carry the same signature in both: prevalence is similar
between arms but the carriage level is ~55-60x higher in cases (*P. micra* among carriers:
median 0.078 vs 0.006, p = 1.1e-6; *G. morbillorum*: 0.054 vs 0.004, p = 0.0011). Their raw
Wilcoxon p-values (0.006, 0.014) survive correction across the five pre-specified targets
and under CLR/log transforms, which are better suited than raw relative abundance to this
zero-inflated compositional data; only BH across all 267 species pushes them to
q = 0.085/0.141. Taken together the enrichment call is well supported.

## 3. The *F. prausnitzii* trap

With all 126 samples included, the naive comparison reports *F. prausnitzii* as
significantly **depleted** in CRC (3.59 vs 4.95%, q = 0.025). That entire signal comes from
the 12 oral samples, which are structurally zero for it. Among real stool samples the
difference is 4.41 vs 4.95% (q = 0.15) — not significant. The same contamination inflates
*F. nucleatum*'s CRC mean 24-fold (0.555 -> 0.023) and roughly triples the apparent number
of significant species (32 -> 12), mostly spurious depletions of gut commensals
(*Eubacterium rectale*, *Ruminococcus bromii*, *Anaerostipes hadrus*, *Streptococcus
thermophilus*, ...).

The four oral-pathobiont enrichments (*F. nucleatum*, *P. stomatis*, *P. micra*,
*G. morbillorum*) are robust to the exclusion and match the established CRC faecal
signature; the *F. prausnitzii* depletion is an artefact of the mislabelled samples.

## Files
- `explore.py`, `explore2.py`, `explore3.py`, `explore4.py` — body-site QC exploration
- `analysis.py` — QC rule + differential abundance, naive vs clean
- `borderline.py` — per-value inspection and method comparison for the borderline species
- `confirm.py` — final CLR / age-adjusted confirmation, writes `answer.json`
- `qc_table.tsv`, `sample_profile*.tsv`, `da_all.tsv`, `da_clean.tsv` — intermediates
