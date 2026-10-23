# CRC metagenomics case-control analysis — notes

## Data
- `data/abundance_species.tsv`: 798 species x 126 samples, each column sums to 100% (relative abundance).
- `data/sample_metadata.tsv`: 126 samples, 65 CRC / 61 control, with age and gender. Sample IDs match exactly between the two files; no missing values.

## QC: 12 samples are not stool
All samples are *recorded* as stool, but 12 have a composition incompatible with faecal material. I scored every sample on
(a) the summed abundance of obligate gut taxa (Bacteroides/Phocaeicola, Faecalibacterium, Roseburia, Blautia, Alistipes,
Ruminococcus, Eubacterium, Parabacteroides, Bifidobacterium longum/adolescentis, Akkermansia, Clostridium, etc.) and
(b) the summed abundance of oral/dental taxa (Streptococcus, Neisseria, Rothia, Actinomyces, Capnocytophaga, Veillonella,
Porphyromonas, Prevotella denticola/nigrescens/oris, Corynebacterium matruchotii, Lautropia, Atopobium, Leptotrichia,
Tannerella, Kingella, Cardiobacterium, plus all `*_oral_taxon_*` species).

The gut-taxa score separates with a large, unambiguous gap:

| sample | core gut % | group |
|---|---|---|
| C020, C041, C100 | 0.00 | CRC |
| C023 | 0.01 | CRC |
| C092 | 0.01 | CRC |
| C016 | 0.48 | CRC |
| C060 | 0.66 | CRC |
| C036 | 1.17 | CRC |
| C083 | 1.24 | CRC |
| C075 | 1.92 | CRC |
| C057 | 2.79 | CRC |
| C051 | 11.08 | CRC |
| **gap** | | |
| next lowest (C042) | 39.4 | CRC |

All remaining 114 samples have >= 39% core gut taxa. The 12 outliers are instead dominated by oral/supragingival-plaque
organisms, e.g. C023 = *Streptococcus sanguinis* 48% + *Neisseria* sp. oral taxon 014 38% + *Capnocytophaga* spp.;
C100 = *Rothia dentocariosa* 22% + *S. gordonii* 11% + *Rothia aeria* + *Corynebacterium matruchotii* + *Lautropia mirabilis*;
C092 = *Actinomyces naeslundii* 21% + *C. matruchotii* 15% + *Capnocytophaga ochracea*. Several also carry near-zero
richness (26-29 species) consistent with a low-diversity plaque/saliva community rather than stool.

**Origin: `oral`.** Not skin (no *Cutibacterium*/*Staphylococcus* dominance), not vaginal (no *Lactobacillus*), not nasal
(no *Corynebacterium propinquum*/*Dolosigranulum*/*S. aureus*), not soil or a reagent blank (these are a coherent
human-commensal community, not environmental taxa or sparse noise).

Critically, **all 12 are CRC cases** — this is a group-correlated technical artefact that will inflate any oral-taxon
enrichment in CRC. They are therefore excluded from the case-control comparison (n = 114: 53 CRC, 61 control).

C107 and C013 were checked and kept: both are plainly faecal (*Bacteroides dorei*, *Ruminococcus bromii*, *Blautia*),
their elevated *Streptococcus* is *S. thermophilus*, a dairy/food organism normal in stool.

## Differential abundance
- Test set: species present (>0) in >= 10% of the 114 included samples -> 267 species tested.
- Test: two-sided Mann-Whitney U (Wilcoxon rank-sum) on relative abundances — appropriate for zero-inflated,
  non-normal compositional data.
- Multiple testing: Benjamini-Hochberg across all 267 tested species; significance at q < 0.05.
- 12 species significant. Enriched in CRC: *Peptostreptococcus stomatis*, *Fusobacterium nucleatum*,
  *Dialister pneumosintes*, *Porphyromonas asaccharolytica*. Depleted in CRC: *Eubacterium hallii*, *E. eligens*,
  *Anaerostipes hadrus*, *Ruminococcus bicirculans*, *Blautia wexlerae*, *Roseburia intestinalis*,
  *Bifidobacterium catenulatum*, *Streptococcus salivarius*.

### The five requested species (114 samples, correct analysis)

| species | median CRC | median control | mean CRC | mean control | q | call |
|---|---|---|---|---|---|---|
| Fusobacterium_nucleatum | 0 | 0 | 0.023 | 0.000 | 8.1e-5 | enriched_in_CRC |
| Peptostreptococcus_stomatis | 0 | 0 | 0.162 | 0.000 | 8.1e-5 | enriched_in_CRC |
| Parvimonas_micra | 0.006 | 0 | 0.721 | 0.012 | 0.085 | not_significant |
| Gemella_morbillorum | 0 | 0 | 0.087 | 0.002 | 0.141 | not_significant |
| Faecalibacterium_prausnitzii | 3.14 | 3.94 | 4.41 | 4.95 | 0.500 | not_significant |

### Why the QC step matters
With all 126 samples the same pipeline calls 32 species significant and **all five** focus species significant
(F. nucleatum q=3e-6, P. stomatis q=6e-4, P. micra q=0.021, G. morbillorum q=0.031, F. prausnitzii depleted q=0.025).
The 12 oral CRC samples alone manufacture the *P. micra*, *G. morbillorum* and *F. prausnitzii* results: *P. micra*
reaches 9.6% and *F. nucleatum* 18.7% in C051 while *F. prausnitzii* is ~0 in all 12 oral samples, which drags the
CRC "gut" distribution down. After removal only the two genuinely faecal CRC markers survive.

## Sensitivity
- Flag-set variants give the same calls for all five species: keeping borderline C051 (11 excluded) -> P. micra q=0.070,
  G. morbillorum q=0.107; additionally excluding C107 (13) -> q=0.079 / 0.111. Only an under-removal that keeps the
  clearly oral C041/C051/C057 (9 excluded) pushes P. micra back to q=0.031, which is exactly the contamination artefact.
- Age is imbalanced after exclusion (CRC mean 66.8 vs control 60.6, p=0.007); gender is balanced. An age- and
  gender-adjusted logistic regression on log10 abundance, BH-corrected over the same 267 species, agrees: only
  *P. stomatis* is significant (q=0.007), with P. micra q=0.25, G. morbillorum q=0.26, F. prausnitzii q=0.56.
  (*F. nucleatum* is degenerate in that model — it is zero in every control, so the logistic fit separates perfectly;
  the rank-sum result is the valid one for it.) So the two positive calls are not driven by age.

## Files
- `work/analyze.py` — main analysis, writes `work/answer.json`
- `work/explore.py`, `work/explore2.py`, `work/explore3.py` — QC / contamination screening
- `work/sens.py` — sensitivity to flag set and covariate adjustment
- `work/results_filtered.tsv` — full per-species results, 114 samples (the correct analysis)
- `work/results_naive.tsv` — full per-species results, all 126 samples (for contrast)
