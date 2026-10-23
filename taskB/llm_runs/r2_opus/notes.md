# CRC shotgun metagenome case/control analysis

Data: `data/abundance_species.tsv` (722 species x 114 samples, percent, every column sums to 100.0),
`data/sample_metadata.tsv` (53 CRC / 61 control, all recorded as stool).
Scripts: `work/qc*.py` (QC), `work/final.py` (final analysis). Outputs: `work/da_final.csv`, `work/answer.json`.

## 1. QC / sample screening

No missing samples, no all-zero or duplicated profiles, all columns already compositional (sum 100).

I scored every sample against marker panels for the candidate contamination sources:

| source | panel species found in table | max per-sample load |
|---|---|---|
| skin (Cutibacterium, Staph. epidermidis/hominis, Corynebacterium spp.) | 7 | 1.3 % |
| vaginal (Lactobacillus crispatus/iners/jensenii, Gardnerella, Atopobium vaginae) | 5 | 2.4 % |
| nasal (Moraxella, Dolosigranulum, S. pneumoniae) | 2 | 0.08 % |
| soil (Streptomyces, Bacillus, Arthrobacter, Sphingomonas, ...) | 0 | 0 % |
| reagent blank (Ralstonia, Bradyrhizobium, Burkholderia, Delftia, Methylobacterium, ...) | 4 | 0.1 % |
| **oral / dental plaque** | **156** | **39.3 %** |

Only the oral panel shows anything. The oral score used obligate oral / periodontal taxa that are
essentially absent from healthy stool: *Rothia*, *Neisseria*, *Actinomyces*, *Capnocytophaga*,
*Corynebacterium matruchotii*, *Lautropia*, *Aggregatibacter*, *Leptotrichia*, *Treponema*,
*Tannerella forsythia*, *Porphyromonas gingivalis/endodontalis*, oral *Prevotella*
(intermedia, nigrescens, denticola, melaninogenica, oris, ...), viridans streptococci
(mitis, oralis, sanguinis, gordonii, mutans, vestibularis), *Haemophilus parainfluenzae*,
*Filifactor*, *Selenomonas*, *Mogibacterium*, *Desulfomicrobium orale*, *Dialister pneumosintes*,
*Gemella*, *Granulicatella*, *Abiotrophia*.

The distribution of this score is clearly bimodal:

- 91 samples: 0.005 – 6.1 % (median 0.15 %, 90th pct ~1.2 %) — normal background
- 23 samples: 14.4 – 39.3 % — a separate cluster, with an empty gap from 6.1 % to 14.4 %

Any threshold between ~7 % and ~14 % selects the same 23 samples. Their composition is
unambiguously oral rather than a plausible gut community, e.g.
- P046: *P. gingivalis* 15.8 %, *T. forsythia* 3.5 %, *T. denticola* 2.9 %, *Desulfomicrobium orale* 2.7 %, *P. endodontalis* 2.7 % (the periodontitis "red complex")
- P013: *Neisseria sicca* 18.4 %, *Rothia aeria* 2.3 %, *N. macacae* 2.2 %, *Capnocytophaga sputigena* 0.7 %
- P022: *Corynebacterium matruchotii* 6.2 %, *Actinobaculum* sp. oral taxon 183 5.3 %, *Capnocytophaga granulosa* 4.2 %, *Rothia aeria* 2.9 %
- P037: *Neisseria* sp. oral taxon 014 11.0 %, *N. elongata* 10.9 %, *N. sicca* 6.2 %, *Lautropia mirabilis* 3.0 %

A panel-free check agrees: ranking samples by the total abundance of cohort-rare species
(prevalence < 20 %) puts 21 of these 23 in the top 23 positions, and the oral score is the
dominant correlate of PC2 of the CLR ordination (rho = 0.61, 4.9 % of variance).

**Flagged: 23 samples** (P005, P009, P013, P016, P022, P029, P030, P031, P037, P040, P043, P044,
P046, P047, P049, P056, P062, P094, P099, P107, P110, P111, P114). **Origin: oral.**
22 of the 23 are CRC cases, so the contamination is strongly confounded with case status — which
is exactly why it has to be removed before testing.

## 2. Case/control test

On the 91 retained samples (31 CRC, 60 control):
- species present (> 0) in >= 10 % of the 91 samples -> 266 species tested
- two-sided Mann-Whitney U per species on relative abundance
- Benjamini-Hochberg across all 266 p-values, significance at q < 0.05

6 species pass: *Peptostreptococcus stomatis* (q = 3.3e-9), *Fusobacterium nucleatum* (q = 1.9e-7),
*Dialister pneumosintes* (q = 5.0e-7) and *Parvimonas micra* (q = 4.0e-3) enriched in CRC;
*Blautia wexlerae* (q = 0.014) and *Eubacterium hallii* (q = 0.042) depleted in CRC.

### Requested species (after exclusion)

| species | p | q | mean CRC | mean control | call |
|---|---|---|---|---|---|
| Fusobacterium_nucleatum | 1.4e-9 | 1.9e-7 | 0.031 | 0.000 | enriched_in_CRC |
| Peptostreptococcus_stomatis | 1.2e-11 | 3.3e-9 | 0.191 | 0.0002 | enriched_in_CRC |
| Parvimonas_micra | 6.0e-5 | 4.0e-3 | 1.050 | 0.013 | enriched_in_CRC |
| Rothia_dentocariosa | 0.50 | 0.80 | 0.0013 | 0.0012 | not_significant |
| Eubacterium_rectale | 0.017 | 0.198 | 1.78 | 3.48 | not_significant |
| Streptococcus_salivarius | 0.032 | 0.254 | 0.371 | 0.845 | not_significant |

### Why exclusion changes the answer

With all 114 samples, 30 species pass and *Rothia dentocariosa* looks enriched in CRC
(q = 0.018) while *E. rectale* (q = 0.013) and *S. salivarius* (q = 0.011) look depleted.
All three are artefacts of the contaminated samples: *R. dentocariosa* is a dental-plaque
organism injected with the oral material, and the 15-40 % oral fraction mechanically dilutes
every native gut taxon in those (almost all CRC) samples, so core commensals appear depleted.
Removing the 23 samples removes all three signals. The genuine CRC oral-pathobiont markers
(*F. nucleatum*, *P. stomatis*, *P. micra*) survive and in fact get stronger, confirming the
exclusion is not simply a loss of power.

### Robustness

- Identical calls for all six species with exclusion thresholds of 2 %, 6 % (24 samples) and 14 %
  (23 samples) on the oral score, so the exact cluster boundary does not matter.
- Cases are older than controls (mean 66.8 vs 60.6, p = 0.007). An age- and gender-adjusted
  logistic regression on CLR-transformed abundances agrees on the three enriched taxa and on
  *R. dentocariosa* / *E. rectale* being null; *S. salivarius* is nominally negative there
  (p = 0.003) but does not pass the prescribed BH criterion, so it is reported not_significant.
