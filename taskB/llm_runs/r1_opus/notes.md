# CRC metagenomics case/control analysis — notes

## Data
- `data/abundance_species.tsv`: 722 species x 114 samples, percent relative abundance (every column sums to 100.0).
- `data/sample_metadata.tsv`: 114 samples, 53 CRC / 61 control. Sample IDs match the abundance table exactly (no missing/extra).
- CRC cases are older (mean 66.8 vs 60.6 y, MWU p = 0.007); gender roughly balanced. No covariate adjustment was used for the reported calls (see "Robustness").

## QC: finding unsuitable samples
Checked all plausible non-stool / contaminant origins by summing marker taxa per sample:

| candidate origin | markers checked | max per-sample total |
|---|---|---|
| vaginal | *Lactobacillus crispatus/gasseri*, *Gardnerella vaginalis*, *Prevotella bivia* | 2.4 % (one sample), else <0.15 % |
| skin | *Cutibacterium acnes*, *Brevibacterium*, *Staph. epidermidis* | 1.3 % |
| nasal | *Staph. aureus*, *S. pneumoniae*, *Moraxella* | 0.08 % |
| soil | *Pseudomonas*, *Bacillus*, *Streptomyces*, *Bradyrhizobium*… | 0.09 % |
| reagent blank | *Ralstonia*, *Burkholderia*, *Delftia*, *Acinetobacter*, *Comamonas*, *Serratia* | 0.05 % |
| **oral / saliva** | obligate oral + periodontal taxa (see below) | **39.6 %** |

Only the oral signature is anomalous, so `flagged_sample_origin = oral`.

Also checked and found nothing: no duplicate/replicate columns (minimum pairwise Bray–Curtis = 0.36), no implausible library compositions beyond the oral issue, no samples with zero gut core.

### Oral score definition
Summed the relative abundance of species that are obligate oral-cavity / periodontal organisms and are **not** gut colonisers: all *Rothia*, *Neisseria*, *Haemophilus*, *Actinomyces*, *Capnocytophaga*, *Leptotrichia*, *Aggregatibacter*, *Granulicatella*, *Lautropia*, *Treponema*, *Tannerella*, *Selenomonas*, *Kingella*, *Eikenella*, *Cardiobacterium*, *Oribacterium*, *Corynebacterium matruchotii/durum*, *Porphyromonas gingivalis/endodontalis*, oral *Prevotella* (*intermedia*, *nigrescens*, *oris*, *melaninogenica*, *denticola*, …), oral *Streptococcus* (*mitis*, *oralis*, *sanguinis*, *gordonii*, *parasanguinis*, …), *Veillonella atypica/dispar/rogosae*, *Filifactor*, *Desulfomicrobium orale*, *Atopobium parvulum/rimae*, oral *Eubacterium* spp., etc.

Deliberately **excluded** from the score: (a) *Fusobacterium nucleatum*, *Parvimonas micra*, *Peptostreptococcus stomatis*, *Gemella morbillorum*, *S. anginosus* — oral-origin but genuinely gut/tumour colonising in CRC, so including them would be circular; (b) *Streptococcus salivarius*, *Veillonella parvula*, *Dialister invisus* — real stool commensals.

### Result
The score is sharply bimodal with a wide gap:

- 24 samples at 7.7–39.6 % (23 CRC, 1 control = P009), each with **2–10 different obligate-oral species individually above 1 %** (e.g. P013 *Neisseria sicca* 18.4 %, P040 *S. sanguinis* 18.9 %, P043 *Rothia dentocariosa* 16.3 % + *R. aeria* 12.9 %, P046 *P. gingivalis* 15.8 % + *T. forsythia* 3.5 % + *T. denticola* 2.9 %, P031 *P. intermedia* 11.3 %).
- gap from 7.7 % down to 4.3 %
- the remaining 90 samples are at 0.04–4.3 % and essentially all have at most **one** oral species above 1 %.

A whole periodontal community at 8–40 % of a stool metagenome is not biology — stool normally carries <1–2 % obligate oral taxa. These 24 samples carry saliva/oral material and were excluded.

Not flagged, deliberately: P053 (*S. salivarius* 13.3 %), P092 (10.7 %), P052 (6.1 %), P102 (*V. parvula* 8.3 %). These are single-species blooms of organisms that do colonise the gut, with no accompanying oral community (obligate-oral score 3.5–4.3 %) — the signature of a real bloom rather than saliva carry-over. Their inclusion/exclusion does not change any of the six calls.

## Differential abundance
On the 90 retained samples (30 CRC, 60 control):
- species present (>0) in >= 10 % of the 90 samples -> **279 species tested**
- two-sided Mann–Whitney U on relative abundance
- Benjamini–Hochberg across all 279, significance at q < 0.05

**5 significant species:** *Peptostreptococcus stomatis* (q = 7.2e-09, up), *Fusobacterium nucleatum* (q = 4.3e-07, up), *Dialister pneumosintes* (q = 1.2e-06, up), *Parvimonas micra* (q = 8.0e-03, up), *Blautia wexlerae* (q = 2.5e-02, down).

Requested six species:

| species | mean CRC | mean control | p | q | call |
|---|---|---|---|---|---|
| *Fusobacterium nucleatum* | 0.031 | 0.000 | 3.1e-09 | 4.3e-07 | enriched_in_CRC |
| *Peptostreptococcus stomatis* | 0.190 | 0.0002 | 2.6e-11 | 7.2e-09 | enriched_in_CRC |
| *Parvimonas micra* | 1.083 | 0.013 | 1.1e-04 | 8.0e-03 | enriched_in_CRC |
| *Rothia dentocariosa* | 0.0014 | 0.0012 | 0.46 | 0.76 | not_significant |
| *Eubacterium rectale* | 1.79 | 3.48 | 0.017 | 0.20 | not_significant |
| *Streptococcus salivarius* | 0.377 | 0.845 | 0.025 | 0.25 | not_significant |

## Why the QC step matters
Because 23 of the 24 contaminated samples are cases, the contamination manufactures a case–control difference for every oral taxon. Analysing all 114 samples gives 30 "significant" species including three false positives among the six asked about:

- *Rothia dentocariosa*: q = 0.018, "enriched_in_CRC". Its entire signal lives in the contaminated samples (CRC mean 0.78 % with them, 0.0014 % without). Outside that set it is ~0 in both groups.
- *Streptococcus salivarius*: q = 0.011, "depleted_in_CRC" — driven by contaminated/bloom controls.
- *Eubacterium rectale*: q = 0.013, "depleted_in_CRC". Dilution by the spiked-in oral community depresses all gut commensals in the contaminated cases; on clean samples the trend persists (median 0.29 vs 2.16 %) but does not survive BH correction (q = 0.20).

The three canonical oral-origin CRC colonisers (*F. nucleatum*, *P. stomatis*, *P. micra*) are the only ones that stay significant, and they get *stronger*, not weaker, after exclusion (e.g. *P. stomatis* q 2.2e-05 -> 7.2e-09), which is the expected behaviour for a real signal whose noise has been removed.

## Robustness
The three "enriched" and three "not_significant" calls are identical at every exclusion threshold that removes samples with more than ~10–15 % obligate-oral content (tested exclusion sets of 21, 23, 24, 25, 28, 30 and 32 samples under several marker definitions). The calls only flip if fewer than ~10 samples are removed, i.e. if obviously contaminated samples such as P022 (34 % oral) are left in. A logistic regression on log-abundance adjusted for age and gender agrees on direction for all six and keeps *P. stomatis* significant; it is unreliable for *F. nucleatum* and *P. micra* because those are perfectly/near-perfectly separated (zero in all retained controls), so the rank-based test is reported.

## Files
- `work/explore.py`, `work/oral.py`, `work/contam.py`, `work/detect.py` — QC and contamination screening
- `work/da.py`, `work/final.py`, `work/robust.py` — differential abundance and threshold sensitivity
- `work/sample_scores.csv`, `work/scores2.csv`, `work/oral_detect.csv` — per-sample scores
- `work/answer.json` — final answer
