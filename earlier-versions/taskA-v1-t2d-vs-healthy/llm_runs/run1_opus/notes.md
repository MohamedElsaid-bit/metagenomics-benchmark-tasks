# T2D vs healthy species differences: analysis notes

Data: MetaPhlAn3 species rows (s__ level, no t__), 1831 samples. T2D = disease contains "T2D" (n=792); healthy = "healthy" (n=234).
Species with prevalence >= 10% in the T2D+healthy set were tested on log10(rel. abundance + 1e-3). Wilcoxon tests and OLS models were used, with BH-FDR correction.

Scripts: work/analysis.py (main), work/explore.py. Full-model table: work/full_model.tsv.

## Findings
- **Naive T2D vs healthy:** 89 species differ (q<0.05). E. coli is strongly increased. F. prausnitzii, R. hominis and A. muciniphila are decreased, which matches the "classic" signature.
- **Medication confounding:** controls take essentially no drugs. 51.8% of T2D participants (410/792) take metformin. The other drugs are statins 56%, aspirin 40%, PPI 19%, DPP4i, SU and insulin.
- **Within T2D, drugs mutually adjusted** (plus age, BMI, country and gender): metformin is associated with 81 species at q<0.05, PPI with 11, and statin, DPP4i, SU, insulin and aspirin with 0.
  - Metformin effects: E. coli increased (b=+1.00 log10, q=1e-12), I. bartlettii decreased (q=6e-17), F. prausnitzii decreased (q=4e-8), R. hominis decreased (q=2e-6), Akkermansia not changed (q=0.90).
- **Metformin-naive T2D vs healthy:** the number of differing species falls from 89 to 30.
  - After adjusting for age, BMI, country and gender, none of the five focal species differ: E. coli p=0.11, I. bartlettii p=0.29, F. prausnitzii p=0.67, R. hominis p=0.51, Akkermansia p=0.67.
  - The unadjusted Akkermansia depletion (q=0.03) disappears once BMI and the other covariates are included.
- **Full T2D + healthy model:** adding the drug indicators to the T2D-status model cuts the significant T2D-status species from 53 to 4. None of the focal species remain significant; E. coli has q=0.21.

## Conclusion
Metformin is the main driver of the apparent T2D signature (more E. coli, fewer butyrate producers and Intestinibacter). This agrees with the published MetaCardis analysis (Forslund et al.).
- **Effect of T2D itself:** no robust difference for any of the five species.
- **Akkermansia:** its apparent depletion is explained by BMI and other covariates, not by T2D or metformin.
