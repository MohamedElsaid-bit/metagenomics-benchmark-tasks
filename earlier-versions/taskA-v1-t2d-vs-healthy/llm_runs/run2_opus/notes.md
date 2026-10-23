# T2D vs healthy: species-level drivers

Groups: T2D = disease contains "T2D" (n=792), healthy = disease == "healthy" (n=234). Species = MetaPhlAn3 `s__` rows (no strain rows), prevalence >=10% in these samples, log10(rel. abundance + 1e-4). OLS per species with BH FDR. Scripts: `explore.py`, `analyze.py`.

- Medication use differs sharply between the groups: 86% of controls take nothing; among T2D, 51.8% take metformin, 55.6% statins, 33% DPP-4 inhibitors, 19% PPIs.
- T2D vs healthy, unadjusted: 84 species at q<0.05. E. coli up, F. prausnitzii and R. hominis down, which matches the published signature. Adjusting for age, BMI, sex and country still leaves 56 species.
- Adding the medication indicators (metformin, statin, DPP-4 inhibitor, sulfonylurea, insulin, PPI, aspirin) to the model drops this to 4 species. None of the five target species stays significant (E. coli nominal p=0.02, q=0.21).
- Metformin-naive T2D vs healthy: 8 species at q<0.05, and none of the five targets. Fully untreated T2D (n=64) vs healthy: 0 species.
- Within T2D, testing each drug while adjusting for the other drugs and the covariates: metformin 75 species at q<0.05; PPI 11; statin, DPP-4 inhibitor, sulfonylurea, insulin and aspirin 0 each.
- Metformin within T2D: E. coli up (q~1e-14), Intestinibacter bartlettii down (q~1e-16), F. prausnitzii down (q~1e-7), R. hominis down (q~1e-7), Akkermansia no difference (q=0.90). The Mann-Whitney tests agree.
- Akkermansia: the raw medians are lower in T2D, but there is no difference after adjusting for BMI and the other covariates.

Conclusion: the "T2D signature" in this cohort comes mostly from metformin treatment, not from the disease. Once medication is accounted for, none of the five species differs robustly between T2D and controls.
