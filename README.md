# Analysis code 

Analysis code for *"Extreme Temperature Events and Hospitalizations for 59 Autoimmune
Diseases in Korea: A Nationwide Outcome-Wide, Time-Stratified, Case-Crossover Study."*

Repository: https://github.com/dajeongham/Extreme-Temperature-Events-and-Hospitalizations-for-59-Autoimmune-Diseases

This repository contains the **analysis code** (exposure definition, statistical model,
multiplicity correction, subgroup analysis). The National Health Insurance Service (NHIS) claims data underlying this study cannot be
shared, owing to data-sharing restrictions imposed by the NHIS of Korea. 

## Software
- SAS 9.4 - exposure definition, conditional logistic regression, subgroup analyses
- R 4.4.1 - multiplicity correction (package: dplyr)

## Files
| File | Description |
|------|-------------|
| `01_exposure_definition.sas` | District-specific temperature percentiles (warm 95/97.5/99th, cold 1/2.5/5th) and single-/2-/3-consecutive-day heat-wave and cold-spell indicators. |
| `02_conditional_logistic_regression.sas` | Time-stratified case-crossover model: conditional logistic regression (PROC LOGISTIC, STRATA), humidity adjusted with a cubic B-spline (df=3); one model per disease x exposure definition. |
| `03_subgroup_analysis.sas` | Main model re-fitted within sex / age / income / urbanicity subgroups; effect modification tested with Wald tests (Altman & Bland 2003). |
| `04_fdr_bh_by.R` | Benjamini-Hochberg and Benjamini-Yekutieli correction within each exposure definition across the 59 diseases. |
| `05_fdr_nested_sensitivity.R` | FDR re-applied after excluding nested outcomes (Crohn's / ulcerative colitis within IBD). |

