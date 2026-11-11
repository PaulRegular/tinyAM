# Scup SAW 60 translation

The comparison uses 1984–2014 and ages 0–7+ from the 2015 S60_BASE_18 assessment. The four-fleet fit retains commercial landings, commercial discards, recreational landings and recreational discards. Recreational discards contribute at the reported 15% mortality rate. The four age-specific catch series sum to the baseline dead removals in every year and age. All fleets use the source's total-fishery weight-at-age table. The plus-group weight and June spawning midpoint remain approximations.

The trial used fleet-specific age means with an AR1 F process within each fleet and shared process parameters. It held the survey, M, recruitment and observation variance settings fixed. The source uses time-blocked selectivity, additional surveys and separate catch-total and composition likelihoods. [2015 assessment report](https://repository.library.noaa.gov/view/noaa/4975)

Both fits passed numerical checks on SAIO `cecb45ad` and all results cover 31 common years.

| Metric | Aggregate median absolute difference | Four-fleet median absolute difference | Four-fleet terminal difference | Four-fleet trend correlation |
|---|---:|---:|---:|---:|
| SSB | 64.2% | 57.4% | −50.5% | 0.983 |
| Recruitment | 41.4% | 26.2% | −49.5% | 0.836 |
| Abundance | 60.1% | 35.2% | −42.9% | 0.895 |
| Fbar | 150.3% | 53.6% | +85.0% | 0.933 |

The four-fleet fit improved the scale of all four comparisons. Correlations declined for Fbar and F-at-age, and residual diagnostics flagged spread and adjacent-year correlation. Both models use the source's fixed M of 0.20 per year. The common-definition 2014 SSB reference is 1.2% below the published value. The four-fleet fit estimates 43% fewer fish in 2014, with higher F yielding a similar total catch. Age 7+ accounts for 75% of the terminal SSB gap and has seven positive fall survey observations in 31 years. Fitted fall-survey q was 1.27–1.47 times the q implied by source states at ages 0–4 on the same index data, consistent with catchability compensating for lower abundance. The source fit uses 12 survey series, compared with one here. Recovering the remaining age-specific index inputs and likelihood settings would help determine how much they explain the scale difference.

## Curated two-survey model

The curated recipe reads the CTDEEP fall trawl index from the 2015 SAIO record. Table A43 reports fish per tow by age for 30 surveys between 1984 and 2014. The survey did not run in 2010. The reported age rows reconcile to the all-age totals, which are kept in SAIO for source checks. The recipe pools ages 7–10+ into the modeled 7+ group, uses 0.75 for the September–October survey midpoint, and estimates separate age-specific q and one log-SD for each survey. CTDEEP supplies 13 positive observations at age 7+, including 2014. [Official assessment, Table A43](https://archive.asmfc.org/uploads/file/55d238fb2015ScupBenchmarkStockAssessment.pdf)

The two-survey trial added the same CTDEEP observations to SAIO `cecb45ad` before the canonical import. It used the four-fleet trial's catch, biology, M and F settings and converged from defaults. Across 31 years, the median absolute SSB difference fell from 57.4% to 11.6%, with correlation 0.999 and a terminal difference of +11.6%. The 2014 age-7+ SSB moved from 41,969 to 115,422 t, against a common-definition reference of 110,301 t. Median absolute differences for abundance, recruitment and Fbar were 7.5%, 10.5% and 10.9%. Residual spread and adjacent-year correlation remain flagged for the fall surveys. The comparison still uses two of the official survey series and an approximation to the native likelihood.
