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
