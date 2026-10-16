# Southern Gulf spring-spawning herring

Translation revision: 2026-10-14.

## Assessment record

The canonical database represents the official 2024 assessment to 2023, the most recent assessment with detailed numerical tables recovered for this stock.

DFO Research Document 2024/058 is primarily a data and supporting document. It does not provide a complete specification of the fitted statistical catch-at-age (SCA) model.

The official 2026 assessment to 2025 is recorded separately as summary-only because the detailed model documentation and complete numerical input and output tables have not been recovered. No 2026-run data or estimates are mixed into the 2024 record.

DFO Research Document 2022/068 provides the most recent detailed SCA methods description used in this review. Its equations and parameter settings are retained as methodological context, but are not assumed to be independently verified for every component of the 2024 assessment.

Source material is cached under `source_cache/dfo_sgsl_herring_spring_2023/`, including the 2024 supporting document, 2022 methods report, 2026 Science Advisory Report, acoustic biomass CSV and dictionary, historical acoustic unit documentation, and file manifest.

## Inputs and assumptions recovered

### Catch and population structure

The 2024 document provides fixed- and mobile-gear spring catch-at-age and weight-at-age, spring fixed-gear CPUE by age, acoustic age-index values, and January 1 abundance, biomass-at-age and fishing-mortality-at-age estimates.

The assessment covers 1978–2023 and ages 2–11+, with age 11 as the plus group.

The 2023 mobile-gear spring catch-at-age is reconstructed as zero because the annual landings table reports no mobile-gear spring catch that year. This is identified as a reconstruction rather than an independently reported catch-at-age series.

The translation combines the two fishing gears into one total catch-at-age stream. The source SCA observation model uses total catch and age-composition likelihoods rather than the independent age-specific lognormal likelihood used by tinyAM.

### Survey observations and acoustic units

The two survey series are:

- Spring fixed-gear CPUE, available for 1990–2021.
- The spring-spawner 4Tmno acoustic survey, with age-index observations for ages 2–10 in 1994–2023.

The 2024 supporting report describes acoustic trawl samples weighted by acoustic density to allocate the survey signal by spawning component and age. Its Table 15 is an abundance-at-age index, not raw trawl sample counts.

Table 15 does not identify a numerical scale multiplier. DFO Research Document 2016/060, Table 16, explicitly labels the corresponding historical acoustic series in thousands of fish.

For example, the 1994 age-4 observation is 100,087 in the earlier table and 100,062 in the 2024 report; the 2006 age-4 observation is 24,601 in both. These differences are consistent with historical revisions rather than a change in units.

The canonical database therefore retains the published 2024 acoustic index values in thousand fish. The shared converter multiplies them by 1,000 before fitting tinyAM.

A separate biomass-scale check supports this interpretation. The 2023 acoustic age-index values sum to 55,619 thousand fish, while a later official spring biomass estimate of 7,358 t implies an average weight of approximately 0.132 kg per fish. Treating the table as individual fish would imply approximately 132 kg per fish.

The later biomass estimate is used only for validation and is not introduced as an observation in the 2024 model translation.

### Acoustic biomass provenance

The Open Government acoustic biomass CSV and its dictionary report spring biomass in tonnes over 1994–2025.

These values are retained separately from the fitted 2024 assessment. In particular, 2024–2025 observations are not inserted into the 2024 translation.

The latest CSV is not confirmed to reproduce every historical biomass input used by the official 2024 SCA model.

Examples include:

- 2022: 27,209 t in the latest CSV versus 27,268.7 t in the 2024 supporting report.
- 2023: 19,363 t in the latest CSV versus 19,363.8 t in the report, consistent after rounding.
- 2021: 23,147.68 t in the latest CSV versus 37,953.1 t in the 2022 methods report.

The reasons for these revisions and the exact historical aggregate-biomass likelihood inputs remain unresolved.

## Documented SCA model structure

DFO Research Document 2022/068 describes the SCA model used as methodological context for the 2024 assessment.

### Recruitment and abundance

Recruitment enters at age 2 and varies annually.

Older cohorts are propagated through survival, with initial-cohort deviations handled separately from the ordinary annual recruitment process.

The documented SCA does not contain the same independently estimated older-cohort survival deviations used by tinyAM.

### Fishing mortality

Fishing mortality is represented using annual fishing intensity and logistic fishery selectivity.

The documented selectivity periods are:

- 1978–1989.
- 1990–2004.
- 2005 onward.

A schematic representation is:

`F[a,t] = F[t] * S[a,period(t)]`

where the logistic selectivity curve varies among periods.

This imposes a structured age pattern on fishing mortality rather than independently estimating a full time series of F for every age.

The 2024 supporting report does not fully restate these equations or confirm every historical parameter grouping.

### Natural mortality

Natural mortality varies between years in two age groups:

- Ages 2–6.
- Ages 7–11+.

The documented model uses random walks on log M, with fixed annual increment SD of 0.075 and initial-M prior means of 0.2.

The initial-state prior is material: the random walks alone do not anchor the initial mortality level.

The numerical fitted M-at-age surface was not recovered from the 2024 assessment tables.

### CPUE catchability

The 2022 methods describe a single annual fully recruited CPUE catchability coefficient, varying through a random walk on the log scale.

The documented process is:

`log(q[t]) = log(q[t-1]) + epsilon[t]`

with annual increments distributed as:

`epsilon[t] ~ N(0, 0.1^2)`

The 1990 catchability establishes the initial level.

The documented CPUE observation structure also uses age-dependent selectivity, rather than unrelated age-specific annual catchability trajectories.

The temporal catchability process and the age-selectivity pattern are distinct model components.

The fixed increment SD of 0.1 is the last documented value. Its continued use in the official 2024 run has not been independently confirmed.

### Observation likelihoods

The documented SCA uses age compositions and aggregate biomass indices for spring CPUE and the acoustic survey.

The 2022 methods describe a multivariate-logistic age-composition likelihood. Small age proportions are grouped with adjacent ages until they exceed 0.01.

The acoustic biomass likelihood is lognormal, uses ages 4–8, and carries likelihood weight 3. The documented CPUE biomass component uses ages 4–10 and weight 1.

The 2024 supporting document does not confirm that these age groupings and weights remained unchanged in every respect.

These features are not reproduced by the current tinyAM translation, which uses age-specific lognormal survey observations.

## Biological inputs and derived comparison outputs

### Weights and maturity

Maturity is knife-edge between ages 3 and 4: ages 2–3 are immature, and ages 4–11+ are mature.

The documented beginning-year stock-weight construction combines fixed- and mobile-gear weights, then uses a geometric mean across adjacent age-year cells.

The complete fitted stock-weight surface and precise gear-combination rule were not published.

The tinyAM translation constructs a fit-only approximation by combining positive gear weights using catch-number-weighted arithmetic means. When both gear catches are zero, it uses the arithmetic mean of available positive gear weights.

Missing weights are interpolated within age.

For an ordinary interior cell, beginning-year weight at age a in year t is calculated as the geometric mean of the combined gear weight at age a−1 in year t−1 and age a in year t.

Age 2 and the first model year use available same-year values because preceding cells are unavailable.

This procedure is an inference, not a fully verified source calculation.

Its median absolute percentage difference from weights implied by the published biomass and abundance surfaces is 0.073%, compared with 2.63% for the earlier equal geometric gear mean.

Source biomass divided by abundance is used to evaluate the weight approximation, not as a fitted model input.

### Population comparison definitions

The canonical database includes derived total January 1 abundance and biomass, obtained by summing the reported ages 2–11+.

It also includes January 1 mature biomass, calculated from the reported biomass-at-age and knife-edge maturity schedule.

These derived values use the source's published biomass surface, not tinyAM's approximate stock weights. They are identified as derived source outputs, without additional uncertainty estimates.

The source's published biomass-at-age represents January 1 biomass.

The common-definition SSB comparison therefore uses January 1 mature biomass. It is not the official assessment's April 1 SSB, which would require pre-spawning survival calculated from the unavailable fitted M surface.

## tinyAM translation revision (2026-10-14)

The latest translation retains the documented observations and fit-only biological transformations but revises the F and CPUE catchability processes.

The principal changes relative to the earlier baseline are:

- Correlated random-walk fishing mortality replaces stationary AR1 F around age-blocked means.
- CPUE catchability uses an annual random walk rather than three fixed temporal levels.
- A separate survey baseline is estimated so the first CPUE year has an identifiable absolute catchability.
- The CPUE random-walk increment SD is estimated rather than fixed at the documented value of 0.1.
- The acoustic survey is excluded from the CPUE random walk.

No source-derived population or fishing-mortality estimates are used as optimizer starting values in the latest script.

### Model period and population ages

The translated model fits 1978–2023 and ages 2–11+, retaining age 11 as the plus group.

The two commercial fishing gears are combined into one catch stream. Both CPUE and acoustic survey observations are retained on their respective published age ranges.

### Recruitment and abundance

The translation specifies:

`N_settings = list(process = "iid", init = "exp")`

Age-2 recruitment follows tinyAM's default random-walk recruitment process.

Initial older-age abundance is obtained through the exponential survivorship initializer.

Older-cohort survival deviations are estimated as independent process errors.

This is a more flexible population model than the official SCA's documented deterministic cohort-survival structure. The additional survival process is retained as part of the numerically workable approximation.

### Fishing mortality

The current translation specifies:

`F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 6:8)`

Fishing mortality follows age-correlated random walks on the log scale.

For age a and year t, annual changes are:

`log(F[a,t]) - log(F[a,t-1]) = eta[a,t]`

with correlation across ages described by an AR1 covariance for the annual increments.

The process permits persistent temporal changes in age-specific fishing mortality while encouraging related ages to change together.

Unlike the official SCA, it does not impose period-specific logistic fishery selectivity or a common annual fully selected fishing-mortality coefficient.

Period-specific logistic F mean curves with AR1 deviations were explored but did not produce a satisfactory converged fit. They are not retained in the latest translation.

The script still constructs `fishery_period` using the documented 1978, 1990 and 2005 boundaries, but the final `F_settings` does not reference that covariate. It has no effect on the fitted F process.

The previous age-blocked mean-F covariate is likewise not used in the final `F_settings`.

Fbar is calculated over ages 6–8.

### Natural mortality

The current translation specifies:

`M_settings = list(process = "ar1", mu_form = NULL, mu_supplied = ~M_process_center, age_breaks = c(2, 7, 11), first_dev_year = 1978L)`

The supplied mortality centre is 0.2.

Natural mortality is represented by two age-blocked, stationary AR1 processes on the log scale, corresponding to ages 2–6 and 7–11+.

The processes fluctuate around the supplied mortality centre rather than following unrestricted random walks.

The temporal and age-block correlation parameters and process variation are estimated.

This differs from the documented SCA in three important ways:

- The source uses random walks rather than stationary AR1 processes.
- The source fixes the M increment SD at 0.075.
- The source has explicit initial-M priors that are not reproduced by the translation.

The AR1 process is retained because it provides a stationary mortality reference and has been more numerically manageable in this translation.

Source M estimates are not used as fixed inputs or comparison observations.

### CPUE catchability

The current translation uses a separate baseline and age-selectivity curve for each survey, with an annual random walk restricted to CPUE.

Its catchability formula is:

`q_form = ~0 + survey + mono(age, by = survey) + rw(cpue_year, by = cpue_active)`

The link is logarithmic.

The ordinary survey terms estimate separate absolute catchability levels for the CPUE and acoustic series.

The `mono(age, by = survey)` term estimates a non-decreasing age-catchability curve for each survey. These flexible age curves can plateau but do not enforce the logistic shape described in the SCA methods.

The CPUE random-walk term allows a single annual catchability multiplier shared across CPUE ages.

For CPUE, the fitted log-catchability structure is:

`log(q[a,t]) = alpha_CPUE + g_CPUE(a) + u[t]`

where the first CPUE random-walk deviation is zero and subsequent deviations follow:

`u[t] - u[t-1] ~ N(0, sigma_q^2)`

The random-walk increment SD, `sigma_q`, is estimated rather than fixed.

A fixed SD of 0.1, as documented in the 2022 methods, was tested but did not produce a converged translation fit. Estimating the SD was more numerically successful and is retained.

This is a deliberate simplification and departure from the last documented source penalty, not confirmation that the official SCA estimates the same parameter.

The first CPUE year defines the random-walk baseline. Its absolute catchability is estimated through the survey intercept and age curve; it does not depend on a separately estimated first-year random-walk deviation.

The implementation uses a numeric CPUE indicator. The indicator equals one for CPUE observations and zero for acoustic observations, so the annual random walk does not affect acoustic catchability.

The acoustic survey retains its separate baseline and non-decreasing age curve, constant through time.

The previously fitted three-period CPUE catchability terms are not retained.

### Survey observation likelihoods

The translation uses:

`index_settings = list(q_form = ~0 + survey + mono(age, by = survey) + rw(cpue_year, by = cpue_active), q_link = "log", sd_form = ~0 + survey)`

Separate observation SD parameters are estimated for CPUE and acoustic observations.

The acoustic index is fitted using its documented thousand-fish scale converted to fish through the canonical data converter.

The CPUE and acoustic age-index observations are fitted using independent age-specific lognormal errors.

This differs substantially from the documented SCA age-composition and aggregate-biomass likelihoods. In particular, the official acoustic biomass weighting is not reproduced.

CPUE timing is approximated at 0.25 of the year and the acoustic survey at 0.75. The latter reflects the late September to early October survey window.

Neither survey timing is claimed to reproduce the exact native SCA timing convention.

In the supplied latest script, `fill_missing` is omitted from the catch and index settings. Under the current tinyAM defaults, this activates missing-observation filling for both components. The earlier translation explicitly used `fill_missing = FALSE`.

This implementation choice should be made explicit after verifying that the default was intended. Filled observations are latent model quantities, not additional measured data.

### Catch observations

The translated catch model retains summed fixed- and mobile-gear catch numbers-at-age, converted from thousand fish to individual fish.

One estimated catch observation SD is used:

`catch_settings = list(sd_form = ~1)`

The source's combined total-catch and age-composition likelihood is not reproduced.

### Model interpretation

The revised model is a relatively flexible age-structured state-space approximation to the official SCA.

Its key structural features are random-walk recruitment, independent older-cohort survival deviations, age-correlated random-walk F, stationary two-block AR1 M, and an estimated CPUE catchability random walk.

These processes can compensate for one another, especially because the translation does not reproduce the original aggregate survey biomass and age-composition likelihoods.

Consequently, numerical convergence and agreement with official population estimates should be evaluated separately from exact structural equivalence.

The latest model is retained as a workable approximation rather than characterized as a complete reproduction of the official assessment.

## Earlier refinement trials (historical results)

The following experiments were documented on 2026-10-10 using earlier translation configurations. Their results provide development history and should not be interpreted as current-model performance.

| Trial | Converged / positive Hessian | Median absolute % difference: total N | January 1 mature biomass |
|---|---|---:|---:|
| Earlier revised-input baseline | Yes / yes | 42.2 | 17.0 |
| Earlier log-q and unresolved-scale configuration | Yes / yes | 42.2 | 17.0 |
| Restrict CPUE to ages 4–10 and acoustic to 4–8 | Yes / yes | 32.8 | 22.6 |
| M random walk | Yes / yes | 49.7 | 27.5 |
| AR1 rather than IID older-age N | Yes / yes | 53.6 | 27.1 |
| Three-degree-of-freedom spline for CPUE year | Yes / yes | 50.5 | 31.8 |
| No older-age N process | No / no | — | — |
| Free initial abundance | No / no | — | — |
| Random initial abundance | No / no | — | — |
| Earlier final fit with revised gear-weight construction | Yes / yes | 42.2 | 19.0 |

The former three-period CPUE baseline also underwent an annual CPUE catchability RW sensitivity test. That earlier candidate converged and improved some trajectory correlations, but mean January 1 mature-biomass error increased from 24.9% to 67.1%, and fishing-mortality error increased from 62.8% to 192.3%.

Those results were obtained before the current survey-baseline parameterization and other model revisions. They should not be used as estimates of the present configuration's performance.

The earlier final model had optimizer code 0, objective 1173.766 and a positive-definite reported Hessian. Its maximum projected gradient after the acoustic unit correction was 0.0000538.

The earlier model's mean absolute percentage difference in January 1 mature biomass was 24.9%, with a terminal-year difference of −33.8%.

Its terminal abundance and recruitment remained substantially below the official source estimates.

These values are retained as historical benchmarks only. The latest CPUE-RW and correlated-F configuration requires its own numerical convergence and comparison results before updated performance claims can be recorded.

The earlier M and initial-state sensitivity tests do not establish that alternative biological assumptions are scientifically inappropriate; they describe the outcomes of the tested translation configurations and starting conditions.

## Remaining translation limitations

The principal limitations of the current translation are:

- The documented age-composition and aggregate-biomass likelihoods are replaced by direct age-specific lognormal observations.
- The acoustic biomass likelihood weight is not reproduced.
- CPUE and acoustic selectivity are modelled with monotone age curves rather than the documented logistic forms.
- The CPUE catchability random-walk SD is estimated rather than fixed at 0.1.
- The F process is an age-correlated random walk without the source's period-specific logistic fishery selectivity.
- The M process is stationary AR1 rather than the documented random walk with fixed increment SD and initial-state priors.
- IID older-cohort survival deviations are added to the documented SCA survival structure.
- The stock-weight surface remains a closely checked but inferred approximation.
- The comparison uses January 1 mature biomass rather than native April 1 SSB.
- The official age-specific M estimates and uncertainty surfaces remain unavailable.
- Exact observation-model settings for the official 2024 assessment are not fully recoverable from the supporting document.

These differences explain why the model should be evaluated as a translation of the major population dynamics and available observations, not as an exact reproduction of the official statistical likelihood.

## Sources

- [DFO Research Document 2024/058](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41256384.pdf)
- [DFO Research Document 2022/068](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41091589.pdf)
- [DFO Research Document 2016/060, Table 16: acoustic units](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/365860.pdf)
- [DFO Science Advisory Report 2026/028](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-6/Fs70-6-2026-028-eng.pdf)

The relevant source documents, cached numerical files and manifest are retained under `source_cache/dfo_sgsl_herring_spring_2023/`.

The stock does not yet have a verified Charbonneau identifier in the local crosswalk, so no identifier is assigned.

## Record summary

The canonical record represents the official 2024 Southern Gulf spring-spawning herring assessment, with detailed inputs and outputs through 2023.

The latest tinyAM translation retains the published age-specific catch, CPUE and acoustic observations and applies a carefully checked approximation to beginning-year stock weights.

It uses random-walk recruitment, IID cohort-survival variation, correlated F random walks, stationary two-block AR1 natural mortality and a CPUE-specific random walk in catchability.

This model is more consistent with the documented temporal CPUE catchability process than the earlier three-period fixed-effect approximation, although its estimated process SD differs from the documented source setting.

The translation is retained as a practical approximation of a more complex SCA assessment. Its principal unresolved differences involve the original observation likelihood, selectivity constraints and mortality-process priors.

Source completeness remains partial, and the latest numerical performance should be recorded separately from the historical refinement results.
