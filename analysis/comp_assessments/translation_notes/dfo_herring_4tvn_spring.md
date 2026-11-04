# Southern Gulf spring-spawning herring: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_herring_4tvn_spring.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_herring_4tvn_spring.md).

## Expanded-library review

### Current tinyAM translation

The recipe models 1978–2023 and ages 2–11+, pooling the two commercial gears.
RW age-2 recruitment and IID older-cohort survival departures with exponential
initialization approximate the source's deterministic survival. F follows
age-correlated log random walks; Fbar uses ages 6–8. Period-specific logistic
fishery selectivity is not imposed.

M follows stationary AR1 processes for ages 2–6 and 7–11+ around a supplied
centre of 0.2, beginning in 1978. This differs from the source M random walks,
fixed increment SD 0.075 and initial-M priors. Source M estimates are not fixed
inputs or comparison observations.

CPUE and acoustic surveys have separate baselines and non-decreasing age-q
curves. Only CPUE has an annual log-q random walk: a numeric indicator is one
for CPUE and zero for acoustic observations. The first temporal effect is zero;
the survey intercept identifies its absolute baseline. The increment SD is
estimated, rather than fixed at the source value 0.1. Acoustic q is constant
through time. Source-derived population/F estimates are not optimizer starts.

Catch and indices use independent age-specific lognormal errors, with one catch
SD and separate survey SDs. This replaces the source total-plus-composition
likelihoods and acoustic biomass weighting. Model catch and count-based indices
use the converter's thousand-fish scale; missing or zero observations remain
excluded, with complete predictions available. CPUE timing 0.25 and acoustic
timing 0.75 are explicit approximations. The acoustic survey's documented unit
and the distinction between January 1 mature biomass and official April 1 SSB
are retained in the definitions above. Population processes and q can compensate
for one another; convergence does not establish source-model equivalence.

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

## Record summary

The canonical record represents the official 2024 Southern Gulf spring-spawning herring assessment, with detailed inputs and outputs through 2023.

The latest tinyAM translation retains the published age-specific catch, CPUE and acoustic observations and applies a carefully checked approximation to beginning-year stock weights.

It uses random-walk recruitment, IID cohort-survival variation, correlated F random walks, stationary two-block AR1 natural mortality and a CPUE-specific random walk in catchability.

This model is more consistent with the documented temporal CPUE catchability process than the earlier three-period fixed-effect approximation, although its estimated process SD differs from the documented source setting.

The translation is retained as a practical approximation of a more complex SCA assessment. Its principal unresolved differences involve the original observation likelihood, selectivity constraints and mortality-process priors.

Source completeness remains partial, and the latest numerical performance should be recorded separately from the historical refinement results.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/dfo_herring_4tvn_spring.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The translation combines the two fishing gears into one total catch-at-age stream.

The translation uses age-specific lognormal survey observations rather than the source SCA's age-composition and aggregate-biomass likelihoods.

The translation constructs a fit-only approximation by combining positive gear weights using catch-number-weighted arithmetic means. When both gear catches are zero, it uses the arithmetic mean of available positive gear weights.

Missing weights are interpolated within age.

For an ordinary interior cell, beginning-year weight at age a in year t is calculated as the geometric mean of the combined gear weight at age a−1 in year t−1 and age a in year t.

Age 2 and the first model year use available same-year values because preceding cells are unavailable.

This procedure is an inference, not a fully verified source calculation.

Its median absolute percentage difference from weights implied by the published biomass and abundance surfaces is 0.073%, compared with 2.63% for the earlier equal geometric gear mean.

Source biomass divided by abundance is used to evaluate the weight approximation, not as a fitted model input.

The common-definition SSB comparison therefore uses January 1 mature biomass. It is not the official assessment's April 1 SSB, which would require pre-spawning survival calculated from the unavailable fitted M surface.
