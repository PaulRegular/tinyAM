# Northern Shelf (North Sea) haddock: 2026 source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_north_sea.md).

## Expanded-library review

### Current tinyAM translation

The recipe retains 1972–2026 and ages 0–8+, RW recruitment, IID survival
residuals, and integrated older initial abundance. `N_sd_key` in population
covariates shares survival SD over destination ages 1–7, with a separate 8+
SD. First-year recruitment remains an estimated fixed effect; SAM initial-state
integration is not claimed to be identical.

F follows `cor_rw` with one shared SD. Supplied annual M, stock weights,
maturity and zero pre-spawning mortality fractions are retained. The shared
age-correlation structure approximates SAM; first-state and integration
conventions may still differ.

Catch uses one lognormal SD. Both delta-GAM survey series retain native scales,
ages, timing and relative precision factors, with survey-by-age q and one
estimated SD level per survey. Missing observations are excluded from the
likelihood. Source density-dependent q powers at Q1 age 1 and Q3+Q4 ages 0–1
are omitted, and exact observation-variance sharing differs. Ordinary q does
not reproduce a power relationship between abundance and an index.

## Initial-state approximation

The translation uses `init = "latent"`: older first-year log-abundances are
integrated without an initial density. First-year recruitment remains an
estimated fixed effect, so this is not complete SAM initial-state equivalence.
