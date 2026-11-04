# North Sea sole source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_sole_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sole_north_sea.md).

## Expanded-library review

### tinyAM translation and fit

The translation fits 1957–2024 at ages 1–10+, because source catch and F-at-age begin at age 1. It sums landings and discards by year and age where both are reported; before the discard series begins in 2002, the model input contains landings only. It fits the two age-specific Q3 indices at timing 0.625, groups q by survey and age 1–4 versus 5+, and omits the closely related aggregate indices. The index values remain on their reported relative scales. Fixed M, annual stock weights, maturity, and start-of-year SSB timing are retained.

tinyAM estimates an IID cohort-survival process, a random walk for log recruitment, and a random-walk F process with mean ages 2–6. The source's fixed-steepness Beverton–Holt curve cannot be represented with steepness fixed in tinyAM. An initial attempt to estimate the tinyAM Beverton–Holt parameters had a non-positive-definite Hessian and undefined standard errors; switching to the simpler log-recruitment random walk resolved that numerical failure without changing observations or other model settings. This omits the source stock-recruitment curve and should be treated as an approximation.

The final fit passed `check_tam()`: optimizer code 0, maximum adjusted gradient about 0.001, positive-definite Hessian, and successful uncertainty estimation. Residual diagnostics still flag some standardized residual spreads outside 0.5–1.5 and some residuals above 3 SD. The cached dashboard is local at `results/cache/ices_sole_north_sea_2025/dashboard.html`.

Common-definition agreement over 1957–2024 is mixed. SSB has 11.7% mean absolute percent difference, −41.8% terminal difference, and 0.957 trend correlation. Recruitment is 29.1%, +5.9%, and 0.662; N-at-age is 27.2%, −41.1%, and 0.913; F-at-age is 31.7%, +106.5%, and 0.731. Total abundance has 16.2% mean absolute percent difference and 0.693 trend correlation. M is identical because it is fixed to the same source input, not independently estimated. These numbers describe this simplified fit and are not evidence of equivalence to the SS3 likelihood.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_sole_north_sea.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The recovered BTS and COAST observations are aggregate biomass indices; the current supports age-specific index observations or aggregate biomass indices paired with age compositions. No matching 2026 age-composition series was recovered. Assigning these totals to ages or deriving compositions from fitted abundance would invent observations, while fitting catch alone would omit the assessment's survey likelihood.
