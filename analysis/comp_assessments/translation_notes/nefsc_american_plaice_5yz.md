# Gulf of Maine/Georges Bank American plaice (5YZ): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/nefsc_american_plaice_5yz.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/nefsc_american_plaice_5yz.md).

## Expanded-library review

### tinyAM translation

The translation fits the full Run29F4 years and ages, retains source stock,
maturity, catch, and survey weight inputs, and uses M = 0.3. It uses IID
recruitment, IID cohort-survival deviations with free initial older ages, and
an AR1 F process around a logistic age mean. q is estimated separately for
each survey and age. Source aggregate index log SDs are applied to every
age-specific index reconstructed for the matching year and survey.

These choices approximate the source model rather than reproduce it. In
particular, tinyAM fits age-specific lognormal catch and survey observations
instead of the source aggregate and composition likelihoods; it does not
propagate age-composition uncertainty, reproduce WHAM's joint random-effects
structure or selectivity constraints, or represent the source's separate age-1
and ages-2+ N process SDs. The later 2024 and 2026 summaries are
not used as input or output references for this fit.

The package default start did not pass checks because its fixed-effect Hessian
was not positive definite (optimizer code 0, adjusted maximum gradient
0.00197). The single permitted retry used source N-at-age and F-at-age aligned
to the model grid. It passed `check_tam()` (optimizer code 0, adjusted
maximum gradient 0.00353, positive-definite Hessian, and successful SD report).
The diagnostics still flag 48 pairs of highly correlated fixed parameters,
some residual groups with unusual spread, residuals above three standard
deviations, and adjacent-year residual correlation. The fit is numerically
usable but these findings warrant caution.

The dashboard's common-definition comparisons over 1980–2019 are: SSB mean
absolute percent difference 19.4%, terminal difference −36.6%, correlation
0.954; total abundance 34.8%, −46.8%, 0.743; N-at-age 38.3%, −40.9%, 0.930;
recruitment 46.1%, −60.2%, 0.852; F-at-age 48.1%, −22.5%, 0.841; and
N-weighted F over ages 6–9 33.4%, +7.9%, 0.956. These statistics describe the
candidate-run comparison and do not establish that Run29F4 was the production
advice model. Fixed M is identical by construction and is only an input check.
