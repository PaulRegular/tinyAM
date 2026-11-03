# Icelandic haddock (had.27.5a): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_iceland.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_iceland.md).

## Expanded-library review

### tinyAM translation revision (2026-10-14)

The translation was revised to improve numerical performance and better represent time-varying fishing mortality.

The previous configuration used IID cohort-survival deviations and a stationary AR1 F process around estimated age-specific mean fishing mortality.

The revised configuration retains the IID survival process and replaces stationary AR1 F with correlated random walks.

#### Population and recruitment

The revised model specifies:

`N_settings = list(process = "iid", init = "free")`

This uses independent deviations from expected cohort survival, with freely estimated initial abundance for older ages.

Recruitment remains a separate random-walk process under tinyAM's default recruitment settings. Specifying `process = "iid"` for N does not imply IID recruitment.

The official SAM model also estimates abundance and recruitment as latent states, but its precise process covariance, variance sharing and initial-state conventions have not been fully reconstructed in the canonical record.

The translation therefore retains a simpler population-process specification.

#### Fishing mortality

The principal revision is:

`F_settings = list(process = "cor_rw", mu_form = NULL)`

This replaces the earlier configuration:

`F_settings = list(process = "ar1", mu_form = ~factor(age))`

The previous AR1 process described fishing mortality fluctuating around estimated age-specific means, with deviations tending to return toward those means.

The revised correlated random walk allows annual changes in log fishing mortality to accumulate over time. The annual increments are correlated across neighbouring ages through an estimated AR(1) age-correlation parameter.

This provides a flexible representation of changing fishing mortality and selectivity-at-age without requiring an estimated stationary mean F curve.

The official 2025 benchmark SAM likewise allows fishing selectivity to vary over time. The revised tinyAM structure is therefore a plausible approximation of the official assessment's temporal flexibility.

However, the precise official SAM F-process covariance, variance-sharing keys and age-state sharing have not been independently verified for this translation. The use of `cor_rw` should not be described as an exact reconstruction of those settings.

#### Natural mortality

The translation retains:

`M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)`

Natural mortality is supplied as a fixed value of 0.2 across ages and years, matching the official assumption.

No additional M process is estimated.

#### Catch observations

The translation uses:

`catch_settings = list(sd_form = ~1)`

Commercial catch numbers-at-age are retained and fitted with one estimated lognormal observation SD.

Missing catch observations are not filled.

This is a simplified observation likelihood rather than a complete reproduction of the official SAM catch-error specification.

#### Survey observations and catchability

The translation retains both the spring IS-SMB and autumn IS-SMH age-specific survey series.

The observation settings are:

`index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey)`

Catchability is estimated separately for the translated survey-age groups, while observation SDs are estimated separately for the spring and autumn surveys.

The model does not impose monotone catchability or density-dependent catchability powers.

The survey indices remain on their published native scales. Their calibrated physical units were not recovered, and differences in index scale are absorbed through the estimated catchability parameters.

The translation uses approximate sampling fractions of 0.20 and 0.80. Exact timing and complete observation-error covariance equivalence with the official model remain unresolved.

#### Biological inputs

The official annual stock weights, catch weights and maturity-at-age inputs are retained.

The source's early-year stock-weight and maturity substitutions are preserved.

No alternative biological surfaces are introduced as part of the revised model configuration.

#### Spawning-stock biomass comparison

The official assessment calculates spawning-stock biomass after accounting for mortality before spawning.

Its published definition uses fractions of 0.4 for fishing mortality and 0.3 for natural mortality.

The current tinyAM translation does not apply those fractions when calculating mature biomass.

Consequently, the comparison summary uses a common-definition mature biomass calculated from the official population numbers and the same translated stock weights and maturity inputs used by tinyAM.

This quantity supports comparisons of estimated mature stock abundance but is not identical to the official assessment's native spawning-time SSB.

The official length-based reference biomass for fish of at least 45 cm is retained as a native output. It is not directly compared with tinyAM's total age-based biomass.

#### Interpretation of the revision

The revision replaces a stationary fishing-mortality process with a correlated random walk, allowing more persistent temporal variation and age-correlated annual changes.

The N process, biological assumptions, survey selection and observation-error specifications remain unchanged.

The revised configuration improves the representation of temporally varying fishing mortality without adding stock-specific biological processes or extra observation streams.

The numerical fix and any improved agreement with official outputs should be documented using the newly generated convergence diagnostics and comparison summaries. No new numerical performance statistics are asserted in this source review.
