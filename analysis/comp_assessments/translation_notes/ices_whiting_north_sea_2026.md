# North Sea whiting source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_whiting_north_sea_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_whiting_north_sea_2026.md).

## Expanded-library review

### tinyAM translation revision (2026-10-14)

The translation was revised to better approximate the official SAM process structure.

The previous configuration used random-walk cohort-survival deviations and independent fishing-mortality random walks. The revised configuration uses:

- `N_settings = list(process = "iid", init = "exp")`
- `F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:5)`
- `M_settings = list(process = "off", mu_supplied = ~M_assumption)`
- `catch_settings = list(sd_form = ~1)`
- `index_settings` with survey-by-age catchability, survey-specific baseline SDs, supplied relative SD factors, and no filling of missing observations.

The recruitment process remains a random walk under the default tinyAM recruitment settings, separate from the `N_settings` process. The `iid` N process describes independent deviations from expected cohort survival rather than a random walk in survival deviations.

The `latent` setting integrates older initial abundances without an initial
density. First-year recruitment remains an estimated fixed effect.

The F process now uses `cor_rw`, allowing AR(1) correlation across ages in annual log-F increments. This represents an important component of the official SAM fishing-mortality process that was absent from the earlier independent-RW translation.

The revised F process uses one shared increment SD. It is therefore a closer approximation to the source process but is not necessarily identical to the official model's complete variance specification.

Catch-at-age observation errors remain independent and use one shared lognormal SD, rather than the source's age-specific observation variance groups.

Survey observations retain their native timings and survey-by-age catchability groups. The revised model estimates separate baseline observation SDs for Q1 and Q3, modified by the relative SD factors recovered from the official SAM input. Unlike SAM, tinyAM does not model AR(1) correlation among survey residuals at adjacent ages.

The numerical M observations through 2022 are supplied directly. The final available M values are held fixed for 2023–2026 to provide the complete mortality surface required by the translation. Unlike SAM, tinyAM does not estimate a mortality GMRF.

An earlier attempt to estimate M in tinyAM resulted in non-finite optimizer gradients and was not retained. The fixed-M extension is a translation-only assumption and does not modify the official assessment database.

Official N and F estimates initialize the tinyAM population and fishing-mortality states. Initial survey catchability coefficients are derived from the observed indices and official N, F, and M estimates at the survey sampling times. These are starting values only, not constraints on the fitted states.

The revised process settings substantially change the translated model relative to the earlier baseline. Agreement and convergence statistics should be interpreted for the new configuration rather than carried forward from the previous translation.

## Initial-state approximation

The translation uses `init = "latent"`: older first-year log-abundances are
integrated without an initial density. First-year recruitment remains an
estimated fixed effect, so this is not complete SAM initial-state equivalence.

## Initial-state approximation

The revised translation remains a proof-of-concept approximation of the official model. Numerical convergence and improved agreement do not establish equivalence of the underlying likelihoods or assessment uncertainty.
