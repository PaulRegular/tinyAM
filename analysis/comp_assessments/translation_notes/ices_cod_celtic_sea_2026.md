# Celtic Sea cod, 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_cod_celtic_sea_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_celtic_sea_2026.md).

## Expanded-library review

### tinyAM fit diagnostics

The package-default initialization was tried first, without `start_par`. Its
optimizer stopped with an NA/NaN gradient. One retry kept the same observations,
formulas and model settings, using the official SAM N-at-age and F-at-age
surfaces to initialize tinyAM's N, recruitment and F states. The retry returned
optimizer code 1 (`false convergence (8)`), objective 308.9221 and maximum
absolute gradient 0.00298. Standard-error estimation succeeded and the Hessian
was positive definite, but the fit does not pass tinyAM's convergence check.
The attempt records are cached locally in
`results/cache/ices_cod_celtic_sea_2026/initialization_attempts.rds`.

An exploratory fit with independent F random-walk increments and package
defaults also failed with an NA/NaN gradient. That process change was not
retained. Since neither attempt passed the convergence check, no comparison
summary or dashboard was produced.

The final SAM output in 2025 is SSB 98.94 tonnes, total stock biomass 251.39
tonnes, age-1 recruitment 183.20 thousand fish and Fbar ages 2-5 of 1.4348
per year. These values agree with the rounded values in WGCSE Table 7.9.
No tinyAM agreement summary is available until a fit passes its numerical
convergence checks; matching definitions alone do not imply matching fitted
models.
