# Central Baltic herring: ICES WGBFAS 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_ices_wgbfas_harengus_bs_25_29_32.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_ices_wgbfas_harengus_bs_25_29_32.md).

## Expanded-library review

### tinyAM translation and fit

The translation fits 1974–2025 and ages 0–8+, preserving the source catch-at-age, BIAS index, weights, maturity, and the medium fixed-M scenario. It omits zero age-0 acoustic values from the lognormal observation table. BIAS timing is set to 0.8 as an explicit approximation for October. The model uses an exponential initial age structure, no separate N survival process, Beverton–Holt recruitment with AR1 residuals and fixed process SD 0.5, an AR1 F process with a logistic age mean, and age-specific BIAS q. Catch and index observation SDs are estimated as common groups. No source-based starting values were used.

The default fit converged: optimizer code 0, objective 157.39, maximum absolute gradient 1.18e-6, positive-definite Hessian, and successful standard-error estimation. `check_tam()` also reports a stock-recruit parameter interval wider than tenfold, an AR1 correlation interval wider than 0.5, some residual groups with standardized spread outside 0.5–1.5, and some standardized residuals above 3 in absolute value. The fit is numerically stable, while recruitment shape and persistence remain weakly identified.

| Metric | Comparison | Mean absolute percent difference | Median absolute percent difference | Terminal percent difference | Trend correlation |
|---|---|---:|---:|---:|---:|
| SSB | Approximate | 4.01% | 3.26% | -15.16% | 0.991 |
| Recruitment | Approximate | 9.17% | 7.66% | 2.30% | 0.985 |
| M-at-age | Non-equivalent input check | — | — | — | — |
| F-at-age, N-at-age | Unavailable | — | — | — | — |
| Fbar | Non-equivalent definitions | — | — | — | — |

SSB compares the official aggregate output with tinyAM SSB over ages 0–8; the official output is the three-scenario ensemble, while tinyAM uses the medium-M component. Recruitment uses age 0 in the same calendar year, but its 2025 official value is a stock-recruit prediction. The M surface matches because it is supplied from the same medium-M scenario; it is not evidence that tinyAM independently estimated M. The official Fbar is an unweighted age-3–6 mean, while tinyAM reports an abundance-weighted mean, and no source F-at-age or N-at-age surface was recovered. The comparison does not establish likelihood equivalence.

The dashboard is generated from the Official output reference and the converged tinyAM fit. Missing official age-specific N/F/M values remain blank in those panels.
