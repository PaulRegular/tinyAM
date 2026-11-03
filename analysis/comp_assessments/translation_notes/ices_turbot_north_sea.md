# North Sea turbot source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_turbot_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_turbot_north_sea.md).

## Expanded-library review

### tinyAM translation and fit

The recipe `ices_turbot_north_sea_2025_benchmark.R` translates the 1981–2023, age 1–8+ series. It keeps both surveys and their 0.705 timing, matches the q sharing and observation-SD groups where possible, uses IID cohort-survival errors and random-walk recruitment, fits correlated random-walk F with four age-based process-SD groups, and supplies fixed age-varying M. It uses the package's default initialization and no custom starting values. The differences from SAM's age covariance, tied F ages 7/8, BTS+COAST residual correlation and weight-state estimation are described in the recipe background.

The M9 output includes fitted stock weights for the complete 344-cell year-age grid, while the canonical observed input series has 242 cells. The 102 missing cells are supplied to the analysis-only observation object from `exp(logSW)` in the official M9 fit. The original observed values are retained unchanged, and the supplement is stored in `outputs.csv` as a derived source output rather than misclassified as an input. tinyAM treats the complete weight surface as fixed.

The fit converged without a starting-value retry: optimizer code 0 (`relative convergence (4)`), objective 472.75, maximum absolute gradient 0.000502, and positive-definite Hessian with standard errors available. Residual diagnostics still identify large standardized residuals, residual spread outside the reference range, and adjacent-year residual correlation; these are retained as model-fit caveats.

The dashboard compares Official and tinyAM outputs. Biomass and SSB use common definitions: official N is combined with the completed weight and maturity surfaces used by tinyAM, rather than compared with SAM native biomass from its estimated latent weight surface. Common-definition mean absolute percent differences are 21.1% for SSB (terminal difference +26.7%, trend correlation 0.920), 24.9% for N-at-age (terminal +25.8%, mean within-age correlation 0.890), 20.4% for F-at-age (terminal −23.4%, mean within-age correlation 0.877), and 13.4% for recruitment (terminal +11.5%, correlation 0.902). Total abundance differs by 10.6% on average (terminal +16.3%, correlation 0.951), total biomass by 11.7% (terminal +19.7%, correlation 0.979), and common-definition Fbar by 10.2% (terminal −20.2%, correlation 0.950). Fixed M is consistent with the supplied input; it is not an estimated fit comparison. Biomass-at-age is unavailable in the comparison. Derived aggregate quantities have no propagated uncertainty.
