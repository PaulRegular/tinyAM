# Irish Sea whiting: 2025 WKBNSCS benchmark: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_whiting_irish_sea_2025.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_whiting_irish_sea_2025.md).

## Expanded-library review

### tinyAM approximation

The translation keeps all three index identities, source timing, age-specific q-sharing groups, the available NIGFS relative precision weights, catch-at-age, fixed source M and the annual stock weight and maturity surfaces. It uses an IID N process with exponential initialization and a random-walk recruitment process, approximates the SAM age-correlated F process with tinyAM's correlated random walk, and uses one catch SD plus one index SD per survey. The SAM AR residual correlation across NIGFS ages is not represented in tinyAM. The two missing NIMIK observations remain absent. Fitting first uses tinyAM's normal initialization; only if numerical checks fail does the recipe make one documented, year/age-aligned source-start retry. Source estimates do not constrain tinyAM estimates.

The default tinyAM attempt returned optimizer code 1 (`false convergence (8)`), objective 706.296, maximum absolute gradient 0.000167 and a non-positive-definite Hessian. Its estimated F-increment age correlation was 0.9803. The one source-start retry returned optimizer code 0 (`relative convergence (4)`) at the same objective, but its maximum absolute gradient was 0.00329 and the Hessian remained non-positive-definite with non-finite fixed-effect covariance. Neither attempt passed `check_tam()`, so no comparison summary or dashboard was produced. This remains a numerical limitation of the current translation; the fit has not been rescued by changing the data or model assumptions.
