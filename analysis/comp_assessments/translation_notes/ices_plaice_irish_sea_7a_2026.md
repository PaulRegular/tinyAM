# Irish Sea plaice (`ple.27.7a`): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_plaice_irish_sea_7a_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_plaice_irish_sea_7a_2026.md).

## Expanded-library review

### tinyAM translation

The translation retains 1981–2025 and ages 1–8+, fits the UK BTS series, and uses the official fixed M, weights, and maturity. It uses tinyAM random walks for N and recruitment, correlated-RW F, a shared UK BTS observation SD, and q by age. This approximates SAM's process covariance, boundary integration, and observation model; the exact SAM likelihood and age-state sharing are not reproduced. The aggregate NIGFS indices are omitted because they lack age composition. SSB, recruitment, abundance, and N/F-at-age can be compared on the common model grid. Fixed M is shown as an input, not treated as evidence of fitted agreement. Native arithmetic Fbar is not like-for-like with tinyAM's abundance-weighted Fbar.

The default-start fit converged (optimizer code 0, objective 302.221, maximum absolute gradient 0.000117, successful standard-error estimation, and a positive-definite Hessian). No starting values were supplied. The fit emitted a support warning because F-increment age correlation was near a boundary. Trend correlations were high for N-at-age (0.966), F-at-age (0.992), SSB (0.965), recruitment (0.974), and total abundance (0.958), but scale differences were large: mean absolute percent differences were 91.6%, 38.3%, 71.5%, 84.1%, and 59.1%, respectively. Fixed M is identical by construction and is not counted as evidence of agreement. These results indicate good trajectory shape but substantial scale disagreement under this simplified translation; no settings were tuned to the official estimates.
