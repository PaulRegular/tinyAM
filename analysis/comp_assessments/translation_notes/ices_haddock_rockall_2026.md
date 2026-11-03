# ICES Rockall haddock, 2026 assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_rockall_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_rockall_2026.md).

## Expanded-library review

### tinyAM fit and comparison review

The retained fit uses package-default starting values and converged at database revision `6b33b16`: optimizer code 0 (`relative convergence (4)`), objective 679.281, maximum absolute gradient 0.001647, successful standard-error calculation, and a positive-definite Hessian. It has 18 fixed effects and 552 random effects. No explicit starting parameters were used. An initial preparation attempt exposed an unused age-1 N-SD category under exponential N0 initialization; the final model uses one shared N-process SD for the active age-2–8 survival states. A diagnostic advisory remains: F-increment age correlation is within 0.05 of a boundary, so its uncertainty and age-correlation structure need caution.

The dashboard is cached locally at `results/cache/ices_haddock_rockall_2026/dashboard.html`. The most relevant comparisons show high temporal correlations for F (0.949), Fbar (0.966), N-at-age (0.961), and approximate common-definition SSB (0.954). Recruitment is approximate and has 27.1% mean absolute percent difference (correlation 0.908). Total abundance has 18.9% mean absolute percent difference, but its terminal difference is 60.7% and correlation is 0.881. Approximate SSB's terminal difference is 73.5%; this uses the catch-weight proxy because the full official stock-weight surface is unavailable. These differences are visible in the dashboard and should be reviewed; fixed-M agreement is exact by construction and is not evidence that the fit recovered M.
