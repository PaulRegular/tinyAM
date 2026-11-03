# Southern Gulf of St. Lawrence 4T American plaice: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_dfo_platessoides_4t.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_dfo_platessoides_4t.md).

## Expanded-library review

### tinyAM translation

The first translation uses the source years 1976–2015 and ages 4–16+, with commercial catch, both abundance indices, source weight and maturity surfaces, and estimated-M age groups. tinyAM's population and fishing-mortality processes and its catch likelihood simplify the VPA; it does not reproduce the source's VPA catch equation, M random-walk priors, age-specific q and observation-variance parameters, or MCMC integration. The 2003 RV index remains omitted. The 2003 age-16+ zero weight is preserved as supplied when preparing observations and will be evaluated as a potential source limitation rather than silently changed.

The first fit used package-default initialization, with no explicit starting values. It passed the numerical convergence checks (optimizer code 0, maximum absolute gradient 0.00139, successful standard-error calculation, positive-definite Hessian). Comparing the cached fit with the newly recovered outputs gives mean absolute percent differences of about 661% for N-at-age, 66% for total abundance, 89% for age-4 recruitment, and 90% for the digitized native SSB series; trend correlations are 0.865, 0.947, 0.944, and 0.962, respectively. The grouped mortality comparison uses official-N-weighted means of tinyAM ages. F has about 70% mean absolute percent difference, a -91% terminal difference, and 0.914 trend correlation; M has about 66%, a -76% terminal difference, and 0.453 trend correlation. These figure-based rate comparisons are approximate and F coverage is incomplete. The fit is numerically usable, but its mortality and population scale do not closely reproduce the VPA.
