# Sandeel in Area 3r: HAWG 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_sandeel_area_3r.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sandeel_area_3r.md).

## Expanded-library review

### Retained independent-increment trial

Replacing correlated F increments with independent random-walk increments
removes the weak correlation parameter, leaving observations, years, biology,
q/SD sharing and all other assumptions unchanged. This trial passed numerical
checks from package defaults in 3.63 seconds: optimizer code 0, objective
373.6773, maximum absolute gradient 0.00000711 and a positive-definite
Hessian (15 fixed, 194 random effects). No source-start retry is permitted
merely because agreement is poor, and none was used.

**The result is not scientifically usable.** Estimated F is only
0.00000000093–0.000000031 per year over the entire fitted surface, while
abundance is inflated by many orders of magnitude. Recruitment/SSB and q
log-scale SEs are about 5,058, giving intervals from zero to infinity.
There are 36 fixed-parameter correlations above 0.95 and residuals above
three SDs. A positive-definite Hessian and small gradient do not establish
that the abundance/q scale is adequately identified. The exact source effort
relationship has been relaxed, and annual aggregation changes the seasonal
survival/observation equations; these are plausible contributors, not a
verified attribution of cause solely to the closure or process option.

| Approximate metric | Mean absolute % difference | Terminal % difference | Trend correlation |
|---|---:|---:|---:|
| SSB | 5.365 billion | +3.832 billion | 0.574 |
| January total abundance | 2.716 billion | +3.460 billion | 0.931 |
| January age-1 recruitment | 2.202 billion | +3.501 billion | 0.939 |
| N at age | 6.473 billion | +2.870 billion | 0.793 |
| Annual F at age | approximately 100 | unavailable (official zero) | 0.793 |

Every comparable core metric is flagged. Native July recruitment and January
age-1 recruitment have different definitions; the table compares the latter
under a common definition. Fixed supplied M matches every cell, an input
check only. Official standalone biomass-at-age and joint uncertainty for
derived aggregates are unavailable. No failed or missing comparison is
treated as a pass.

The reproducible recipe and **exploratory** dashboard are retained at
`results/cache/ices_sandeel_area_3r_2026/dashboard.html`. The earlier failed
correlated trial remains in `cor_rw_failed/`, with its dashboard at
`dashboard_exploratory.html`. Current fit/initialization diagnostics and the
review result are cached alongside the dashboard. A non-finite SSB interval
produces a prominent dashboard warning without overriding the numerical
convergence flag or triggering extra initialization trials.

This case counts as verified source recovery and an actual translation/fit
attempt, not a successful official–tinyAM replication. It should not be used
as evidence of estimated abundance or mortality agreement. No source values,
terminal years, tinyAM mathematics or aggregate library results were changed
to conceal the problem.
