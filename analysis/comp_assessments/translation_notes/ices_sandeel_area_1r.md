# Sandeel in Area 1r: HAWG 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_sandeel_area_1r.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sandeel_area_1r.md).

## Expanded-library review

### Retained annual tinyAM translation

The stock recipe fits 1983–2025, January ages 1–4+, without changing the
canonical source rows. It sums the two half-year catch counts at each age,
converts million to thousand fish, and adds the integrated seasonal M hazards.
January weights and maturity define the stock and spawning biomass. Catch
reporting uses observed-catch-weighted seasonal catch weights; it excludes
age-0 catches rather than comparing an ages-1+ prediction with all-age SAG
landings. The Dredge age-0 index is omitted because its July recruitment state
does not exist in this annual model. Dredge age 1 uses the midpoint 0.9375 of
the source averaging interval; RTM ages 1–3 retain configured time zero.

Both indices retain their published numerical values. Their source units
identify scaled indices with unstated absolute multipliers; analysis labels
make this explicit so the shared converter does not normalize RTM values as
absolute fish counts. An initial unit-review fit did normalize RTM to thousands;
the corrected fit retains source scale. As expected with freely estimated
survey-age q, population results are essentially unchanged. Both fits remain
in the ignored local review cache.

Older-cohort N has no process error and initial January abundance is free.
Age-1 recruitment is IID about a fitted mean, a simplification of weakly
penalized July recruitment and its hockey-stick expectation. F has IID states
about period-by-age-group means, with source changes in 1999 and 2022 and a
fitted log annual-effort covariate. Ages 3–4+ share their mean, not their latent
F states. This approximates the source seasonal effort relationship; it is
not an exact proportional constraint. Catch SDs share ages 1–2 and 3–4+;
survey q is separate by age and survey, with RTM ages 2–3 sharing an SD.

The corrected model converged from package defaults, with optimizer code 0,
objective 251.8868, maximum absolute gradient 0.0001363 and a positive-definite
Hessian. Standard errors were obtained. There are 26 fixed and 214 random
effects. No aligned source-start retry was needed. These numerical results
do not remove the source gaps or establish equivalence to the seasonal model.
tinyAM integrates recruitment/F states as random effects. That does not
reproduce the native effort-coefficient construction or certify the source's
exact likelihood integration from the incomplete public snapshot.

Comparisons use the same calendar-year January N at age 1 as recruitment,
not a shifted July age-0 series. Annual F is compared once, excluding duplicate
seasonal surfaces. Fbar is recalculated with January N weights over ages 1–2.
Native July recruitment and native aggregate outputs remain source context.
Shared-biology biomass/SSB and derived abundance receive no invented uncertainty.

| Metric | Mean absolute % difference | Median absolute % difference | Terminal % difference | Trend correlation |
|---|---:|---:|---:|---:|
| SSB | 12.24 | 7.51 | +17.99 | 0.968 |
| January total abundance, ages 1–4+ | 14.65 | 9.83 | +29.59 | 0.956 |
| January age-1 recruitment | 22.70 | 16.44 | +37.88 | 0.965 |
| N at age | 22.06 | 13.41 | +20.07 | 0.969 |
| Annual F at age | 8.33 | 6.75 | +19.94 | 0.988 |
| Common-definition Fbar | 6.03 | 5.99 | +3.75 | 0.993 |

All these comparisons are approximate. Recruitment, N-at-age and terminal
total abundance exceed the existing review thresholds and remain flagged;
they were not tuned to reach those thresholds. Fixed annual M agrees in
every cell, an input-consistency check rather than evidence of M estimation.
Standalone official biomass-at-age is unavailable. Age-specific correlation
is the mean within-age temporal correlation, not a correlation across the
age gradient.

The local dashboard is
`results/cache/ices_sandeel_area_1r_2026/dashboard.html`. Fit, initialization,
comparison and validation records are retained alongside it. Shared aggregate
result tables and plots have not been refreshed from this single-stock run.
