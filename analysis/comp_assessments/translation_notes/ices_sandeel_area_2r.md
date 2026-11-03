# Sandeel in Area 2r: HAWG 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_sandeel_area_2r.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sandeel_area_2r.md).

## Expanded-library review

### Retained annual tinyAM translation

The translation fits calendar years 1983–2025 and January ages 1–4+.
It sums seasonal catch numbers, converts million to thousand fish, and adds
the two integrated M hazards. Annual M is 0.98, 0.92, 0.71 and 0.64 at
ages 1–4+, respectively. First-half weights and maturity define January
biomass and SSB, with spawning time zero. Catch reporting uses observed
catch-weighted seasonal weights for ages 1–4+, without comparing that subset
to all-age SAG catches. Canonical values are unchanged by preparation.

The fit retains the 22 Dredge age-1 values, on their original numerical scale,
at the midpoint 0.9375 of the source prediction interval. It omits the
22 July age-0 observations because this annual model starts at age 1.
There is one constant q and one observation SD. Catch SDs share ages 1–2
and 3–4+; seasonal SD distinctions are not reproduced after summing catch.

Initial abundance is free, older-cohort N has no process error, and age-1
recruitment is IID about a mean. This changes the July hockey-stick
construction and its likelihood weight. F has IID states about period/age
means and a fitted log annual-effort effect. The four source periods and
shared age-3/4+ mean are retained; exact proportionality to seasonal effort
and identical latent F states are not. No observations are adjusted using
source fitted parameters.
tinyAM integrates recruitment/F states as random effects; the incomplete
public snapshot does not establish equivalence of likelihood integration.

The fit converged from package defaults in about four seconds: optimizer
code 0, objective 247.6649, maximum absolute gradient 0.00001853,
successful standard errors and a positive-definite Hessian. It has 24 fixed
and 214 random effects. No source-start retry was needed. Numerical
convergence does not establish equivalence to the seasonal assessment.
There is one residual advisory: an age-4 catch residual reaches about 3.16
standard deviations. It is retained without changing observations or tuning
the model. No structural redundancy or numerical failure is detected.

| Metric | Mean absolute % difference | Median absolute % difference | Terminal % difference | Trend correlation |
|---|---:|---:|---:|---:|
| SSB | 10.10 | 9.74 | +21.95 | 0.990 |
| January total abundance, ages 1–4+ | 13.28 | 10.51 | +20.73 | 0.993 |
| January age-1 recruitment | 19.77 | 14.34 | +20.18 | 0.994 |
| N at age | 19.79 | 14.11 | +23.98 | 0.993 |
| Annual F at age | 8.87 | 7.51 | +8.89 | 0.961 |
| Common-definition Fbar | 3.24 | 2.19 | +3.42 | 0.991 |
| January biomass | 11.72 | 8.88 | +21.01 | 0.990 |

All are approximate comparisons. SSB, recruitment, total abundance and
N-at-age remain flagged for terminal differences above 20%; the secondary
biomass comparison also exceeds that threshold. The model was not tuned to
the thresholds. Age-specific correlations are means of within-age temporal
correlations. Fixed M agrees in every cell, an input check rather than an
estimated-M result. A small Mbar difference arises solely from different
N weights. Official standalone biomass-at-age is unavailable.

Recruitment comparisons use same-calendar-year January N at age 1, without
shifting the July age-0 series. Native recruitment and aggregate estimates
remain source context. Common-definition aggregates and January age-1
recruitment have no invented uncertainty. The local dashboard is
`results/cache/ices_sandeel_area_2r_2026/dashboard.html`; initialization,
fit, comparison and validation records are cached alongside it. Shared
aggregate results have not been refreshed from this single-stock run.
