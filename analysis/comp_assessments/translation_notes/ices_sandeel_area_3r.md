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

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_sandeel_area_3r.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

A finite latent log-F state cannot reproduce the structural zero in 2025.
The annual translation must preserve this terminal year and disclose that
limitation; adding small effort/catches would change source observations.
The source's zero F also has no defined percent-difference denominator, so
an unavailable terminal F statistic cannot be interpreted as agreement.

The first annual translation retains 1987–2025, ages 1–4+, both surveys
(Dredge age 1 and Acoustic ages 1–4+), the native q/SD sharing among retained
observations, January biology and supplied annual M. Seasonal catch counts
are summed and converted from million to thousand fish. July age-0 Dredge
observations are omitted. Catch weights are observed-catch-weighted seasonal
means; the four zero-catch cells in 2025 use equal-duration means for yield
reporting only. The zeros themselves are unchanged.

Free initial abundance, deterministic older-cohort survival and IID age-1
recruitment approximate the source's January adults and weak July recruitment
penalty. Correlated random-walk F increments were first used to approximate
common annual fishing variation without taking the log of zero effort.
This does not reproduce exact effort proportionality or shared F states.

| Initialization | Optimizer code | Objective | Maximum absolute gradient | Positive-definite Hessian |
|---|---:|---:|---:|---|
| Package defaults | 1, false convergence | 362.1409 | 503,245 | No |
| Aligned source starts | 1, false convergence | 361.8999 | 124,547,007 | No |

Both failed numerical checks. An sdreport object was returned, but its
covariance was not finite/usable. The estimated F increment age correlation
was near its boundary; some residuals exceeded three SDs. A single
source-start retry used positive official F cells and aligned N states;
the unrepresentable zero-F cells retained package log-F starts of zero
(F = 1), not an invented positive source observation. Settings/data were
identical between attempts.

The resulting exploratory comparison has large differences and must not be
treated as a converged result. Terminal official F is zero in all four ages,
so those percent differences remain unavailable while absolute differences
remain visible. The failure diagnostics, source-start attempt and exploratory
dashboard are cached. This motivates testing independent annual F increments
as the simpler existing alternative identified by the correlation diagnostic;
it does not establish that the closure alone caused the failure.

Nine targeted scripts pass: the three regional source and translation
checks, reference alignment, plus groups and default-first initialization.
They protect the original survey values, q/SD sharing, annual transformations,
time-varying supplied M, mixed recruitment dates, terminal zeros and missing
uncertainty. Dashboard metadata were refreshed from the cached fit only after
checking that all five fitted data tables were unchanged. The rendered HTML
contains the scientific-review warning and the structural-zero limitation.

The full package suite passes **2,929 checks, zero failures, warnings or
skips**. Its first attempt could not locate Pandoc; rerunning with the bundled
renderer configured resolves those environment-only dashboard failures.
`R CMD check --no-manual --no-tests` returns **Status: OK**; tests ran
separately as above. Logs are retained beside the dashboard. The broader
assessment-analysis sweep was not repeated for this milestone; its previously
documented unrelated baseline failures are not claimed to be resolved.
