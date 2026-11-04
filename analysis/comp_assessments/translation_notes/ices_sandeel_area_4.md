# Sandeel in Area 4: detailed 2025 assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_sandeel_area_4.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sandeel_area_4.md).

## Expanded-library review

### Annual translation and fitting outcome

The recipe uses January ages 1–4+ in 1993–2024, free initial abundance,
deterministic older-cohort survival, IID age-1 recruitment and independent
random-walk annual F increments. It sums seasonal M hazards and catch
numbers, converts million to thousand fish, and uses January biology with
spawning time zero. This preserves annual older-cohort survival for given
annual F/M, but simplifies the within-year catch/survey equations and does
not reproduce effort-proportional/shared F or exact structural zeros.
The survey SD floor and profiled seasonal catch-SD construction are not
reproduced. No fitted values alter source observations.

The juvenile survey is retained through an observation-mean re-expression,
without changing any index values. For July recruitment R0, constant supplied
second-half M0 and zero age-0 F, define

`A(M0) = mean(exp(-M0*t)), t in [0.75, 1]`.

The source prediction is `q0 * R0 * A(M0)`. Following January abundance is
`N1 = R0 * exp(-M0)`, so the same prediction is `q_star * N1`, where
`q_star = q0 * exp(M0) * A(M0)`. Because M0 is constant and each era/age has
its own free q, this absorbs a constant into q and preserves the conditional
log-observation likelihood exactly. It does not require source fitted N/F/q
to transform observations. Tests verify both expected values and likelihoods.

Original age-0 observations are linked to the next January's age-1 state at
model time zero. Their actual survey year and age remain `source_year` and
`source_age` in the data/dashboard. Original age-1 observations retain their
year at autumn midpoint 0.9375. Four series preserve era/source-age q and SD
sharing. There are 42 numerical observations: only the final 2024 juvenile
index is omitted because its January survivor state is outside the fitted
window; the original missing 2003 age-1 cell is not filled. This mapping is
conditional on constant M0 and zero age-0 F and is guarded accordingly.
It is **not** applied to reference recruitment: comparisons still use official
January N1 in the same calendar year, with native July recruitment separate.

Catch reporting uses observed-count-weighted seasonal weights for positive
cells and equal-duration weights for zero cells, only for predicted-yield
reporting. Zero-coded age arrays with positive SAG landings do not become
zero observed biomass: those six annual biomass observations remain missing.
Verified zero-landings years 2009/2024 retain zero. Source conditional Yield
is not compared as though it were unrestricted tinyAM yield.

Both package-default and the single aligned source-start attempt failed
with **NA/NaN gradient evaluation**, taking about 3.34 seconds together in the
earlier attempt. A recheck on 2026-10-26 against database revision `b73231f`
failed in both attempts again: 1.47 seconds from package defaults and 2.02
seconds from the aligned source start. Both retained missing optimizer codes,
objectives and gradients, and neither produced a successful `sdreport`.
No optimizer result, final objective/gradient, usable sdreport or `tam_fit`
object was returned. Source starts used aligned January N and positive
annual F cells; unrepresentable zero F retained package log-F starts of zero
(F = 1). The data and settings were identical between attempts. No further
starting-value or assumption tuning was used to force convergence.

The exact cause of the non-finite gradient has not been isolated. Missing
support for the source's structural-zero/shared F representation and survey
SD floor are documented differences, not a demonstrated causal diagnosis.
This case is a verified import and tested runnable translation with actual
failed fitting attempts, **not a successful replication**. SSB, recruitment,
N/F-at-age and abundance agreement statistics are unavailable, not passes;
fixed M is an input consistency check, not estimated-M agreement.

The local `results/cache/ices_sandeel_area_4_2025/dashboard_official_only.html`
shows recovered Official outputs and translated inputs, explicitly without
tinyAM estimates. It supplies no invented parameters, standard errors or
optimizer. The missing yield comparison is labelled unavailable. Two narrow
reporting fixes permit this: a missing-yield guard in the dashboard template,
and complete projection/missing-q metadata for standalone references.
Cached initialization attempts, failure diagnostics, the reference object
and original-to-model index mapping remain available. No core mathematics
or public API was changed and shared aggregate tables were not regenerated.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_sandeel_area_4.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The four regional source checks, four translation checks, reference alignment,
plus-group handling and default-first initialization checks all passed (11
assessment-analysis scripts). The dashboard tests passed 203 assertions;
the full package suite passed 2,936 assertions with no failures, warnings or
skips. A package-only `R CMD check --no-manual --no-tests` returned Status: OK;
the complete package tests were run separately. The full assessment-analysis
suite was not repeated for this milestone; its previously documented unrelated
failures remain outside this source/translation review. The actual Official-only
Area 4 dashboard rendered successfully.
