# Multi-fleet validation

This development sequence targets tinyAM 0.2.0. Fleet support is optional;
the development version is 0.1.1.9001. Scoped formulas, new curves, smooths and
the planned recovery studies follow the fleet foundation on separate branches.

## Fleet foundation

The implementation adapts the independent-fleet calculations from this
repository's parked `multi-fleet` branch (`e975b58`, `66b169d`). It retains
current-main recruitment, F blocks, typed weights and evaluation workflows;
the parked branch was not merged. No settled stock hypotheses were changed.
F covariates now belong to catch inputs, including single-fishery models.

Deterministic tests cover fleet Baranov catches, independent densities for
all five F process options, per-fleet blocks and structural zeros, shared and
separate correlation groups, dense covariance checks, finite derivatives,
typed fleet catch weights, retained plus ages and catch-target projections.
One labeled fleet is compared with the original single-fishery calculation.

An independent comparison with main at `47100da7` found exactly zero
differences in joint objectives, gradients and predictions for off, IID, RW,
correlated-RW and AR1 F examples after moving their F covariates. The R package
code at that revision equals the starting GitHub main (`3d1ac69a`).

The pinned native-SAM `nscodbiopro` baseline passes all 33 existing comparisons:
maximum matched joint-objective difference 2.28e-13, scaled gradient difference
6.81e-13, scaled Hessian difference 1.60e-12 and Laplace difference 6.83e-13.
Optimized, perturbed-start and two retrospective comparisons pass their existing
tolerances. These are numerical-equivalence checks, not bitwise-identity claims.
Source revisions and test definitions remain in `analysis/sam_equivalence`.

All 30 unblocked stock recipes prepare without fitting. Five remain blocked
by existing source limitations: BSAI arrowtooth, AI cod, Strait of Georgia
herring, BC petrale sole and BC Pacific ocean perch. Southern Gulf cod and
Atlantic mackerel F predictors were moved to catch without changing values.

`multi_fleet_validation.R core` uses an explicitly artificial two-fleet split
of packaged catch to exercise fitting, uncertainty, warm starts, simulation,
catch-conditioned projections, hindcasts and dashboard rendering. Its SCAA
fit converged (maximum gradient 0.000291, positive-definite Hessian); the
one-fold hindcast was scored successfully, and both projected yield targets
were reached. This is a workflow check, not a multi-fleet recovery result.

An additional active-F rendering check used the same artificial split with
different fleet blocks and separate correlated-RW SDs/correlations. It failed
numerical checks (false convergence, maximum gradient 0.437, non-positive-definite
Hessian). Its warnings and fit were retained; its dashboard rendered, with
unavailable joint-sample residuals identified. No tuning was applied.

The full package suite in `R CMD check --no-manual` passed 3,943 checks, with no
failures, warnings or skips; the final package check returned **OK**. Focused
F-block, fleet and process-display checks passed 428 assertions.
Six additional hand-calculated fleet-hindcast assertions passed after that run,
checking actual removals, missing targets and separate training-weight averages.
The 69
comparative-assessment scripts returned 59 passes and the same ten failures
documented on starting main: age subsets, Atka, BC canary, coastal hake, EBS
cod/pollock, GOA cod/pollock, Pacific M priors and recipe blocks. No canonical
assessment exports were retained from testing.

## Scoped formulas

`for_fleet()` and `for_survey()` compile each inner formula on its selected
rows before combining it with shared global terms. Tests cover contrast and
formula environments, disjoint Gaussian mean effects, different fleet curves,
sparse survey support, abundance powers, projection preparation and warm starts.
The targeted suite passed 525 assertions with no failures or warnings.
The full package check returned **OK** with 3,995 passing assertions and no
failures, warnings or skips (`scoped-model-formulas`, commit `6b545481`).
The comparative suite returned 59 passes and the same ten starting-main
failures listed above. Two worker-based checks needed the ordinary local
worker connections that were blocked by the restricted test runner; they
passed when rerun with those connections available.

`multi_fleet_validation.R scoped` compares equivalent ordinary and scoped
fixed designs on the artificial catch split and two copies of the survey.
Both fits converged with positive-definite Hessians (maximum gradients 0.00229
and 0.00251). Both had 55 fixed parameters, no random coefficients, and
AIC/cAIC 1162.469; their cAIC difference was 7.49e-9. The scoped one-fold
hindcast was scored, and explicitly requested OSA residuals left the supplied
fit's parameter values unchanged. The two-model dashboard rendered without
warnings. These duplicated observations test software behavior, not recovery
or the benefit of collecting an additional survey.

Browser-based visual inspection was attempted but the computer-use tool stopped
because it could not verify the current Chrome URL. HTML generation succeeded;
interactive visual inspection remains outstanding.

## Selectivity curves

The new curve tests independently check warped half-selection, slope, peak,
tail, Gaussian/asymptotic endpoints, first and second AD derivatives, and
WHAM's raw four-parameter double-logistic values and derivatives. WHAM source
is pinned to `ffcbfd093161e8210aa20716933975191b661807`; its `get_selAA()`
age counter is translated explicitly when comparing biological ages.
No peak normalization or midpoint ordering is added. Tiny positive rates stay
on the log scale for likelihood calculations, without a positive floor.
The 123 new curve assertions and 391 workflow assertions pass without warnings.
Rechecking the five one-fleet baseline processes found zero objective changes,
maximum gradient difference 7.11e-15 and prediction difference 1.11e-16.

`multi_fleet_validation.R curves` generates a 20-year, two-fleet example with
logistic and warped F, IID recruitment and supplied M. It uses default starts
and optimizer settings, seed 47281, and no rescue tuning. The initial fit hit
the 1,000-function-evaluation limit (945 iterations, maximum gradient 5.70),
despite a positive-definite Hessian. Both hindcast attempts failed numerical
checks and supplied no forecast scores. The attempted fit, warnings, truth,
projections, simulations and rendered dashboard are retained. A read-only
derivative review agreed with independent central differences; the largest
absolute discrepancy for the six checked curve/initial parameters was
7.52e-5 (the corresponding gradient was 5.70).

This attempt demonstrates an unresolved estimation limitation with many annual
fleet coefficients and nonlinear curves. Close parameter recovery in one
unconverged fit is not evidence of reliable estimation. The planned recovery
experiments must retain this difficulty rather than tune it away.

The final curve package check returned **OK** with 4,118 passing assertions,
no warnings/failures/skips, after fixing log-residual extraction from cached
tables without SD columns. That regression and the underflow case pass 161
focused assertions together. Comparative assessments returned the same 59
passes and ten documented failures. The curve milestone is pushed in
`c21ee870`, `6ea3c2fc` and `a3ad9881`.

## Penalized smooths

Univariate thin-plate/cubic bases, categorical-by smooths and three-penalty
age-year tensors agree with mgcv's matrices. Independently constructed dense
Gaussian covariances reproduce integrated objectives and gradients. Tests also
cover retained-basis predictions, training-only forecasts, compatible warm
starts, scopes, coherent simulation, signed null-space coefficients and full
function uncertainty. The expanded targeted run passed 581 assertions without
failures or warnings, including curve regression and warm-start checks.

`multi_fleet_validation.R smooths` generates and fits a separate, fully simulated
smooth model (seed 47281). It converged with a positive-definite Hessian and
maximum gradient 0.00329. Two smoothing SD intervals span more than tenfold;
these limitations remain flagged. The fit has 52 fixed and 27 random parameters.
Its cAIC accounting reconciles: recruitment EDF 18.72865 plus smooth/other EDF
3.948775 gives random EDF 22.67742; adding fixed EDF 47 gives EDF 69.67742.
This is an accounting check, not a model ranking.

Catch-target projections, coherent simulations, explicit OSA extraction and
dashboard rendering completed. The one hindcast peel failed because its Hessian
was not positive definite and its fixed covariance was unavailable/non-finite;
it was retained and excluded from scoring. The unpeeled reference passed.
No rescue tuning or automatic component removal was applied.

The smooth milestone's full package check returned **OK**, with 4,218 passing
assertions and no failures, warnings or skips. Comparative assessments again
returned 59 passes and the same ten pre-existing failures listed above.

## Recovery experiments

`multi_fleet_recovery.R` runs 100 isolated replicates for each of 13 curve/smooth
configurations and 30 full assessments for each of 13 scenarios: **1,690
attempts**, all retained. `postprocess` adds N-at-age uncertainty to a private
copy of each successful objective without optimizing; `metrics` produces
replicate-level recovery/coverage and Monte Carlo SE tables. The original
manifest retains generation hashes, seeds and runtime versions. Subsequent
changes added only reporting/uncertainty extraction and corrected documentation,
not the data generator, likelihood or fitting controls. A fresh results
directory is required if the script's provenance changes; existing results are
never silently reused under changed sources.

Full assessments use 12 training years, ages 1–10, two fleets, known M, distinct
catch weights, and two held-out years. The SCAA has logistic/double-logistic
fleet curves; SAM scenarios use shared/separate age correlations and fleet F
blocks (one unfished young age). Warped shapes are asymptotic, intermediate and
Gaussian, each with full or sparse older-age training coverage. Smooth scenarios
use age/year means with off/IID/AR1 residuals; the weak-separation scenario has
more temporal basis functions and stronger residual persistence. Penalized
smooth SDs are 0.2/0.15 on the retained mgcv basis scale, so their signal can be
weak even when the residual process is off. Survey q is logistic throughout.

No starts, optimizer limits or model assumptions were adjusted after seeing
results. Forecast truth holds terminal fleet F; it is a conditional
status-quo-F experiment, without new F innovations, not a test of
process-integrated forecast uncertainty.

| Full scenario | Numerically successful / 30 |
|---|---:|
| SCAA: logistic/double logistic | 30 |
| SAM: shared correlation | 30 |
| SAM: separate correlations | 30 |
| Warped: asymptotic, full / weak coverage | 0 / 0 |
| Warped: intermediate, full / weak coverage | 3 / 6 |
| Warped: Gaussian, full / weak coverage | 0 / 11 |
| Smooth: off / IID / AR1 | 7 / 11 / 21 |
| Smooth: weak mean/AR1 separation | 20 |

Thirty successes out of thirty has an exact 95% binomial interval of
88.4–100%; zero of thirty has 0–11.6%. The observed 100% rates do not establish
a population success rate above 90%. For smooth off/IID/AR1, the plug-in
success-rate MCSEs are 0.077/0.088/0.084; weak separation has 0.086.
There were no uncaught study errors. Numerical failures include optimizer
limits, indefinite Hessians and unavailable/non-finite covariance. Scientific
advisories and all warnings are stored separately.

Recovery summaries below **condition on successful fits**. Bias is the mean
relative error across historical cells within each replicate, then across
replicates. The ± values are between-replicate Monte Carlo SEs, in percentage
points; they do not count ages or years as independent replicates.

| Scenario | F bias % | N bias % | SSB bias % | Recruitment bias % | 95% SSB coverage % |
|---|---:|---:|---:|---:|---:|
| SCAA | −0.16 ± 0.23 | 0.31 ± 0.24 | 0.31 ± 0.24 | 0.22 ± 0.28 | 99.2 ± 0.6 |
| SAM shared | 0.25 ± 0.69 | 0.19 ± 0.62 | −0.02 ± 0.63 | −0.07 ± 0.80 | 93.1 ± 1.7 |
| SAM separate | −1.04 ± 0.66 | 1.31 ± 0.61 | 0.83 ± 0.64 | 2.77 ± 0.76 | 94.7 ± 1.7 |
| Smooth off (7) | 0.19 ± 0.44 | −0.02 ± 0.26 | 0.10 ± 0.29 | 0.06 ± 0.41 | 98.8 ± 1.2 |
| Smooth IID (11) | 0.29 ± 0.20 | −0.01 ± 0.33 | −0.15 ± 0.34 | 0.41 ± 0.33 | 93.9 ± 3.9 |
| Smooth AR1 (21) | 0.48 ± 0.18 | −0.18 ± 0.12 | −0.32 ± 0.11 | 0.01 ± 0.32 | 98.0 ± 1.4 |
| Smooth weak (20) | 0.99 ± 0.26 | −0.48 ± 0.31 | −0.53 ± 0.34 | 0.28 ± 0.26 | 91.7 ± 3.2 |

Warped recovery is based on only 3–11 retained fits and cannot offset the poor
success rates. Their full surface, shape, scale and interval results remain in
the generated tables. Among the three successful full-coverage intermediate
fits, mean c was 0.3461 (MCSE 0.0018; truth 0.35); six weak-coverage fits gave
0.3486 (0.0034). These selected estimates are not an unconditional recovery claim.

SAM correlation recovery on the fitted atanh scale had bias 0.041 ± 0.041
(shared), 0.038 ± 0.041 (Gillnet) and −0.112 ± 0.066 (Trawl), with interval
coverage 100%, 96.7% and 93.3%. Smooth/AR1 separation was poor: smoothing
log-SDs often shrank toward −9 despite generating values around −1.6 to −1.9.
Very wide intervals could still cover the truth. In the deliberately weak
scenario, age/year AR1 logit-parameter biases were −0.718 ± 0.093 and
−0.537 ± 0.101; coverage was only 65%/70%. Good population recovery among
converged fits does not establish identification of its component processes.

Isolated fits had 100/100 successes for logistic, both double-logistic midpoint
orders, intermediate warped shapes, and tp/cr smooths; asymptotic warped shapes
had 99/100, factor-by/tensor smooths 97/100, and estimated Gaussian-endpoint
warped shapes 0/100. One hundred successes has exact 95% interval 96.4–100%;
zero has 0–3.6%. Function coverage was 95.4% for ordinary double logistic but
85.3% for reversed midpoints; tp/cr/factor-by/tensor smooth coverage was
87.0%/93.9%/89.2%/92.0%. Local Gaussian/Wald intervals can under-cover.

A read-only Gaussian failure review found a flat fitted curve with a50=38.8
and log scale=23.16, objective 44,138.76, compared with −39.73 at the known
truth (a50=3.5, scale=1.5). Independent formula evaluation agrees at both
points; the truth has positive curvature, while the fitted Hessian is nearly
singular/indefinite. This is an unresolved default-start/local-optimization
failure, not evidence that Gaussian shapes are intrinsically unidentified.
No alternative start was fitted.

One-/two-year survey log-RMSE averaged 0.130/0.167 for SCAA, 0.187/0.214 for
SAM shared, and 0.171/0.204 for SAM separate (MCSE 0.005–0.010). These are
different generators, not a model ranking. Catch-proportion RMSE averaged
0.0172/0.0240/0.0229 respectively (MCSE 0.0007–0.0008). MASE is particularly
unstable with one origin per replicate and nearly zero persistence errors:
SCAA's one-year mean was 359.7 ± 358.1. All individual cases/denominators and
coverage remain available; this statistic should not drive conclusions here.
Predictive interval coverage and process-integrated scoring remain future work.

## Historical herring trial and regression checks

`multi_fleet_validation.R herring <database-directory>` uses the recovered
2020 west-of-Scotland/Ireland inputs: 1957–2019, ages 1–9+, both catch fleets,
their own catch weights, stock biology, fixed M and three surveys. It retains
the established IID recruitment/N assumptions, adds one correlated F random
walk per fleet with fleet-specific SD/correlation parameters, and estimates
fleet catch SDs. Survey q/SD groups follow the recovered tables; no rescue
tuning or settled-recipe changes were made.

The fit passed numerical checks: positive-definite Hessian, maximum gradient
0.00397. Source data, warnings, fit and a self-contained dashboard with sparse
Official references are cached under `results/multi_fleet/herring/`.
`herring-dashboard` rerenders without refitting. Area-specific source F is
explicitly mapped to its corresponding catch stream and remains an approximate
comparison; unlabeled stock-wide F is not invented from those components.
The reference converter now preserves fleet F and 3-D report identities, with
regressions for missing cells and retained plus ages.

Source limitations remain: inferred catch units, unresolved survey units,
unknown detailed source F bindings, survey-observation covariance unsupported
by tinyAM, and no executed native source fit. The trial is neither exact
equivalence nor evidence of guaranteed improvement. Local source metadata had
literal encoded CSV record separators; the trial decoded those in memory only.
Unavailable index cells stay missing. The original files and canonical exports
were not modified. Preparation/render errors are retained in the cached history.

No new official native multi-fleet SAM comparison is claimed: this recovered
source has unsupported observation covariance and no verifiable executed fit.
Existing supported native SAM regressions were rerun: all 132 likelihood/fit
checks and six SAMjr printed-output comparisons pass. The haddock fit's
objective difference is 5.00e-12, with convergence and positive-definite
uncertainty. Printed-output agreement is not bitwise identity.

## Remaining validation

Visual inspection, recovery experiments and source-supported multi-fleet
assessment trials are recorded as they finish on the subsequent branches.
No recovery or superiority claim is made from the artificial workflow example.
Generated fits, warnings, dashboards and numerical tables stay in ignored
`results/multi_fleet/` and `analysis/sam_equivalence/results/` directories.
