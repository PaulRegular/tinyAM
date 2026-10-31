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

## Remaining validation

Visual inspection, recovery experiments and source-supported multi-fleet
assessment trials are recorded as they finish on the subsequent branches.
No recovery or superiority claim is made from the artificial workflow example.
Generated fits, warnings, dashboards and numerical tables stay in ignored
`results/multi_fleet/` and `analysis/sam_equivalence/results/` directories.
