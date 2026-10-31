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

## Remaining validation

Visual inspection, recovery experiments and source-supported multi-fleet
assessment trials are recorded as they finish on the subsequent branches.
No recovery or superiority claim is made from the artificial workflow example.
Generated fits, warnings, dashboards and numerical tables stay in ignored
`results/multi_fleet/` and `analysis/sam_equivalence/results/` directories.
