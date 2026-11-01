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

Logistic and WHAM-style double-logistic curves remain supported for F means
and catchability. Independent tests compare the double-logistic raw product
and AD derivatives with WHAM revision
`ffcbfd093161e8210aa20716933975191b661807`, including its age-counter offset
and reversed midpoints. No peak normalization or midpoint ordering is added.
Tiny selectivities stay on the log scale without artificial positive floors.
Shared tests cover scoped fleets/surveys, simulation and projections.

### Retired warped-normal exploration

`warped_normal()` was retired before release because default fitting was
unreliable. Diagnostics found poor optimizer solutions, nearly flat curves at
extreme widths and unresolved boundary uncertainty. The complete code, results
and diagnosis are recorded at [revision 417e2636](https://github.com/PaulRegular/tinyAM/tree/417e2636/analysis/formula_expansion).
Original caches, warnings, profiles, manifests and benchmark summaries remain
unchanged under ignored `results/multi_fleet/recovery/`, `curves/` and
`warped_diagnostics/`. No replacement family or compatibility alias is provided.

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

`multi_fleet_recovery.R` retains seven isolated configurations (100 replicates
each) and seven full-assessment scenarios (30 each): **910 configured attempts**.
The historical results below are the supported scenarios from the original
1,690-attempt benchmark; they were not refitted or relabeled during retirement.
Seeds, generating models and optimizer controls are unchanged.

Fresh runs default to `results/multi_fleet/recovery_supported/`, keeping the
original benchmark immutable. Source-fingerprint guards remain active.
`postprocess` adds N uncertainty to private objectives; `metrics` reports
recovery, coverage and between-replicate Monte Carlo SEs. Removing scenarios
does not authorize overwriting or mixing the historical checkpoints.

Full assessments use 12 training years, ages 1–10, two fleets, known M, distinct
catch weights, and two held-out years. The SCAA has logistic/double-logistic
fleet curves; SAM scenarios use shared/separate age correlations and fleet F
blocks (one unfished young age). Smooth scenarios
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

SAM correlation recovery on the fitted atanh scale had bias 0.041 ± 0.041
(shared), 0.038 ± 0.041 (Gillnet) and −0.112 ± 0.066 (Trawl), with interval
coverage 100%, 96.7% and 93.3%. Smooth/AR1 separation was poor: smoothing
log-SDs often shrank toward −9 despite generating values around −1.6 to −1.9.
Very wide intervals could still cover the truth. In the deliberately weak
scenario, age/year AR1 logit-parameter biases were −0.718 ± 0.093 and
−0.537 ± 0.101; coverage was only 65%/70%. Good population recovery among
converged fits does not establish identification of its component processes.

Isolated fits had 100/100 successes for logistic, both double-logistic midpoint
orders and tp/cr smooths; factor-by/tensor smooths had 97/100. One hundred
successes has exact 95% interval 96.4–100%. Function coverage was 95.4% for
ordinary double logistic and 85.3% for reversed midpoints; tp/cr/factor-by/tensor
smooth coverage was 87.0%/93.9%/89.2%/92.0%. Local Gaussian/Wald intervals can
under-cover.

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

## Final validation and delivery

The final focused run passed **395 assertions**, without warnings or failures.
The clean package build passed `R CMD check --no-manual`: **OK**, **4,218
assertions**, zero failures/warnings/skips. Comparative assessments returned
**60 passes and the same ten pre-existing failures** listed in the core section;
the additional fleet-reference test passed. No canonical aggregate exports
were retained from targeted runs.

The SCAA/SAM/NCAM cheatsheet fits, two-fold multi-model hindcast, SSB figure,
companion HTML and tracked PDF were regenerated through `cheatsheet/examples.R`,
with its source/provenance guards intact. The SAM example had one failed
hindcast attempt; differing validation coverage remains flagged. Both PDF pages
were rendered and visually checked for clipping and legibility. The PDF's
relative `tinyAM_dashboard.html` link and the adjacent self-contained HTML were
verified. The reference now also shows the existing cohort-covariate term and
optional latent recruitment boundary convention.

Core, scoped, curve, smooth, herring and multi-model example dashboards rendered
with `open_file = FALSE`. Interactive browser inspection was blocked when the
computer-use tool could not verify the current Chrome URL; it was not bypassed.
Manual control/rotation/resize checks therefore remain unverified. Programmatic
widget, orientation, scope, interval and rendering checks passed.

Keep the default-start smooth estimation failures visible. This work adds
modeling capabilities, not a recommendation for routine use of every flexible
combination or a claim of improved/exact official assessment reproduction.
Generated fits, warnings, dashboards, numerical summaries and source hashes
remain in ignored `results/multi_fleet/` and `analysis/sam_equivalence/results/`.
Only source code, concise reports, documentation and the designated cheatsheet
PDF are tracked. Version is `0.1.1.9001`, targeting a future 0.2.0 release.

## Fisheries dashboard refinements

On `fleet-dashboard-refinements`, Fisheries replaces the separate Fishery and
Fleets entries. Fleet pages retain reported yield, catch-at-age, F-at-age and
available dotted means, plus all five observation-residual views. Aggregate
fishery follows the fleets and shows total yield; unlabeled models are never
allocated among fleets. With no fleet labels it retains the existing complete
single-fishery diagnostics. Incomplete reported fleet yields remain missing.

Translated reference count units now separate explicit inference annotations
from recognized count labels, including both herring annotation forms. The
original statements remain in notes and raw source records; values, SEs and
intervals are unchanged. Unknown/descriptive survey units stay unchanged and
incompatible known units still fail overlay checks. Herring's Population
patterns colorbar reads `Abundance<br>(thousand fish)`.

Focused checks pass **425 assertions** without failures or warnings: 34 new
Fisheries checks, 204 dashboard regressions and 187 unit/matrix/retrospective
checks. Standalone reference-unit tests also pass. The comparative suite returns
**61 passes and the same ten pre-existing failures** listed in the core section.
Its preview-render test initially lacked the installed Pandoc path; it passes
with that path set. Both attempts are retained. The temporary comparative
checkout and its generated readiness table were removed; canonical exports
were not updated.

Herring, three-model cheatsheet and single-SCAA dashboards rendered from
existing fits with `open_file = FALSE`. The single-fit check confirms unchanged
parameters and RTMB objective state. The cheatsheet generator reused its fit
and hindcast caches, regenerated its companion HTML and two-page PDF, and
preserved provenance checks. Both PDF pages and the relative companion link
were inspected. Interactive Chrome inspection succeeded for this refinement:
fleet order/selection, residual panels, F age selection and both orientations,
heatmap/surface switching, the concise unit label and resizing. This supersedes
the earlier browser-access limitation for those checks, without claiming an
exhaustive manual check of every widget combination.

The first full package run passed 4,251 assertions and failed one old test
expecting the retired Fishery heading. Its expectation was updated to Aggregate
fishery; all 204 dashboard regressions then passed. That check log is retained.
Curve, parameterization, starting values, optimizer controls, scientific
specifications and original recovery success rates are unchanged.

The fresh final build passes `R CMD check --no-manual`: **OK**, with **4,252
passing assertions**, zero failures, warnings or skips. Version remains
`0.1.1.9001`. Separate feature commits are `070ed90c` (diagnosis), `48c3ee2e`
(Fisheries navigation), `0d7227b2` (reference units), and `0b80a0e1` (updated
heading regression). This branch is intended for review, without a main merge.

## Retirement validation

On `retire-warped-normal`, six supported model fixtures retain their parameter
lists, objectives, gradients and predictions within a 1e-12 comparison tolerance.
The seven retained full generators, seeds and smooth workflow example also
match their recorded pre-cleanup baselines at that tolerance. No new curve,
initial-value policy or fitter change was introduced.

One smoke attempt per retained configuration produced seven isolated successes
and six full-assessment successes. The smooth-off assessment failed numerical
checks; its fit and warnings were retained without rescue tuning. All 14 runners
completed without uncaught study errors. These are workflow checks, not a new
recovery benchmark. Their separate output is
`results/multi_fleet/retirement_validation/smoke/`.
