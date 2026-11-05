# tinyAM 0.1.1.9001

* Add concise website guides for comparing Northern cod model structures and
  fitting multiple fisheries and surveys. Link both from a Guides menu alongside
  the model-structure examples and cheatsheet. Include an evaluated NCAM
  retrospective SSB comparison and a companion simulated-cod dashboard.

* Add simulated, 3Ps-inspired cod inputs with three active fleets and two
  surveys of coastal and offshore areas. Keep compact generating truth and
  monitoring provenance separate from the estimated observations; SimSurvey
  is required only by the data-generation script.
  Use a compact, static-availability teaching scenario with known mean length
  and a reference M baseline; compare fixed and AR1 M with exact cohort survival.

* Organize plots by purpose and colocate population getters. Share fitted test
  fixtures, remove redundant renders and unused implementations, and retire the
  `check_obs()` and `tidy_mat()` aliases in favor of `check_tam_data()` and
  `tidy_array()`. Old saved-fit migrations are no longer maintained.

* Open dashboards on Background when supplied, followed by Model structure
  with Diagram and Function calls tabs. Use compact, responsive HTML boxes,
  readable typography and cleaner SVG connectors without a new dependency.

* Restore catch-only preparation, projections, simulation and tidy outputs,
  retaining empty survey tables without unused catchability or survey-SD
  parameters. Render reporting references with missing yield, observations or
  age-year surfaces without inventing unavailable quantities.

* Add `plot_model()` and a Model structure dashboard page showing biological
  connections, process switches, formulas, initialization and input predictors.
  Compare stored specifications with an independent model selector; link
  components to maintained equations without calculating or changing fits.

* Offer Standardized / OSA residual controls without bubble outlines. Start with
  retained OSAs when available; overlays use a common residual type and model
  frames choose separately. Preserve partial gaps and identify the OSA method;
  plotting never calculates residuals. Explicit row-wise hover values stay
  aligned through model frames and missing residuals.

* Select any supported OSA algorithm through `method` and retain diagnostics
  with `add_osa_residuals()`, without refitting. Explicit computation requests
  refresh catch and survey residuals together; cached extraction stays cheap.
  Retain method/settings/coverage metadata through tidying and model comparisons.

* Add `get_ssb()`, `get_recruitment()`, `get_fbar()`, `get_mbar()`, `get_f()`
  and `get_m()` for common population tables, retaining uncertainty, units,
  projections and comparison labels without calculating diagnostics.

* Print a concise input-validation summary with `check_tam_data()`; internal
  preparation remains quiet. Install dashboard R dependencies automatically and
  check for Pandoc before rendering.

* Preserve smooth warm-start functions across equivalent mgcv basis signs and
  rotations, while rejecting changes to predictor grids or Gaussian penalties.
  This fixes platform-dependent retrospective warm-start checks.

* Run comparative assessments from an explicitly refreshed local SAIO snapshot,
  recording its commit in a tracked manifest. Separate fitting, dashboard
  rendering and overview figures; honor the caller's future plan. Keep curated
  and preliminary membership and revision-specific result caches separate.

* Group fishery diagnostics under Fisheries: one complete page per fleet,
  followed by Aggregate fishery. Unlabeled catches remain stock-wide; aggregate
  residuals are not invented for multi-fleet fits.
* Separate inferred count-unit labels from source notes in translated reporting
  references, keeping herring abundance labels concise without changing values.
* Add mgcv thin-plate/cubic `s()` and age-year `t2()` smooths to F means and
  catchability formulas, including fleet/survey scopes. Estimate normalized
  Gaussian smooth penalties jointly, retain training bases for prediction,
  and report complete smooth functions with joint uncertainty. Flag weak
  separation from IID/AR1 F residuals; intrinsic F random walks are excluded.
* Retire experimental `warped_normal()` support before release. Its code and
  diagnostic history remain at revision `417e2636`; historical caches are retained.
* Add WHAM-style `double_logistic()` curves for F and q. Report maximum realized
  fleet F without normalizing selectivity peaks.
  Keep likelihood calculations stable when selectivity is extremely small.
* Add `for_fleet()` and `for_survey()` for group-specific F means and
  catchability relationships within one formula, alongside shared effects.
  Scoped curves, Gaussian effects and abundance powers retain their own groups.
* Add optional catch fleets with independent F surfaces on one shared population,
  shared or fleet-specific correlation groups and F age blocks, and typed
  fleet-specific catch weights. Total F drives population survival.
* F mean, process-SD and correlation covariates are supplied in the catch table;
  recruitment, N and M covariates remain in the population table.
* Expand reproducible curve/smooth recovery and multi-fleet assessment trials.
  Default-start failures and weak smoothing identification remain
  documented limitations; these flexible combinations are experimental.

# tinyAM 0.1.1.9000

* Add optional F age blocks with explicit closing boundaries. Blocks share
  actual F; younger ages below the first boundary have exactly zero F.
  M accepts the same boundaries while preserving all legacy specifications.
  Pinned SAMjr examples and native SAM provide numerical equivalence checks.

* Partial OSA requests now map returned residuals to their original observations;
  unrequested and conditioning rows retain missing residuals.

* `power(N, by = group)` in `q_form` estimates positive abundance exponents
  using abundance at survey time. Missing group labels retain p = 1 without
  removing observations. Effective and baseline catchability are reported separately.
  A focused Northern Shelf haddock comparison reproduces the pinned 2026 SAM
  likelihood and fit; see `analysis/sam_equivalence` for checks and limitations.

* `rec_settings$init = "latent"` integrates unpenalized recruitment boundary
  states without an initial density or SD. The default `"free"` preserves
  fixed boundary parameters, independently of the subsequent process and N initialization.

* Optional `weight$type` distinguishes complete stock and catch weights.
  Stock weights determine biomass and SSB; catch weights determine yield,
  catch targets and observed-catch hindcasts. Untyped and stock-only inputs
  retain their existing behavior. `W_catch` reports effective catch weights.

* `vis_tam()` now accepts retrospective and hindcast results directly. Select
  one model's reference and successful peels; omitted attempts are explained
  on Background. Forecast hover text identifies each peel's origin and horizon.

* The assessment translation batch honours its dashboard switch, reports stock
  progress, and saves per-stock checkpoints immediately. Its script applies a
  configurable subprocess deadline, recording timeouts without changing model
  assumptions or blocking the remaining assessments.

* Add reproducible Gaussian joint-sample process residuals to `tidy_tam()`,
  alongside fitted deviations, with explicit uncertainty status and a local seed.
  New standalone process plots and the Random dashboard separate Recruitment,
  N, F, M, fitted deviations and joint-sample residuals; signed matrices retain
  unique M state blocks. Model likelihoods and fitting assumptions are unchanged.
  Each process plotting helper returns one Plotly widget selected by `component`;
  flexdashboard supplies the Recruitment/N/F/M dropdowns. Residual bubbles use
  N/F/M buttons and a model animation slider, matching ordinary residual plots.

* `check_tam()` accepts multiple fitted models or a named `model_list`, with
  labeled CLI output and separate diagnostics for each model. Single-model
  results and all diagnostic calculations are unchanged. Supply `grad_tol`
  and `detailed` by name; unnamed arguments now represent additional models.

## Example data

* Simplified `cod_obs` by removing the optional `smith_sound_effect` and
  `log_capelin_per_cod` predictors. Existing observations, mortality blocks,
  units, and the allometric natural-mortality baseline are retained.

* The cheatsheet now demonstrates the case-study SCAA, SAM, and NCAM structures,
  simple multi-model visualization and hindcasting, and tidy outputs for static
  report figures. Its cached SSB comparison uses the package palette.

## Evaluation refinements

* Deterministic recruitment with an annual `factor(year)` effect now carries its
  terminal historical effect into projections and hindcasts, without estimating
  unsupported future-year coefficients. Numeric time trends are unchanged.

* `residuals(fit, type = "osa")` calculates missing OSA residuals on request,
  reusing cached results when present. Calculations use a temporary fitted
  objective without refitting or changing the supplied fit.

* Hindcasts expose their one-year, observed-catch default in Usage. They accept
  actual removals or one F multiplier; numeric catch scenarios and horizon vectors
  remain available only in ordinary projections. Printed statistics and forecast
  conditioning assumptions are documented explicitly.

* `rec_settings$process = "off"` supports deterministic recruitment formulas
  and BH/Ricker curves, retaining free recruitment boundary states and derived
  uncertainty without extra annual states or process parameters.
  Joint uncertainty simulation falls back to fixed-parameter draws when no
  random effects remain. Recruitment warm starts respect changes in fitted bases.

* F mean formulas are trained on historical support and carry terminal design
  rows into projections, fixing SCAA hindcasts with `factor(year)`. Cohort
  means share the mean-formula compiler while retaining destination-age checks.

## Evaluation workflows

* `fit_retro()` and `fit_hindcast()` accept multiple fits and return `tam_retro`
  objects with a concise print method, detailed fold diagnostics and retained failures.
* Hindcasts support multiple forecast horizons, biological averaging settings,
  persistence-relative MASE, survey age/aggregate scores and catch-composition scores.
* **Breaking:** forecast controls now belong in `proj_settings`. Fitted folds are nested by model and terminal year.
  Scores are in `$scores`, replacing the former single `$hindcast_rmse` value.
* `start_par` accepts a fitted model, including when uncertainty estimation failed.

## Diagnostic displays

* `residuals()` extracts keyed observation residuals without refitting.
  OSA calculation is requested explicitly with `type = "osa"`.
* Random parameters distinguish latent N/F/M states from their deviations.
  Recruitment deviations retain signed annual trends; N/F/M use bubbles.
  The same point tables are available in `tidy_tam()$deviations`.

## Cohort covariates

* `N_settings$mu_form` accepts ordinary fixed covariates for proportional
  apparent cohort gains/losses after survival, indexed by destination year/age.
  Recruitment and initial abundance remain separate. Constant directions under
  an intrinsic N random walk are rejected rather than silently estimated.

## Formula-driven F

* `F_settings$process = "off"` uses an explicit ordinary, monotonic or logistic
  mean without F latent states or residual variance/correlation parameters.
  Existing terminal-F and catch-conditioned projection controls remain available.
