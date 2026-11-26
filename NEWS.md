# tinyAM 0.1.1.9001

* **Breaking change:** Use `R_settings` in place of `rec_settings` to configure
  recruitment, matching the model's notation. Update model calls and refit
  models saved with the old setting name.

* Reject explicitly labelled age compositions in catch and survey inputs.
  `check_tam_data()` warns when repeated annual sums suggest age proportions.

* Add `"logit-link"` (default) and `"poisson-link"` types for
  `delta_lognormal()`, `delta_gamma()`, and `delta_gengamma()`. Add `tweedie()` with
  an optional fixed power for catch and survey observations. Tweedie models
  require the suggested `tweedie` package for residuals and positive medians.

* Replace the previous Poisson-link-only constructors with the `delta_*()` interface. Shared
  observation-family help describes prediction scales, dispersion, and zeros.
  The Northern cod guide compares four families in a common assessment.

* Pool complete native-age observations for delta and Tweedie plus groups.
  Keep incomplete tail observations missing and preserve native counts for
  observed-catch hindcasts.

* Default to quantile residuals from a conditional Gaussian draw in state-space
  fits. Keep fitted-state quantiles available through `method = "mle-eb"` and
  standardized residuals for lognormal fits without random effects.

* Add Gamma and generalized gamma observations. Use `stats::Gamma(link = "log")`
  or `gengamma()` in catch and index settings.

* Calculate conditional AIC for all observation families. AIC and cAIC can
  compare families that use identical observations and units, including sampled
  zeros where applicable.

* Delta families fit one fixed or estimated soft detection scale per fleet or
  survey. The overall arithmetic mean uses a log link.
  Outputs include fitted zero probabilities and randomized quantile residuals.

* Rename the catchability marker `power()` to `q_power()` to avoid masking
  `stats::power()`. Use `q_power(N)` in catchability formulas.

* Add symmetric dome-shaped `normal()` curves for fishing mortality and
  survey catchability. Document all four selectivity markers together.

* Allow positive within-year age AR1 observation correlation separately for
  catch fleets and surveys with `cor_struct`. Independent errors remain the default.

* Allow penalized mgcv `s()` and age-year `t2()` terms in catch and index
  observation SD formulas. Estimate smoothness jointly with the assessment.

* Add concise website guides for comparing Northern cod model structures and
  fitting multiple fisheries and surveys. Link both from a Guides menu alongside
  the model-structure examples and cheatsheet. Include an evaluated NCAM
  retrospective SSB comparison and a companion simulated-cod dashboard.

* Add simulated, 3Ps-inspired cod inputs with three active fleets and two
  surveys of coastal and offshore areas. Keep compact generating truth and
  monitoring provenance separate from the estimated observations; SimSurvey
  is required only by the data-generation script. Use age-based selectivity
  and annual fleet F targets with a trawl reduction after year 10. Compare
  fleet-specific and aggregate assessments with an age-14 plus group.
  Show aggregate catch fits and residuals alongside fleet pages in mixed
  dashboards.

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

* Run comparative assessments from committed data in a local SAIO checkout
  and record its commit with the results. Separate fitting, dashboard
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

* `q_power(N, by = group)` in `q_form` estimates positive abundance exponents
  using abundance at survey time. Missing group labels retain p = 1 without
  removing observations. Effective and baseline catchability are reported separately.
  A focused Northern Shelf haddock comparison reproduces the pinned 2026 SAM
  likelihood and fit; see `analysis/sam_equivalence` for checks and limitations.

* `R_settings$init = "latent"` integrates unpenalized recruitment boundary
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

* `R_settings$process = "off"` supports deterministic recruitment formulas
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
