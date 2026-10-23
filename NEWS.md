# tinyAM 0.1.1.9000

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
