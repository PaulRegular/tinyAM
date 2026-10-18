
<!-- badges: start -->

[![R-CMD-check](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/PaulRegular/tinyAM/graph/badge.svg)](https://app.codecov.io/gh/PaulRegular/tinyAM)

<!-- badges: end -->

# tinyAM <img src="man/figures/logo.png" align="right" height="138"/>

**tinyAM** (*tiny assessment model*; aka TAM) is an experimental R package for fitting age-structured stock assessment models using [RTMB](https://github.com/kaskr/RTMB).

The package is designed to be **flexible, modular, and lightweight**, with an emphasis on making alternative model structures and process assumptions easy to specify, inspect, simulate, and compare.

Rather than providing a collection of predefined assessment models, tinyAM uses a common age-structured state-space framework whose components can be modified to represent different hypotheses about population and fishery dynamics.

## Documentation and examples

The [package website](https://paulregular.github.io/tinyAM/) provides the function
reference, [cheatsheet](https://paulregular.github.io/tinyAM/cheatsheet/tinyAM_cheatsheet.pdf)
and [example dashboard](https://paulregular.github.io/tinyAM/cheatsheet/).
The website is published from `main` by GitHub Actions.

To reproduce the cheatsheet locally:

```sh
Rscript cheatsheet/examples.R prepare
Rscript cheatsheet/examples.R render
```

The first command reuses valid fits and hindcasts, then writes a combined
main-fit dashboard and individual hindcast dashboards to `cheatsheet/output/`.
The individual dashboards share one local `output/assets/` directory; keep it
with their HTML files. The combined dashboard is self-contained.
The second command renders the PDF without fitting. Current results live in ignored
`cheatsheet/cache/`; temporary rendering files are removed automatically.
Use `Rscript cheatsheet/examples.R dashboard` to build just the combined dashboard.

To preview the package site, run
`withr::with_envvar(c(CI = "true"), pkgdown::build_site(preview = FALSE))`.
This skips pkgdown's optional external favicon service.
The **Package website** workflow builds the PDF link and combined dashboard
on pushes to `main`, with a manual rebuild option under Actions. GitHub Pages
uses **GitHub Actions** as its publishing source; the `github-pages` environment
permits deployment from `main`. The PDF itself is the reviewed tracked copy;
website builds do not rerun hindcasts or compile LaTeX.
The workflow deploys generated artifacts; site files and HTML are not committed.

Historical cheatsheet trials are recorded in Git, including archival commit
`9920c9a`; retired experiments are absent from the current checkout.

## Philosophy

tinyAM is intended primarily as a **model-development and experimentation tool**.

The general philosophy is that assessment models can be viewed as hypotheses about where variation and mismatch between population dynamics and observations arise. For example, unexplained variation might be represented through:

- abundance or cohort process error;
- fishing mortality;
- natural mortality;
- recruitment;
- survey catchability; or
- observation error.

tinyAM makes many of these assumptions configurable while retaining a common model and output structure. This makes it possible to change one structural assumption at a time and evaluate its consequences using the same fitting, diagnostic, simulation, retrospective, and hindcast tools.

The goal is not to make the underlying models statistically trivial. Instead, **“tiny” refers to keeping the modelling vocabulary small**: model complexity should arise from combining a limited number of understandable components rather than from maintaining many special-purpose model implementations.

## Model structure

The core model is an age-structured state-space model with:

- recruitment varying through time;
- cohort survival governed by total mortality;
- a terminal plus group;
- fishing mortality at age;
- fixed or time-varying natural mortality;
- Baranov catch-at-age predictions;
- one or more survey indices with within-year sampling times; and
- lognormal observation models for catch and survey indices.

Departures from expected log abundance or mean log mortality can follow three processes:

- `"iid"`: independent departures in each year and age or age block;
- `"rw"`: departures accumulate from year to year, independently across ages;
- `"ar1"`: departures are correlated across adjacent years and ages and return toward a mean.

F also supports `"cor_rw"`: annual log-F changes are correlated between ages,
while their effects accumulate over time. N, F and M settings accept ordinary
age-based `sd_form` formulas, such as `~ age_group`, to share process SDs.
These formulas use columns in `data$pop_covariates` and do not share latent states.
Their default `~ 1` retains a common SD.

`data$maturity$spawn_time <- 0.25` accounts for survival to spawning at
one quarter of the year. The same SSB supplies stock–recruit relationships.
Omit the column for start-of-year SSB. Timing can vary by year and age;
projections retain terminal values at each age. See `?"tinyAM-model"` for equations,
SD interpretations and examples.

Optional `unit` columns make plot labels explicit without changing input values:
use `"thousand fish"` for catch and comparable indices, and `"kg"` for
weight. This combination produces biomass and SSB in tonnes. Survey units may
differ by survey; native labels such as `"fish per tow"` are retained. Keep
units consistent within each series and call `check_tam_data()` after editing.
The included `cod_obs` carries these labels. See `?as_tam_data` for details.

Settings lists can be partial. Initial fits inherit package defaults; updates
retain the original settings. For example,
`update(fit, N_settings = list(init = "free"))` changes only initialization.

A random walk has no penalty on its starting level. Its SD describes annual increments; an AR1 SD is an innovation scale, not the marginal variability of the states. Recruitment defaults to a temporal random walk; `rec_settings` separates its process from mean formulas for covariates and Beverton–Holt or Ricker stock–recruit relationships.

Mean structures for quantities such as fishing mortality, natural mortality, catchability, and observation error can be specified using familiar R formulas and design matrices.

Initial abundance is specified independently of the subsequent N process through `N_settings$init`: `"exp"` (default) uses parsimonious survivorship from fixed first-year recruitment (`log_r0`), `"free"` estimates fixed older-age `log_n0` states, `"latent"` integrates older initial states without a density, `"mono"` constrains abundance to decline from r0, and `"capped"` bounds older ages by r0 without age ordering. All retain fixed r0; both constraints include the plus group and allow ties. N process states begin in year 2. Recruitment starts there by default; stock–recruit models retain additional free early recruitment states when parent SSB precedes the modeled period. All initializers require at least two ages. Latent initialization is closer to SAM’s usual older-age boundary, but r0 remains fixed. Process simulations condition on initial states; joint uncertainty draws can vary latent initial abundance.

For an active M process, the default normally shares one state across all ages except the youngest, starting in year 2 (see `?prepare_tam` for very short age ranges). Earlier M and excluded ages retain their supplied or mean values. F is estimated in every historical year. These boundary assumptions can affect initial abundance and catchability estimates; they are not guarantees of identifiability. The full equations and conventions are in `help("tinyAM-model", package = "tinyAM")`.

## Input data

Use a named list of `catch`, `index`, `weight`, and `maturity` tables, with an
optional `pop_covariates` table. Population covariates and supplied M assumptions
belong in `pop_covariates`; observation-SD and q covariates stay in catch/index.
Weight, maturity and spawning timing stay in the biological tables.

```r
data <- as_tam_data(cod_obs)
data                              # compact previews of all tables
check_tam_data(data)
```

The class is optional; plain lists also work. Population covariates have one
row per input year and age. Missing values are allowed where unused, and
projections carry terminal values forward by age. See `?as_tam_data`.

## Recruitment formulas

Use `rec_settings` to describe recruitment at the youngest modeled age:

```r
rec_settings = list(process = "rw", mu_form = NULL) # default
rec_settings = list(process = "iid", mu_form = ~ 1)
rec_settings = list(process = "iid", mu_form = ~ bh(ssb))
rec_settings = list(process = "ar1", mu_form = ~ ricker(ssb))
rec_settings = list(process = "rw", mu_form = ~ temperature)
```

`N_settings$process` separately controls survival deviations in older cohorts.
Annual covariates belong in `data$pop_covariates` at the youngest modeled age; older-age entries
may be missing. IID/AR1 require an explicit mean, such as `mu_form = ~ 1`.
Supply `sd` or AR1 `phi` in the settings list to fix them; otherwise they are
estimated. Curves describe median recruitment from parent SSB at the supplied
spawning time (zero by default). The lag defaults to recruitment age, and early years without modeled
parent SSB remain free fixed states. Curve-plus-RW residuals are not supported.

Use `check_tam()` to inspect uncertainty as well as numerical convergence.
Limited SSB variation can leave curve shape poorly estimated even in a
converged fit. `tidy_recruitment()` provides parent pairs, predictions and
residuals; `vis_tam()` adds a stock–recruitment page when a curve is fitted.
See `?recruitment_formulas` for equations and safeguards, and
[runnable examples](https://github.com/PaulRegular/tinyAM/blob/HEAD/inst/examples/example_recruitment.R).

## Rising mean F curves

When fisheries progressively select older fish, the long-term mean F can be
specified with a flexible rising curve or a smooth logistic curve:

```r
F_settings = list(process = "iid", mu_form = ~ mono(age))
F_settings = list(process = "ar1", mu_form = ~ logistic(age))
```

Ordinary terms supply baseline F. `mono()` fits non-negative log-F increments
`dmu_f`, including exact plateaus; `logistic()` multiplies baseline F by rising
selectivity. Annual realized F can depart from these means. RW residuals do
not identify static age curves, so these combinations are rejected. See
`?catchability_curves` for equations, grouping and uncertainty limitations, and
[fitting examples](https://github.com/PaulRegular/tinyAM/blob/HEAD/inst/examples/example_f_mean_curves.R).

## Monotonic survey catchability

Catchability formulas can retain an unconstrained shape or impose monotonic changes between ordered blocks:

``` r
~ q_block                                      # unconstrained
~ mono(q_block)                                # non-decreasing
~ survey + mono(q_block, by = survey)           # independent survey curves
```

Use these as `index_settings$q_form`. Numeric levels are sorted increasingly; factor levels follow their declared order. The first represented level is the baseline. Later levels add non-negative log-q increments `dq`, fitted directly with a zero lower bound and initialized at 0.05. Zero increments allow exact plateaus between separate levels; pooled blocks also give exact plateaus. Estimates and SEs for `dq` share the increment scale; Wald inference is only a local approximation at an active boundary. Each `by` group uses its own represented levels (at least two) and independent steps. Ordinary terms supply baselines: omitting `survey` in the last example shares one intercept across surveys. Other ordinary covariates are held constant when interpreting monotonicity.

## Workflow

A typical tinyAM workflow is:

``` r
library(tinyAM)

fit <- fit_tam(
  cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(
    process = "iid",
    init = "exp"
  ),
  F_settings = list(
    process = "rw",
    mu_form = NULL
  ),
  M_settings = list(
    process = "off",
    mu_supplied = ~ I(0.3)
  ),
  catch_settings = list(
    sd_form = ~ 1
  ),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ q_block
  )
)

fit
```

Alternative models can then be constructed by modifying the fitted call:

``` r
fit_ar1 <- update(
  fit,
  F_settings = list(
    process = "ar1",
    mu_form = ~ F_a_block + F_y_block
  )
)
```

Because fitted models share a common structure, they can be compared using the same downstream tools.

``` r
fits <- list(
  "F random walk" = fit,
  "F AR1" = fit_ar1
)

tabs <- tidy_tam(model_list = fits)
vis_tam(fits)
```

## Diagnostics and simulation

tinyAM includes tools for examining both model fit and model behaviour, including:

- standardized and one-step-ahead residual diagnostics;
- observed-versus-predicted comparisons;
- age-year process and population visualizations;
- retrospective analyses and Mohn's rho;
- one-step-ahead hindcasts and hindcast RMSE;
- short-term projections;
- simulation from fitted models; and
- propagation of fixed-effect or joint fixed/random-effect uncertainty.

The likelihood also serves as the model's simulation engine, helping to keep estimation and simulation assumptions consistent.

``` r
check_tam(fit)
fit$opt$message
head(fit$pop$ssb)
plot_trend(fit$pop$ssb, ylab = "Spawning stock biomass")
plot_age_year(fit$pop$F, ylab = "F")
plot_age_year(fit$obs_pred$catch, plot = plot_obs_pred, ylab = "Catch")

# One-year forecasts refitted from historical terminal years
# hindcasts <- fit_hindcast(fit, folds = 3)

# Possible observations conditional on the fitted population history
sims <- sim_tam(fit, n = 10, par_uncertainty = "none",
                redraw_random = FALSE, seed = 1)
```

`redraw_random = TRUE` regenerates process states across the whole modeled period; it is not a future-only forecast conditional on the fitted history. Projections use terminal historical F times a chosen multiplier, recent mean weights and maturities, and continuing recruitment/N/M processes. Compare diagnostics and sensitivities before treating differences among model outputs as biological evidence.

## Outputs

Internally, tinyAM uses matrices and arrays where these naturally reflect the age-year structure of the model.

For analysis, fitted quantities are converted to tidy data frames. A `tam_fit` object contains the underlying RTMB model and optimization objects as well as convenient summaries of:

- fixed and random parameters;
- population abundance and recruitment;
- fishing and natural mortality;
- biomass and spawning stock biomass;
- observations and predictions;
- residual diagnostics; and
- uncertainty estimates.

The underlying model components remain accessible so that unusual or experimental analyses are not hidden behind a high-level interface.

Uncertainty tables put `est`, `lwr`, and `upr` first, then `se` and `se_scale`. Estimates and confidence limits are on the reported scale; SEs stay on the estimation scale. A small log-scale SE, such as 0.10, approximates 10% relative uncertainty (CV). This approximation is poor for large SEs, so use the confidence limits to communicate uncertainty. Observation `sd` describes variation in log measurements and is a different quantity from the SE of an estimate.

Catch and survey predictions are conditional medians of lognormal observations. Their arithmetic means are higher. Zero observations are treated as missing, not as observations from a count or censoring model. SSB is calculated at the start of the year without an extra spawning-time mortality correction.

For practical help, start with `?fit_tam`, `?prepare_tam`, and `?mono`. For mathematical details, use `help("tinyAM-model", package = "tinyAM")`, `?dprocess_ar1`, and `?dprocess_rw`.

## Current status

tinyAM is under active development and should currently be considered **highly experimental**.

- It has **not been validated for operational stock assessment**.
- It is intended primarily for **testing model structures, assumptions, and ideas**.
- Model formulations and the user interface may change substantially.
- Breaking changes should be expected while the package develops.

The package currently uses Northern cod data as an example and development test case, but the modelling framework is not intended to be specific to that stock.

## Installation

The development version can be installed from GitHub:

``` r
# install.packages("pak")
pak::pak("PaulRegular/tinyAM")
```

## Why tinyAM?

Many useful assessment questions involve changing something relatively small:

> What happens if process error is assigned to mortality rather than abundance?

> Does allowing natural mortality to vary improve retrospective behaviour?

> How sensitive are estimates to the assumed correlation structure of fishing mortality?

> Does a model that fits the historical data well also predict withheld observations well?

tinyAM is intended to make questions like these relatively inexpensive to translate into fitted models.

The emphasis is therefore less on providing a single preferred assessment formulation and more on providing a transparent framework in which alternative assumptions can be expressed and tested.

## Next priorities

- Review the latest Pacific assessment examples: diagnostics, comparison
  definitions, scientific approximations and remaining source blockers.
- Plan and implement conditional AIC separately, beginning with its statistical
  specification and validation. No conditional-AIC work is included in the
  repository cleanup or website publication.
