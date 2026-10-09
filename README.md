
<!-- badges: start -->

[![R-CMD-check](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/PaulRegular/tinyAM/graph/badge.svg)](https://app.codecov.io/gh/PaulRegular/tinyAM)

<!-- badges: end -->

# tinyAM <img src="man/figures/logo.png" align="right" height="138"/>

**tinyAM** (*tiny assessment model*; aka TAM) is an experimental R package for fitting age-structured stock assessment models using [RTMB](https://github.com/kaskr/RTMB).

The package is designed to be **flexible, modular, and lightweight**, with an emphasis on making alternative model structures and process assumptions easy to specify, inspect, simulate, and compare.

Rather than providing a collection of predefined assessment models, tinyAM uses a common age-structured state-space framework whose components can be modified to represent different hypotheses about population and fishery dynamics.

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

A random walk has no penalty on its starting level. Its SD describes annual increments; an AR1 SD is an innovation scale, not the marginal variability of the states. Recruitment also follows a temporal random walk.

Mean structures for quantities such as fishing mortality, natural mortality, catchability, and observation error can be specified using familiar R formulas and design matrices.

Initial abundance is specified independently of the subsequent N process through `N_settings$init`: `"exp"` (default) uses parsimonious survivorship from fixed first-year recruitment (`log_r0`), `"free"` estimates fixed older-age `log_n0` states, and `"random"` estimates random states with IID survivorship residuals (`eta_log_n0`) and a separate SD. Recruitment and N process states begin in year 2. Random initialization requires at least two ages and warns below ten ages because its SD may be weakly identified.

For an active M process, the default normally shares one state across all ages except the youngest, starting in year 2 (see `?prepare_tam` for very short age ranges). Earlier M and excluded ages retain their supplied or mean values. F is estimated in every historical year. These boundary assumptions can affect initial abundance and catchability estimates; they are not guarantees of identifiability. The full equations and conventions are in `help("tinyAM-model", package = "tinyAM")`.

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
