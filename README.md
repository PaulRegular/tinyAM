
<!-- badges: start -->
[![R-CMD-check](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/PaulRegular/tinyAM/graph/badge.svg)](https://app.codecov.io/gh/PaulRegular/tinyAM)
<!-- badges: end -->

# tinyAM <img src="man/figures/logo.png" align="right" height="138" alt="" />

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

Age-year deviations in abundance, fishing mortality, and natural mortality can be represented using a common two-dimensional process model, including IID, highly correlated approximate random-walk, and estimated AR(1) structures.

Mean structures for quantities such as fishing mortality, natural mortality, catchability, and observation error can be specified using familiar R formulas and design matrices.

## Workflow

A typical tinyAM workflow is:

```r
library(tinyAM)

fit <- fit_tam(
  cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(
    process = "iid",
    init_N0 = FALSE
  ),
  F_settings = list(
    process = "approx_rw",
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

```r
fit_ar1 <- update(
  fit,
  F_settings = list(
    process = "ar1",
    mu_form = ~ F_a_block + F_y_block
  )
)
```

Because fitted models share a common structure, they can be compared using the same downstream tools.

```r
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

## Current status

tinyAM is under active development and should currently be considered **highly experimental**.

- It has **not been validated for operational stock assessment**.
- It is intended primarily for **testing model structures, assumptions, and ideas**.
- Model formulations and the user interface may change substantially.
- Breaking changes should be expected while the package develops.

The package currently uses Northern cod data as an example and development test case, but the modelling framework is not intended to be specific to that stock.

## Installation

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("PaulRegular/tinyAM")
```

## Why tinyAM?

Many useful assessment questions involve changing something relatively small:

> What happens if process error is assigned to mortality rather than abundance?

> Does allowing natural mortality to vary improve retrospective behaviour?

> How sensitive are estimates to the assumed correlation structure of fishing mortality?

> Does a model that fits the historical data well also predict withheld observations well?

tinyAM is intended to make questions like these relatively inexpensive to translate into fitted models.

The emphasis is therefore less on providing a single preferred assessment formulation and more on providing a transparent framework in which alternative assumptions can be expressed and tested.
