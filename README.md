<!-- badges: start -->

[![R-CMD-check](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaulRegular/tinyAM/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/PaulRegular/tinyAM/graph/badge.svg)](https://app.codecov.io/gh/PaulRegular/tinyAM)

<!-- badges: end -->

# tinyAM <img src="man/figures/logo.png" align="right" height="138"/>

**tinyAM** (*tiny assessment model*; aka TAM) is an R package for fitting age-structured stock assessment models using [RTMB](https://github.com/kaskr/RTMB).

The package is designed to be **flexible, modular, and lightweight**, with an emphasis on making alternative model structures and process assumptions easy to specify, inspect, simulate, and compare.

The development version can be installed from GitHub:

```r
# install.packages("pak")
pak::pak(c("PaulRegular/tinyAM", "rmarkdown", "flexdashboard"))
```

See the [package website](https://paulregular.github.io/tinyAM/) and
[cheatsheet](https://paulregular.github.io/tinyAM/cheatsheet/tinyAM_cheatsheet.pdf)
for guidance and examples.

## Philosophy

Assessment models express hypotheses about where unexplained variation arises:
recruitment, abundance deviations, fishing mortality, natural mortality, catchability
or observations. tinyAM lets you change one or more assumptions and examine
their consequences using the same diagnostic, simulation and forecasting tools.

“Tiny” refers to a small modelling vocabulary: complexity comes from combining
understandable components, rather than choosing among many special-purpose
models.

## 1. Set up tidy inputs

The required tables are `catch`, `index`, `weight` and `maturity`;
`pop_covariates` is optional. The Northern cod example shows all five:

```r
library(tinyAM)
print(cod_obs, n = 2)
```

```text

── tinyAM data ───────────────────────────────────────────────────────────────────────────

── catch: 1512 rows, 4 columns ──

Years 1954–2025; ages 0–20
     year age obs          unit
1    1954   0  NA thousand fish
2    1955   0  NA thousand fish
⋮       ⋮   ⋮   ⋮             ⋮
1511 2024  20   0 thousand fish
1512 2025  20   0 thousand fish

── index: 1512 rows, 7 columns ──

Years 1954–2025; ages 0–20
     year age obs     survey samp_time q_block          unit
1    1954   0  NA Fall 2J3KL       0.8       0 thousand fish
2    1954   1  NA Fall 2J3KL       0.8       1 thousand fish
⋮       ⋮   ⋮   ⋮          ⋮         ⋮       ⋮             ⋮
1511 2025  18   0 Fall 2J3KL       0.8    6-20 thousand fish
1512 2025  19   0 Fall 2J3KL       0.8    6-20 thousand fish

── weight: 1512 rows, 4 columns ──

Years 1954–2025; ages 0–20
     year age          obs unit
1    1954   0 2.965967e-04   kg
24   1955   0 2.965967e-04   kg
⋮       ⋮   ⋮            ⋮    ⋮
1631 2024  20 1.409146e+01   kg
1654 2025  20 1.369174e+01   kg

── maturity: 1512 rows, 3 columns ──

Years 1954–2025; ages 0–20
     year age obs
1    1954   0   0
2    1955   0   0
⋮       ⋮   ⋮   ⋮
1711 2024  20   1
1712 2025  20   1

── pop_covariates: 1512 rows, 6 columns ──

Years 1954–2025; ages 0–20
     year age F_y_block F_a_block allometric_M_assumption M_a_block
1    1954   0 1954-1991       0-1               4.9990199       0-1
2    1954  16 1954-1991     16-17               0.1793953      8-20
⋮       ⋮   ⋮         ⋮         ⋮                       ⋮         ⋮
1511 2025  10 2023-2025     10-11               0.2927334      8-20
1512 2025  20 2023-2025     18-20               0.1889190      8-20
```

The four required tables have numeric `year`, `age` and `obs` columns.
`obs` contains catch numbers, survey indices, individual weights or proportions
mature. Index rows also identify `survey` and `samp_time`: `0.8` means sampling
eight-tenths of the way through the year. `q_block` groups ages sharing
catchability in the example.

Keep q and observation-SD covariates in catch/index tables. The
`pop_covariates` holds population predictors, process-SD groups and supplied M
assumptions; it needs a complete year–age grid, without an `obs` column.
Keep weights in `weight`, and maturity and optional `spawn_time` in `maturity`.

```r
dat <- cod_obs                    # replace with your own named list of tables
check_tam_data(dat)               # validate before fitting
# dat <- as_tam_data(dat)         # optional class for compact printing
```

Catch in **thousand fish** and weight in **kg** give biomass and SSB in **tonnes**.
Optional `unit` columns label values; they do not convert them. Survey units
can differ. Catch/index `NA` values and zeros are excluded from the lognormal
likelihood; weights and maturities must be complete. No missing-observation
imputation or zero-count model is provided.

See [?check_tam_data](https://paulregular.github.io/tinyAM/reference/check_tam_data.html)
for input requirements and
[`?cod_obs`](https://paulregular.github.io/tinyAM/reference/cod_obs.html) for
details on the sources and preparation of the supplied Northern cod data.

## 2. Specify competing models

These are illustrative hypotheses, not package defaults or recommended stock
settings. The baseline lets recruitment and fishing mortality vary while
holding natural mortality at 0.3 per year:

```r
base <- fit_tam(
  dat, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "rw"),
  F_settings = list(process = "ar1", mu_form = ~ 1),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(sd_form = ~ 1, q_form = ~ mono(q_block))
)

N_iid <- update(base, N_settings = list(process = "iid"))
M_ar1 <- update(base, M_settings = list(
  process = "ar1", age_breaks = c(2, 4, 6, 8, 14)
))
```

Partial updates retain omitted settings. Here the alternatives add independent
cohort-survival variation or temporally correlated M shared within age blocks.
The age range `ages = 2:14` uses age 14 as the terminal plus group.

Other possibilities include:

| Component | Example options | Guidance |
|---|---|---|
| Recruitment | `rec_settings`: `"rw"`; `"iid"` with `mu_form = ~ bh(ssb)`; `"off"` with `~ factor(year)` | [Recruitment](https://paulregular.github.io/tinyAM/reference/recruitment_formulas.html): Beverton–Holt/Ricker curves support IID/AR1 residuals or deterministic recruitment, not RW residuals. |
| Cohort survival | `N_settings`: `"off"`, `"iid"`, `"rw"`, `"ar1"`; `mu_form = ~ temperature` | [Population model](https://paulregular.github.io/tinyAM/reference/tinyAM-model.html): apparent gains/losses cannot separate migration from mortality. |
| Fishing mortality | `F_settings`: `"off"`, `"iid"`, `"rw"`, `"ar1"`, `"cor_rw"`; `sd_form = ~ age_group` | [Settings](https://paulregular.github.io/tinyAM/reference/prepare_tam.html): deterministic F requires an explicit mean; `"cor_rw"` correlates annual changes across ages. |
| Natural mortality | `M_settings`: `"off"`, `"iid"`, `"rw"`, `"ar1"`; `sd_form = ~ age_group` | [Settings](https://paulregular.github.io/tinyAM/reference/prepare_tam.html): M can be supplied or estimated; age blocks share actual M. |
| Mean F and catchability | `~ age_group`, `~ logistic(age)`, `~ mono(age)`; q also supports `iid()`, `rw()`, `ar1()` effects | [Curves](https://paulregular.github.io/tinyAM/reference/catchability_curves.html) and [formula effects](https://paulregular.github.io/tinyAM/reference/formula_effects.html): supported combinations have identifiability safeguards; `q_link = "logit"` bounds q between zero and one. |
| Observation SD | `sd_form = ~ age_group` | [Formula effects](https://paulregular.github.io/tinyAM/reference/formula_effects.html): formulas modify log-SD, not observation covariance. |

## 3. Diagnose, visualize and compare

```r
check_tam(base, N_iid, M_ar1)
vis_tam(base, N_iid, M_ar1)
```

`check_tam()` separates numerical checks from structural concerns and
statistical advisories. Convergence alone does not establish that the data
identify every component.
Inspect each check before interpreting or ranking an alternative.

In the dashboard, start with **Fishery** and **Surveys** to inspect fits and
standardized residuals. Review **Tables → Diagnostics** and **Parameters** for
checks, uncertainty and process deviations before interpreting **Trends** for
abundance and mortality. Residuals use fitted states and can be affected by
shrinkage.

Test forecasting skill by repeatedly withholding terminal years:

```r
hind <- fit_hindcast(base, N_iid, M_ar1, folds = 5)
hind
```

The summary shows Mohn’s rho (retrospective change; ideally near zero), survey
log-RMSE (lower is better), MASE (below one beats holding the last observation
constant), and scored/failed-fold counts. Compare
validation coverage before comparing scores; failed fits are retained in
`hind$fold_diagnostics`, not scored.

By default, each one-year forecast conditions on actual aggregate catch.
Held-out surveys and catch compositions test prediction; matching the imposed
catch total is not forecast skill. Use `proj_settings = list(F_mult = 1)` for
status-quo F instead. See [hindcast definitions](https://paulregular.github.io/tinyAM/reference/fit_retro.html).

`AIC()` and `cAIC()` complement predictive checks for numerically acceptable
models fitted to the same observations. Intrinsic random walks and
unpenalized initial states need particular care; read the
[information-criterion assumptions](https://paulregular.github.io/tinyAM/reference/information_criteria.html)
before interpreting rankings.

## Further tools

- [Tidy results](https://paulregular.github.io/tinyAM/reference/tidy_tam.html):
  `tidy_tam(base, N_iid, M_ar1)` collects population, parameter and observation tables.
  Confidence limits use the reported scale; `se_scale` identifies the SE scale.
- [Residuals](https://paulregular.github.io/tinyAM/reference/residuals.tam_fit.html):
  `residuals(base)` extracts standardized residuals; `type = "osa"` explicitly
  requests one-step-ahead calculations without refitting.
- [Simulation](https://paulregular.github.io/tinyAM/reference/sim_tam.html) and
  [projections](https://paulregular.github.io/tinyAM/reference/prepare_tam.html)
  explore model behavior and future scenarios.
- [Model equations](https://paulregular.github.io/tinyAM/reference/tinyAM-model.html)
  explain initialization, process scales, boundaries and spawning-time SSB.
