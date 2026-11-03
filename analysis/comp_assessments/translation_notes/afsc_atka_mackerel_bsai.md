# Bering Sea/Aleutian Islands Atka mackerel: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_atka_mackerel_bsai.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_atka_mackerel_bsai.md).

## Expanded-library review

### Current tinyAM translation (revised October 2026)

The working translation continues to cover **1990–2023 and ages 2–11+**. The
original early catch-at-age rows are incomplete, and catch age 1 is not reported.
Although the database retains the 2024 catch and survey biomass inputs, the
age-index reconstruction cannot include those aggregate-only observations
without additional assumptions about age composition. No missing ages are
manufactured or sourced from official fitted population states.

Biological inputs and observational conversions are unchanged: original population,
catch and survey weights; deterministic cohort survival; fixed M = 0.3;
source survey time 6.5/12; female maturity multiplied by 0.5; and August 1
spawning time 7/12. Initial abundances at older ages are free. Age-2
recruitment is IID around a constant estimated mean with **estimated process SD**,
not the source's age-1 Beverton–Holt relationship with fixed steepness 0.8.

#### Fishing mortality: logistic mean and separable AR1 residuals

The current working fit uses:

```r
F_settings = list(process = "ar1", mu_form = ~ logistic(age))
```

This is **not** the previous `cor_rw` candidate and not a fixed
`logistic(age) + factor(year)` model. The rising logistic function describes
the **mean F-at-age on the log scale**. A stationary Gaussian AR1 residual
surface, correlated across both ages and calendar years, permits realized F to
depart from that mean, including shifting age-specific vulnerabilities.

The accepted AMAK model instead uses annual fishery intensity and penalized,
time-varying nonparametric selectivity, shared across ages 10+. The current
tinyAM approximation does not separately normalize fishery intensity and
selectivity, does not impose source age-10+ equality in F and does not recreate
the accepted penalties. AR1 flexibility is a pragmatic representation of
age/time variation; it is not evidence that any particular source selectivity
change was recovered.

#### Survey q: monotone age curve shared from age 7

The fitted survey observation setup is:

```r
obs$index$q_age <- pmin(obs$index$age, 7)   # Retained but not used in current q_form
obs$index$q_block <- factor(pmin(obs$index$age, 7))

index_settings = list(
  q_form = ~ 1 + mono(q_block),
  q_link = "log",
  sd_form = ~ 0,
  sd_supplied = ~ relative_sd
)
```

`mono(q_block)` adds **non-negative successive increments on the log-q scale**.
Consequently, effective survey q is **non-decreasing across ages 2–7** and is
**identical at ages 7–11+**. This differs substantively from the previous
`rw(q_age)` candidate: there is no age random-walk q process in the current
fit, and the curve cannot turn downward. There is also **no calendar-year q
effect**. The unused `q_age` column is retained in the current uploaded script
but does not drive fitted q; `q_block` does.

The earlier use of shared q at ages 10+ was closer to the accepted assessment's
age-sharing convention but allowed extreme q estimates in this approximation.
The decision to couple from **age 7** was made to stabilize high-age q; it is
an explicit translation choice, **not** a statement that the accepted model
assumes selectivity flattens at age 7. The source survey selectivity is
nonparametric and can be dome-shaped; its scalar q has a prior and is
normalized jointly with selectivity over ages 4–10. The tinyAM monotone
curve cannot reproduce descending source selectivity or identify a directly
comparable normalized scalar q. `q_link = "log"` does not bound q below 1.

#### Observation uncertainty and fitting

Published rounded annual survey CVs are converted to lognormal observation SDs:
`sqrt(log1p(CV^2))`. The `sd_supplied` offset is used **without estimating a
common multiplier** (`sd_form = ~ 0`). The source CV describes the annual
**aggregate survey biomass**, however, and is repeated across reconstructed
age-specific indices. The independent lognormal age likelihood does not
reproduce the source's separate aggregate survey/index and multinomial age
composition likelihoods, their effective sample sizes, or within-year
covariance. Catch-at-age uses an estimated common lognormal SD.

The present R script calls **one** `tinyAM::fit_tam()` model directly. It has
no fallback fit and does not initialize from official N-at-age or F-at-age.
The fitted M = 0.3 is externally fixed, so agreement in M is an input check.

#### Performance and comparison boundaries

The user reports that the current specification reaches a reasonable practical
match and that coupling q at **age 7** prevents extremely high estimates.
No optimizer diagnostic table, current parameter estimates or recalculated
agreement statistics accompanied these two revised files, so this assessment
of the current fit is **qualitative**, not a newly quantified validation.

The previous fit's convergence (objective **511.565**, maximum gradient
**0.000267**, positive-definite Hessian), and historical comparisons (modeled-age
SSB MAPD **25.7%**, terminal difference **−23.6%**, correlation **0.966**)
refer to a **superseded independent-F-RW / independent age-q translation with
source-informed retry**. They do **not** describe the present AR1-F,
monotone-q fit. Additional previously reported differences for N and F also
belong only to that earlier translation.

The dashboard retains native age-1+ SSB and a separate **approximate modeled-age
August 1 female SSB** derived from rounded official N-at-age, F-at-age,
fixed M and recorded biology, with no reconstructed uncertainty. Native
SSB intervals must not be applied to these derived sums. The original
source estimates of N and F are used for comparisons only, never as
observations or as current-fit initialization.

The key remaining structural limitations are the different catch/survey
likelihood, unavailable historical and 2024 aggregate-only age information,
missing native composition effective sample sizes and ageing error, lack of
a separately normalized scalar survey q, the inability of monotone q to
represent dome-shaped survey selectivity, and a simpler recruitment process.
Further fine-tuning can be assessed case by case; the current translation is
presented as a practical approximation rather than a reconstruction of the
accepted AMAK likelihood.

### Earlier age-specific mean AR1 trial

A separate superseded trial used an age-specific mean (`~ factor(age)`) with
AR1 F residuals and q shared from age 10. It converged from default starts
(objective 574.295, maximum gradient 0.000834, positive-definite Hessian and
successful uncertainty estimation). Its approximate SSB mean absolute percent
difference was 22.0%, terminal difference -25.9%, and trend correlation 0.951.
These results describe that trial only, not the current logistic-mean F and
monotone-q specification.
