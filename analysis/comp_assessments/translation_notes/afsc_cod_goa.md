# Gulf of Alaska Pacific cod: official-assessment source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_cod_goa.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_cod_goa.md).

## Expanded-library review

### Current tinyAM candidate translation and fit

The translation fits 2007–2025, when fishery length compositions and conditional
age-at-length samples overlap. Separate fishery-specific pooled age–length keys
are estimated from supplied samples (through 2024), used to convert annual length
compositions and biomass catches into catch-at-age, and then the three fisheries
are aggregated. Uncovered length-bin mass is omitted and reconstructed age shares
are normalized. The official likelihood is based on fleet catches, length
compositions and conditional ages, not these pseudo-observed catch numbers at age.

Only the age-sampled NMFS bottom-trawl survey is translated. Its aggregate index
is allocated across ages with the pooled age–length key; the original aggregate
log SD is repeated across reconstructed age rows. Month-7 observations use an
approximate year fraction of 0.5. The longline survey is not fitted because
age-composition coverage required by this translation is unavailable.

Population weights/maturity are reconstructed from official fitted mean-growth
parameters and the source weight–length and maturity functions. The approximation
does not integrate length variability within age. The female maturity multiplier
is 0.5. Source M estimates (including the 2014–2016 elevated M block) are supplied
and held fixed. SSB in the original assessment includes age 0, whereas the tinyAM
translation starts at age 1, so agreement is approximate by definition.

#### Candidate specification (requires a new fit)

```r
N_settings     = list(process = "off", init = "exp")
# Default tinyAM recruitment RW retained (not source BH recruitment)
F_settings     = list(process = "cor_rw", mu_form = NULL)
M_settings     = list(process = "off", mu_form = NULL,
                      mu_supplied = ~ M_assumption)
catch_settings = list(sd_form = ~ 1)
index_settings = list(
  q_form = ~ logistic(age) + rw(year),
  sd_form = ~ 1,
  sd_supplied = ~ relative_sd
)
```

This paired-specification experiment replaces the earlier **IID F deviations
around age means** (`process = "iid", mu_form = ~ factor(age)`) with correlated
annual F random walks. In contrast to a fixed logistic F-age curve, each age can
evolve through time while annual changes are correlated across ages. The survey q
model now combines a smooth ascending logistic age curve with a shared year random
walk on log catchability. That RW affects the *survey level* and is **not**
equivalent to the original length-based time-varying survey selectivity. Unlike the
source, the small translation has one aggregated fishery and does not reproduce
the conditional age-at-length likelihood.

Earlier candidate exploration suggested that additional **F flexibility** was
important for agreement, but did not establish that `cor_rw` outperforms IID F.
This new candidate has not been fitted or diagnosed as part of this revision.
The GOA translation retains exponential older-age initialization because the
previous free-initial-state configuration was numerically unreliable.

#### Historical comparison from the preceding F-IID candidate

The previously documented `F_settings = list(process = "iid", mu_form = ~ factor(age))`
candidate reported optimizer code 0, relative convergence 4, objective 274.38,
maximum absolute gradient 0.000411, positive-definite Hessian, 25 fixed and 208
random effects. Its approximate 2007–2025 aggregate SSB comparison reported
mean absolute difference 24,893 t (32.10%), terminal difference +7.73%, and
trend correlation 0.945. Those results belong to that earlier fit and should
**not** be used to describe the new `cor_rw` F / logistic-plus-year-RW q model.

Official N-at-age and F-at-age surfaces are now available from the reconstructed
native report. The M comparison is equal by construction because numerical
official M was supplied; recruitment and total biomass remain
definition-non-equivalent because of source age 0.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### Current tinyAM candidate translation and fit

#### Candidate specification (requires a new fit)

#### Historical comparison from the preceding F-IID candidate

Official N-at-age and F-at-age surfaces remain unavailable. The M comparison
is equal by construction because numerical official M was supplied; recruitment
and total biomass remain definition-non-equivalent because of source age 0.

### Historical annual-key sensitivity

A separate earlier fit used annual age–length keys for fisheries where age
samples were available and pooled keys for missing years and bins, while keeping
the survey key pooled. It converged (optimizer code 0, maximum absolute
gradient 0.000137, positive-definite Hessian; objective 303.51). Its approximate
SSB comparison had a 16.37% median and 23.22% mean absolute percent difference,
a +3.97% terminal difference, and 0.928 trend correlation. These results are
for that annual-key sensitivity and the preceding F-IID/q-block specification;
they are not fit results for the current pooled-key, correlated-F, logistic-age
and year-RW q candidate.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/afsc_cod_goa.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The translation uses a pooled age-length key and does not preserve the native conditional-composition likelihood.
