# British Columbia Canary Rockfish: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_canary_bc.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_canary_bc.md).

## Expanded-library review

### Observation translation

The current stock recipe models **female abundance** from 1977–2022 with fitted
ages **1–40+**, retaining source biology and source observations through age
**60+** for tinyAM's internally partitioned older-age group. It uses one
one-year reporting projection. The accepted SS3 model instead carries separate
ages 0–60+ and both sexes over 1935–2022. January spawning and the original
female weight/reproductive-weight functions are preserved, including the
first-mature-age cutoff at 10. Age 0 is immature. There is **no additional 0.5
female fraction**, since female numbers are the model's population states.

`.canary_female_numbers()` calculates each female catch- or index-at-age value
using the **joint-sex** age-composition denominator and the matching original
sex-specific weights before selecting female cells. Neither official fitted
abundance nor fishing mortality is used to construct the observation data.
The sampled age compositions and biology are retained through source age 60;
model ages 40–60+ are aggregated into the model's 40+ population state, with
older-age biology retained internally rather than replacing it by the age-40
weight and maturity.

Pooled catch biomass is allocated using the joint-sex trawl age profile in
each available sampled year (36 samples, 1977–2017). This explicitly assumes
the unsampled non-trawl fleet shares the trawl profile; its total biomass is
about 4% of removals during the modeled period. No missing age samples are
filled. Consequently total catch information in years without trawl age
samples—especially **2018–2022**—does not directly constrain the tinyAM F
likelihood, although the accepted SS3 model conditions its fleet fishing
intensities on the reported catches for those years.

The three surveys with source age compositions (QCS, WCVI, NMFS) are retained.
Commercial CPUE and the HS, WCHG and GIG surveys have no corresponding source
age compositions, so their aggregate index values are not fitted. They remain
available in the canonical assessment inputs. Unlike the SS3 likelihood, the
tinyAM translation fits reconstructed female numbers-at-age as independent
lognormal observations rather than separately fitting aggregate index and
joint-sex Dirichlet-multinomial age-composition likelihoods.

#### Source index uncertainty: retained but not fitted as offsets

The translation brings the accepted run's final year-specific `log_index_sd`
rows into the reconstructed input table. The database converter matches those
values to age-specific index observations by survey and year as
`obs$index$relative_sd`; the recipe checks that the values are positive and
finite. **However, the current fitted `index_settings` omit `sd_supplied`.**
Instead `sd_form = ~ 0 + survey` estimates one lognormal observation SD per
retained survey. Therefore source reweighted SDs are **not used as offsets or
multiplicative SD baselines** in this candidate's likelihood. The source SD
belongs to a whole annual aggregate survey biomass observation, rather than
independent age cells; using it as an age-specific offset would not by itself
recover the native observation covariance.

### Current tinyAM model specification

The exact fitted structure in the accompanying `dfo_canary_bc_2022.R` is:

```r
N_settings = list(process = "off", init = "exp")
rec_settings = list(process = "iid", mu_form = ~ 1, sd = 0.9)
F_settings = list(process = "ar1", mu_form = ~ logistic(age))
M_settings = list(
  process = "off",
  mu_form = ~ 0 + iid(M_prior_group, sd = 0.3),
  mu_supplied = ~ I(0.06)
)
catch_settings = list(sd_form = ~ 1)
index_settings = list(
  q_form = ~ 0 + survey + logistic(age, by = survey),
  q_link = "log",
  sd_form = ~ 0 + survey
)
```

The fitted F surface combines a **static rising logistic mean by age** with a
**correlated residual process across age and year** (`ar1`). This is not the
superseded specification with a year-specific mean RW plus IID residuals; nor
is it simply `logistic(age) + factor(year)` with deterministic F. Unlike the
accepted assessment's two separately catch-conditioned F surfaces, tinyAM
fits one aggregate F surface using Baranov removals. The AR1 departures permit
the effective F-at-age pattern to change over time around its persistent mean.

The possibility that such flexibility approximates **changing relative
contributions of fleets with different fixed selectivity-at-age curves** is
scientifically plausible, but it has not been verified for this case. The
accepted source's trawl and other-fishery curves differ, but the other fleet
accounts for a comparatively small fraction of the total biomass removals;
changes in fleet mixing alone should not be asserted to explain the latest
SSB discrepancy. Distinguishing the effects of fleet mixture, recruitment,
missing terminal catch constraints, and observation likelihood differences
would require targeted sensitivity tests.

The fixed recruitment SD 0.9 approximates source log-recruitment variability
but not Beverton–Holt recruitment, recruitment constraints or bias correction.
The single estimated female log-M effect approximates the source natural-scale
Normal prior (mean 0.06, SD 0.018) with log SD 0.3. No official fitted M value
is substituted as an input or used to tune observations. Survey effective q
has **separate rising logistic age curves** on the **log link**, not a logit
restriction, and is **time invariant**. The base assessment's profiled scalar
q and fixed double-normal selectivity are not reproduced exactly. The
explicitly dome-shaped female selectivities shown in its S09/Run33 sensitivity
are **not the selected asymptotic Run24/24.01 base configuration**.

### Interpretation of the latest fit

The user-provided SSB dashboard shows a materially improved **broad agreement
in absolute magnitude and long-term decline** relative to the previously
reported 60+ candidate. The trajectories nevertheless diverge after roughly
the early 2000s: the official female SSB recovers into the terminal years,
while tinyAM declines further to approximately the mid-2010s and rebounds only
modestly. The user also reports a substantially improved match for **female M**.
These are **qualitative** descriptions of the provided dashboard and user
report, not re-estimated percentage differences, confirmation of convergence
criteria, or independent comparison statistics. The shaded tinyAM interval
is not interchangeable with the official SSB posterior credible interval.

A complete interpretation will require numerical summary statistics and fit
checks from the **current** 40+ AR1 candidate. The official female M posterior
median is approximately 0.09329 (90% limits 0.08094–0.1063); exact current
tinyAM M estimates and confidence intervals have not been provided here.
Official numerical N/F-at-age and annual recruitment series remain unavailable;
missing comparisons are not treated as evidence of agreement.

#### Superseded fit: historical record only

The earlier **60+** translation with free older-age initial states, a different
F process/mean, and source-informed retry was previously reported to converge
after the default initialization failed. Its objective was 1,726.629,
maximum gradient 0.000333, and Hessian positive definite. Historical
comparison statistics were:

| Metric | Superseded comparison | Mean absolute % difference | Median absolute % difference | Terminal % difference | Trend correlation |
|---|---|---:|---:|---:|---:|
| Female SSB | Approximate, 1977–2022 | 27.91 | 23.98 | -59.09 | 0.295 |
| Female M | Shared constant M | 36.51 | 36.51 | -36.51 | Not defined for a constant series |

**None of these numbers describe the current 40+ AR1 model.** They should not
be displayed as current results or reused to support claims about its fit.
Repeated source M-at-age rows correspond to one shared source parameter, not
independent age/year validation.

### Remaining limitations and next diagnostics

The largest structural differences are: omission of catch constraints for
years lacking age samples (notably 2018–2022); omission of aggregate-only
survey and CPUE observations; collapsing the two-fleet, two-sex system into one
female F surface; replacement of native composition/covariance likelihoods
with independent lognormal age-specific observations; shorter initialization
history; and approximate recruitment and M prior models. There is currently
**no identified unique cause** of the residual recent SSB disagreement.

Useful diagnostics, without changing the candidate reported here, would be to
inspect age-specific F deviations and survey residuals over the years when
the curves separate; compare modelled and official annual fleet catch totals;
examine recruitment-cohort contributions to spawning biomass; and run a
well-labelled sensitivity on the F process or original 60+ age boundary.
These are suggestions, not analyses already completed. New quantitative
metrics should be generated from the current fit before replacing any
historical result tables. The canonical source files and observation
transformation are not altered in this documentation update.

### Other superseded fit diagnostics

The 60+ recipe using IID F around `logistic(age) + rw(year)` was also tested
at database revision `00cc998`. Default starts failed the fixed-effect covariance
check; the source-start retry returned optimizer code 1, objective 1934.708,
maximum gradient 3.79e11 and a non-positive-definite Hessian. Its uncertainty
was not usable. These diagnostics do not describe the current 40+ AR1 recipe.
An earlier cached Beverton-Holt fit at revision `442d89` is also not a result
of the current recipe.
