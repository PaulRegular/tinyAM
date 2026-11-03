# Eastern Bering Sea pollock: official-assessment source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_pollock_ebs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_pollock_ebs.md).

## Expanded-library review

### Final tinyAM translation (current candidate)

The attached stock translation retains model years **1964–2024**, ages **1–15**,
and age **15+** as the terminal plus group. It retains the source's fixed
age-dependent natural mortality: M = 0.9 at age 1, 0.45 at age 2 and 0.3 at
ages 3–15. Annual stock weights, source-specific survey weights, and source
maturity are supplied through the translation. The maturity vector is
multiplied by **0.5** for female biomass and SSB uses **`spawn_time = 0.25`**.

The fitted settings in the updated stock recipe are:

```r
N_settings = list(process = "off", init = "exp")
F_settings = list(process = "cor_rw", mu_form = NULL)
M_settings = list(process = "off", mu_form = NULL,
                  mu_supplied = ~M_assumption)
catch_settings = list(sd_form = ~1)
index_settings = list(
  q_form = ~0 + survey + logistic(age, by = survey),
  sd_form = ~0 + survey
)
# rec_settings is omitted: tinyAM's default recruitment random walk applies.
```

**Population process and scale.** Older cohorts follow deterministic survival
through F and M, without IID abundance-process deviations. Initial
older-age abundance uses tinyAM's exponential survivorship parameterization.
This remains different from the source's penalized initial-age deviations,
and tinyAM's default recruitment random walk is not the native recruitment
parameterization. Earlier tinyAM trials included IID N deviations; removing
those deviations was the principal change associated with much better
agreement in scale, according to the subsequent model-development comparison.
This observation is qualitative until the updated diagnostics and comparison
summary are regenerated.

**Fishing mortality.** The current translation uses `cor_rw`, with an
age-correlated innovation vector for annual changes in log F; the model does
not impose a fixed F selectivity curve. In contrast, the official model
factorizes annual fishing intensity and age-specific fishery selectivity,
with selectivity itself varying over time under regularization. Thus the
current correlated F random walk is a parsimonious approximation, not an
exact representation of the native fishery process.

**Catch observations.** The official fit uses total fishery catch biomass
and an effective-sample-size-weighted age-composition likelihood, rather
than direct lognormal catch-at-age observations. tinyAM reconstructs catch
numbers-at-age and fits the positive observations with **one shared lognormal
catch SD** (`sd_form = ~1`), without the earlier quadratic age effect.
Reconstructed catch observations below **0.02** are set to zero to align
with rounded zero cells in the report; those zero cells are not fitted as
positive lognormal observations. The source's missing 2024 fishery
composition is not fabricated.

**Bottom-trawl and acoustic-trawl surveys.** The native BTS model includes
a full across-year VAST biomass covariance, age-composition likelihood,
and time-varying logistic selectivity with a separate treatment of age 1.
Native ATS uses biomass and age-composition likelihoods with time-varying,
penalized age-selectivity coefficients rather than a static logistic curve.
The current tinyAM translation replaces these components with reconstructed
mid-year age-specific indices. Each survey receives its **own catchability
baseline and static logistic age curve** through
`~0 + survey + logistic(age, by = survey)`, and its own independent lognormal
observation SD through `~0 + survey`. There is **no temporal survey-q random
walk** in the final translation. The former survey-by-age q blocks and their
pooled ages 9+ are not used by the current fitted formula.

The source treats ATS age 1 as a separate recruitment index, whereas tinyAM
combines available ATS age-1 values with the main acoustic-trawl age-specific
index. The 2024 age-1 value remains excluded in both representations because
it fails the native uncertainty criterion. Historical fishery CPUE and
acoustic-vessel-of-opportunity indices are omitted from tinyAM's age-specific
translation because they lack the associated age-composition information.

**Comparison definitions.** Reported official N for ages 1–9 is compared
directly with tinyAM, whereas the official reported **10+** abundance is
compared with the sum of tinyAM ages **10–15**; it is not the model's age-15
plus group. Age-1 recruitment, total abundance over ages 1–15, and reported
age-3+ biomass are matched on their corresponding age definitions. Female
SSB now also uses the same **0.25-year spawning fraction** and **0.5 female
maturity multiplier** as the source, so it is no longer labelled an
approximate comparison *because of spawning time*. Differences in fitted
SSB may still arise from different mortality, population, and observation
models. Official F-at-age is reconstructed from the pinned native fit and
source code; F uncertainty remains unavailable. The fixed source M vector
matches by construction.

**Updated numerical validation.** The current default-start fit converged
(optimizer code 0, objective 1,911.58, maximum absolute gradient 0.000831,
positive-definite Hessian, successful standard-error estimation, 13 fixed
parameters and 975 random effects). Over the common fitted period, median
absolute percent differences were 10.72% for SSB, 19.45% for N-at-age,
10.32% for total abundance, 8.42% for matched age-3+ biomass, 14.01% for
recruitment and 28.32% for F-at-age. Their trend correlations were respectively
0.958, 0.940, 0.911, 0.951, 0.917 and 0.696. SSB is a definition-matched
comparison under the stated spawning-time and female-maturity conventions;
the fixed-M comparison is identical by construction.

F-at-age remains the clearest discrepancy: its mean absolute percent
difference is 43.39%, terminal-year difference is +41.31%, and trend
correlation is 0.696. Recruitment has a large terminal-year difference
(+71.13%) despite a 14.01% median difference and 0.917 trend correlation;
total abundance also has a +34.26% terminal difference. These metrics warrant
review and are not hidden by the closer aggregate SSB result.

An AR1 F process around age-specific means was tested as a source-plausible
alternative. It converged (objective 1,900.75, gradient 0.00153, positive-
definite Hessian and successful uncertainty estimation), but F-at-age trend
correlation fell to 0.665 while its median percent difference improved only
slightly (27.83% versus 28.32%). SSB median difference decreased to 9.72%, but
its trend correlation also decreased to 0.954. The alternative does not resolve
the F trend discrepancy and is not retained; the correlated random walk
remains the current translation. These results supersede earlier diagnostic
and comparison statistics for the older specification.

The current candidate is a deliberately **simplified structural
translation**, not a claim of likelihood-level equivalence to the operational
assessment.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### Final tinyAM translation (current candidate)

**Updated numerical validation pending.** The numerical diagnostics and
agreement statistics previously recorded here (including the objective
1,917.084, 27 fixed parameters, 1,815 random effects, and the corresponding
scale differences and correlations) describe an **older tinyAM configuration**
with IID N deviations, independent F random walks, quadratic catch SD, and
survey-by-age q blocks. They should **not** be attributed to this revised
model. Once the updated fit and comparison tables are available, regenerate
`fit_diagnostics.csv`, `comparison_summary.csv` and the relevant figures, and
insert the current convergence, SSB-scale, abundance, recruitment and F-at-age
comparison results here. No updated numerical performance claims have been
introduced without those outputs.
