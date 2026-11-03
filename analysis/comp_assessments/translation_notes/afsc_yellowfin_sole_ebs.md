# Eastern Bering Sea yellowfin sole: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_yellowfin_sole_ebs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_yellowfin_sole_ebs.md).

## Expanded-library review

### tinyAM translation

The recipe represents female ages 2–20+ over 1982–2023. Joint-sex original
age samples and sex-specific weights convert the total biomass observations
to female numbers before filtering ages. Male observations are omitted from
the fit; both sexes remain in the canonical tables. This avoids imposing an
equal adult sex ratio on spawning biomass or averaging fixed female M with
estimated male M. Comparisons of N, total abundance and recruitment are
explicitly female quantities, not combined-sex stock totals.

The recipe retains original female biology, May 1 spawning, fixed female M,
the three centered q covariates and supplied biomass log SDs. Its independent
F random walks, constant-mean age-2 recruitment, rising logistic survey curve
and age-specific lognormal observations approximate the source selectivity,
Ricker and composition likelihoods. The source's different times for biomass
and age proportions remain a limitation. Full native outputs are retained;
modeled-age female SSB is derived without uncertainty and checked against
native SSB using a bound for printed N/F/SSB rounding.

The latest tracked fit diagnostic uses database revision `00cc998` and records
two failed attempts. Package starts failed the fixed-effect covariance and Hessian
checks (gradient 8.06e11). The aligned female N/F source-start retry returned
code 1, objective 9602.126, gradient 6.99e11 and a non-positive-definite
Hessian. An SD-report object exists, but fixed-effect uncertainty is unavailable.
The fit is not interpretable as a successful assessment.

The local cache contains different objectives (12255.88 and 10402.38) from the
tracked diagnostic, and the source review previously cited yet another pair.
Those fit and dashboard artifacts are not reconciled to the latest run and must
not be used as current comparisons. The shared comparison summary and population
series contain no valid yellowfin rows.

A cached exploratory dashboard and historical comparison snapshot remain
locally available for diagnosis only. They are not part of the current shared
comparison tables and do not establish agreement. Potential
limitations include the initial abundance/F boundary, independent F states,
the source's inconsistent composition/index timing and the omitted joint
composition likelihood. These results do not isolate a single cause.

An exponential-initialization sensitivity with package-default starts and all
other settings unchanged was stopped after more than four minutes without a
returned fit. It produced no diagnostics and remains unevaluated; no conclusion
about that initializer is drawn and no model change is retained.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### tinyAM translation

Both permitted initialization attempts failed. Default starts returned code 1
(false convergence), objective 7723.736 and maximum gradient 9.75e11 after
625 seconds. One aligned female N/F retry returned code 1, objective 10453.269
and maximum gradient 1.49e11 after 277 seconds. Both Hessians were not positive
definite and fixed-parameter uncertainty was non-finite. An SD-report object
exists, but that is not successful uncertainty estimation. Initial age-2
abundance inflated to about 1.3e13 fish and recruitment log SD to about 19 in
the returned retry; these estimates are not interpretable as a successful
assessment. No further retries or model changes were made.

An exploratory local dashboard and comparison summaries are retained to expose
the failure. They must be read with the failed numerical diagnostics; even
apparently close metrics do not establish agreement from this fit. Potential
limitations include the initial abundance/F boundary, independent F states,
the source's inconsistent composition/index timing and the omitted joint
composition likelihood. These results do not isolate a single cause.
