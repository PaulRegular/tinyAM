# Bering Sea/Aleutian Islands Pacific ocean perch: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_pop_bsai.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_pop_bsai.md).

## Expanded-library review

### tinyAM translation

The recipe uses 1960–2024 and ages 3–40+, the observed/exported terminal age
group. It reconstructs age-specific catch and survey indices from the original
biomass, compositions and weights. Source length samples, ageing error,
composition likelihoods, spline penalties and priors are omitted. Survey time
0.5 approximates the source's annual-average survival. Recruitment has separate
mean periods before/after 1990, with estimated IID variation; survival is
deterministic, initial abundance is free, F has independent temporal random
walks, and M is an estimated constant without the source prior. Maturity is
held at the official fitted logistic curve, with a female fraction of 0.5
and spawning time 0.25. These choices are explicit in the dashboard background.

Preparation and the default-first/retry, plus-group, reference-alignment and
comparison-definition regression tests pass. The first default-initialization
attempt was manually stopped after more than seven hours without a returned fit;
numerical checks were unavailable. The single permitted source-start retry began
at 17:36 local and used aligned official N and F values. The attempt cache
recorded no returned result by 18:52 local on 2026-10-20, when no R/Rscript
process remained. Exact termination time and cause are unknown, so the objective,
gradient, standard errors, Hessian and convergence are unavailable. No agreement
statistics or dashboard are claimed. A later default-initialization attempt was
started on 2026-10-23. Its current attempt record still says that it started
without returning a result, so its objective, gradient, standard errors, Hessian
and convergence remain unavailable. No further retry has been made while this
attempt is pending, and no data or model assumptions were changed.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### tinyAM translation

Preparation and the default-first/retry, plus-group, reference-alignment and
comparison-definition regression tests pass. The full default-initialization
pipeline was launched at the translation checkpoint but had not returned a fit
after more than seven hours of CPU time. It was terminated as an execution
failure. Objective, gradients, uncertainty and convergence are unknown; no
agreement statistics or dashboard are claimed. The helper then in use retained
attempts only when the complete pipeline returned, so its active initialization
stage and any earlier completed-attempt diagnostics could not be recovered.
No additional retry was launched because that could exceed the single-retry
policy. The helper now saves each attempt before and after execution so later
interruptions cannot lose this history. No data or model assumptions were
changed to rescue this case.
