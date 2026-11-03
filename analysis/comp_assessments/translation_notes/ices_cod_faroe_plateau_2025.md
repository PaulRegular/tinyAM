# Faroe Plateau cod: 2025 ICES SAM assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_cod_faroe_plateau_2025.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_faroe_plateau_2025.md).

## Expanded-library review

### tinyAM approximation

The translation retains the model years and age groups, both survey identities and timing, observed catch and index scales, source age-sharing for q and observation SDs, available annual biology, and fixed source M. It uses a separate random-walk recruitment process and one shared IID SD for older cohort-survival deviations, plus correlated random-walk F with mean F defined over ages 3–7. The recruitment and survival processes each have their own SD, corresponding to SAM's separate age-1 recruitment variance and shared age-2–10 survival variance. These settings do not reproduce SAM's initial-state treatment or full covariance conventions.

tinyAM has independent index residuals, so it does not represent the summer survey's AR(1) age correlation. The source model has no estimated F state at age 1; tinyAM's F model includes age 1, so that boundary is an approximation. The source model's observation slots omitted because they are missing remain missing in tinyAM. Source values are used only as comparison references. The fit tries package defaults first and uses one aligned source-based retry only if numerical checks fail.

### tinyAM fitting result

The 2026-10-22 attempt failed with package defaults, then failed again with one source-based retry. Both `nlminb` calls stopped after NA/NaN function or gradient evaluations; neither returned optimizer diagnostics, a fixed-effect gradient, or an uncertainty report. The initial objective and gradient evaluations were finite for both starting vectors, so the invalid evaluations arose during the optimizer search. No fitted tinyAM model or dashboard was produced, and no additional starting-value retries or model changes were made.
