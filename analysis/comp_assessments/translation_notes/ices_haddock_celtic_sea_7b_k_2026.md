# Celtic Sea haddock in divisions 7.b–k, 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_celtic_sea_7b_k_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_celtic_sea_7b_k_2026.md).

## Expanded-library review

### tinyAM fit diagnostics

The package-default initialization was attempted first and failed during optimization with repeated NA/NaN gradient and function evaluations. One retry kept the same observations, formulas, and model settings, using year- and age-aligned official N-at-age, F-at-age, and index-derived q values as starting values. The initial attempts returned no optimizer status, objective, or standard errors, and no dashboard was produced. Their logs are cached locally at `results/cache/ices_haddock_celtic_sea_7b_k_2026/initialization_attempts.rds`.

Review of the failed attempts found a concrete observation mismatch: the translation passed 35 zero catch cells as observations, while the native SAM data object omits their log observations. These are the 33 age-0 cells and two older-age zeros. The stock translation now marks those exact zero catch cells as missing for fitting, while the original zero values remain unchanged in the canonical input table. The default fit and its single source-start retry were rerun after this correction on 2026-10-26; both still failed with NA/NaN gradient or function evaluations. A direct check of the corrected default starting point produced a finite TMB objective (1103.989) and finite gradient, so the invalid evaluations arise later in the optimizer search. The zero-cell mismatch was real, but it was not the only cause. Neither fit attempt returned a final objective, optimizer status, or standard errors, and the assessment is still not ready for a dashboard. No further rescue settings were tried.

The tinyAM process estimates F at age 0, whereas SAM fixes age-0 F to zero. This remains a model-structure difference and was not changed as part of the zero-observation correction.
