# Faroe saithe 2025 (`pok.27.5b`): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_saithe_faroe_2025.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_saithe_faroe_2025.md).

## Expanded-library review

### tinyAM translation

The translation retains both surveys and their source timing, uses the package's normal initialization first, and permits one retry with year/age-aligned SAM estimates as starting values. It uses IID N survival deviations, random-walk recruitment, a correlated F random walk, and the official fixed M. Survey residual age correlation and shared oldest F states remain approximations.

During the 2026-10-23 source audit, catch and survey observations were found to have been flattened in age-major order against year-major labels. The importer now matches every observation to the native SAM year, age and fleet map before writing; regression checks include off-diagonal source cells and preserve the source zero catches. Survey numerators are corrected with their matching indices; effort, biological inputs and official outputs retain their source values. The importer can replace this assessment's existing rows without changing other assessments.

Earlier tinyAM results used the misaligned observations and are invalidated by this repair. Their local cache is retained for provenance. With corrected observations and unchanged model settings, the normal-start fit converges: optimizer code 0 (`relative convergence (4)`), objective 1296.036, maximum absolute gradient 0.001226, successful SD estimation and a positive-definite Hessian. No source-based starting-value retry was needed. Checks still flag one fixed-parameter correlation above 0.95, uneven residual spread and residual outliers.

The refreshed dashboard is `results/cache/ices_saithe_faroe_2025/dashboard.html`. Numerical convergence does not resolve the remaining disagreement:

| Metric | Mean absolute difference (%) | Median absolute difference (%) | Terminal difference (%) | Trend correlation |
|---|---:|---:|---:|---:|
| Total abundance | 36.80 | 35.76 | 35.39 | 0.949 |
| N-at-age | 1348.23 | 105.01 | 165.11 | 0.924 |
| F-at-age | 54.24 | 53.82 | -50.42 | 0.679 |
| Recruitment | 22.59 | 18.96 | 33.32 | 0.962 |
| SSB under shared biology | 103.23 | 99.64 | 103.80 | 0.853 |

These core metrics exceed the scale-review thresholds; F and SSB also fail the 0.90 correlation screen. The retained approximation omits survey residual age correlation and separately estimates F states at ages 11–15 rather than sharing the official oldest state. These are unresolved structural differences, not reasons to modify source observations or interpret the fit as equivalent to SAM. Fixed M agrees by construction and is not an estimation result.

Biomass and SSB comparisons recompute official N with the same translated biology used by tinyAM: original annual stock-weight inputs through 2024 and the fitted SAM stock weights in 2025. They isolate differences in population estimates under shared biology. Native SAM biomass and SSB use its fitted stock weights and are retained separately. The common Fbar comparison uses population weighting over ages 4–8 for both models; native SAM arithmetic Fbar remains available for source context. Derived comparison aggregates have no invented uncertainty.

For the dashboard only, official predicted indices and q estimates, including their intervals, are multiplied by 0.001 to match tinyAM's thousand-fish index scale; q is correspondingly reported per effort-hour. The canonical database retains the original official values and units. The translation test checks these conversions and verifies that the source outputs remain unchanged.

Validation after the repair: the source and translation regression checks and full package test suite pass, as does `R CMD check --no-manual --no-vignettes` (0 errors, 0 warnings, 0 notes). Canonical validation at the repair checkpoint reports 114 stocks, 99 assessments, 2022 assumptions, 346474 inputs and 184026 outputs. Re-running the importer leaves every canonical file byte-for-byte unchanged; records belonging to other assessments are unchanged by the repair.
