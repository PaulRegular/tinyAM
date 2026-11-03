# Northern Norwegian coastal cod, 2026 ICES assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_cod_northern_coastal_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_northern_coastal_2026.md).

## Expanded-library review

### tinyAM translation

The translation fits 1995-2025 at ages 3-10. It excludes age 2 so the youngest tinyAM age corresponds to official recruitment age 3. The assessment's age-10 plus group remains a terminal modeled group in tinyAM; hidden ages within that group are unavailable.

The fit uses total catch-at-age and the age-specific swept-area index. It retains official stock weight, maturity, and supplied natural mortality. M is not estimated. Recruitment and abundance residuals use random walks, while F uses correlated random-walk increments and the source mean-F ages 4-8. q groups reproduce the reported swept-survey sharing at ages 3, 4, 5-9, and 10+. The fit uses one independent lognormal SD for the age-specific swept observations; SAM's cross-age observation correlation is not reproduced.

The aggregate acoustic biomass index is not fitted because tinyAM supports age-specific abundance-index likelihoods rather than an aggregate biomass index. No source observations are changed to improve agreement. The current SAM run's initial-state treatment, exact variance sharing, shared F states at ages 7-9, and full likelihood are not reproduced.

The dashboard's common-age SSB, abundance, and biomass comparisons use ages 3-10; native SAM aggregates remain separately reported. The fixed M comparison checks that the same source surface reaches tinyAM and is not evidence that M was estimated. The assessment's arithmetic Fbar is retained as a native output; the common-definition comparison uses abundance-weighted F over ages 4-8.

For comparison reporting, biomass-at-age is derived as official N-at-age in thousand fish multiplied by official stock weight in kilograms per fish, giving tonnes. Age 10 retains the source 10+ group. This comparison-only series is not added to the canonical database, and no uncertainty is assigned because age-specific N uncertainty is unavailable. The refreshed fit remains numerically converged (optimizer code 0, objective 154.1478, maximum absolute gradient 0.0003204, successful standard-error report, positive-definite Hessian). The dashboard is regenerated from this fit and includes the derived biomass-at-age comparison.
