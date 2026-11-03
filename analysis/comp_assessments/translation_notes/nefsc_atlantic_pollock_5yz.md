# Atlantic pollock (5YZ): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/nefsc_atlantic_pollock_5yz.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/nefsc_atlantic_pollock_5yz.md).

## Expanded-library review

### tinyAM translation choices

The translation combines catch and dead-discard numbers-at-age into one total-removals series. It uses the four active survey series and source survey timing, but approximates the source's biomass-index/composition likelihood as age-specific abundance indices. The model uses one F process in place of ASAP's six selectivity blocks and estimates survey-by-age q effects in place of the source survey-specific annual q deviations. It also uses tinyAM's own process and observation-error structures. The source totals, age compositions, biological surfaces, and official outputs remain unchanged in the canonical database; transformations for fitting are performed only in the stock translation script.

### tinyAM fit status

The 2026-10-21 fit did not pass numerical convergence checks. The package-default initialization returned optimizer code 1 (`false convergence (8)`), objective 3154.205, maximum gradient 168,801, and a non-positive-definite Hessian with non-finite fixed-parameter covariance. The single policy-permitted retry with year- and age-aligned source N and F starts also failed: optimizer code 1, objective 3140.939, maximum gradient 5.57 billion, NA/NaN function evaluations, and a non-positive-definite Hessian. An exploratory dashboard was rendered from that failed fit for inspection only; its large differences are not a valid agreement result. The fit and attempt diagnostics are cached locally under `results/cache/nefsc_pollock_5yz_2024/`.
