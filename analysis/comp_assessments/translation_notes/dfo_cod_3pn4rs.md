# Northern Gulf cod (3Pn4RS) source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_cod_3pn4rs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_cod_3pn4rs.md).

## Expanded-library review

### Initial tinyAM fit

The Sentinel-only translation converged with package-default starting values
and no source-based starts (optimizer code 0; objective 615.876; maximum
absolute gradient 0.000113; positive-definite Hessian and standard errors
available). The fit used 40 fixed and 509 random effects. The dashboard is
cached locally at `results/cache/dfo_cod_3pn4rs_2025/dashboard.html`. An
advisory reports a recruitment AR(1) correlation interval spanning more than
0.5, so recruitment persistence remains imprecise.

Numerical convergence did not imply close agreement with the official outputs.
Across 1995–2024, mean absolute percent differences were 104.9% for approximate
SSB (terminal difference −15.4%; trend correlation 0.295), 43.8% for approximate
recruitment (−47.4%; 0.898), and 27.8% for total abundance (+48.5%; 0.807).
For age-specific series they were 312.4% for N (terminal +48.5%; mean within-age
trend correlation 0.853), 62.2% for F (+55.3%; 0.743), and 57.0% for estimated
M (−16.7%; 0.738). Fbar differed by 52.1% on average and by +270.2% in 2024,
despite a 0.939 trend correlation. The limited index coverage and simplified
population processes remain substantial differences; retain the fit as a
numerically successful but weakly matching first translation.

The cached source files and checksums are listed in
`source_cache/dfo_cod_3pn4rs_2025/manifest.csv`.
`scripts/database/038_import_dfo_cod_3pn4rs_2025.R` reproduces table imports and
digitizes the fitted biological curves from the cached assessment PDF;
`039_validate_dfo_cod_3pn4rs_2025.R` checks coverage, units, provenance, values,
and biological constraints. `057_cache_dfo_cod_3pn4rs_sources.R` retrieves the
additional public survey archives and maturity report and records their URLs
and checksums. The source record remains partial because most official indices
and some exact likelihood settings have not been recovered.
