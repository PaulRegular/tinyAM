# Norwegian spring-spawning herring: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_herring_norwegian_spring_spawning.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_herring_norwegian_spring_spawning.md).

## Expanded-library review

### tinyAM translation

The current-run approximation is defined in
[preliminary recipe](../scripts/translation/preliminary/ices_herring_norwegian_spring_2026.R). It fits
1988–2025 at ages 2–12, using IID N with free initial abundance, an F random
walk, and the standard age-specific M pattern as fixed. Catch and the four
numerical survey series are fitted with tinyAM observation errors; SAM's
relative-error weights and unavailable annual M deviations are not reproduced.
The fit uses package-default starting values; no source-based retry was needed.

The fit converged (optimizer code 0, objective 1059.441, maximum absolute
gradient 0.00111). The Hessian was positive definite and `sdreport` succeeded;
the fit had 53 fixed and 825 random effects and took 3.6 seconds. Numerical
convergence does not remove the model differences above. Residual checks flag
some groups with unusual standardized-residual spread, residuals above three
standard deviations, and adjacent-year residual correlation.

Over 1988–2025, common-definition SSB had 10.75% mean absolute percent
difference, +0.16% terminal difference, and 0.965 trend correlation.
Recruitment had 46.22%, +50.21%, and 0.743; total abundance had 20.03%,
+8.55%, and 0.798. N-at-age had 21.82% mean absolute percent difference and
0.943 mean within-age trend correlation; F-at-age had 29.27%, -21.82%
terminal difference, and 0.816. Fbar had 13.59%, -0.90%, and 0.923. Fixed M
matches its supplied baseline by construction and is not evidence that tinyAM
estimated or reproduced the source's annual M deviations. These comparisons
describe this documented approximation; they are not fit targets.

The fit, initialization diagnostics, and dashboard are cached locally under
`results/cache/ices_herring_norwegian_spring_2026/`.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_herring_norwegian_spring_spawning.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The published age-12 stock weight and maturity are used as proxies for the
12+ group because no abundance-by-single-age surface above age 12 supports a
defensible weighted biological average. Maturity cohorts are mapped to
calendar year by `cohort = year - age` in the translation.
