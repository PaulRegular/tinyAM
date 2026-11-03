# West of Scotland cod: final standalone 2022 assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_cod_west_scotland_2022.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_west_scotland_2022.md).

## Expanded-library review

### Translation and diagnostic outcome

The initial translation retains the source grid, all five surveys, supplied M and biology, q sharing, observation-SD sharing and relative precision. It uses IID N process errors, RW recruitment and age-correlated RW F. Its F states at ages 4–7 remain separate rather than exactly shared, and it omits SAM's catch-scaling factors. These are explicit approximations; the original catches are unchanged.

The first normal-start trial used fixed free initial abundance. It failed numerical checks: optimizer code 1 (`false convergence (8)`), objective 728.1359, maximum absolute gradient 0.0071505 and a non-positive-definite Hessian. An sdreport object was returned but its covariance contains non-finite entries, so usable uncertainty was not established. This trial is retained locally and in Git history.

The retained source-oriented recipe uses `init = "latent"`, which integrates initial older-age abundance without an initial process density. First-year recruitment remains fixed in tinyAM, whereas SAM integrates that state too. Its default-start attempt fails with an undefined optimizer gradient. One retry, with the same observations, formulas and settings, initializes N, F and q from the official estimates; other starts retain package defaults. Both attempts and warnings are retained in `results/cache/ices_cod_west_scotland_2022/initialization_attempts.rds`.

The retry passes numerical checks (optimizer code 0, objective 662.5856, maximum absolute gradient 0.0001743 and a positive-definite Hessian), but is severely implausible as a replication. It estimates enormous N with F close to zero; 231 fixed-parameter pairs have absolute correlation above 0.95. Residual checks flag biased groups, uneven spread, outliers and serial correlation. This is not a successful biological approximation despite its numerical convergence flag.

| Metric | Mean absolute difference (%) | Terminal difference (%) | Trend correlation | Status |
|---|---:|---:|---:|---|
| Total abundance | 120346800 | 203466700 | 0.939 | Matched definition; severe scale failure |
| N-at-age | 4915381000 | 2552743000 | 0.933 | Matched definition; severe scale failure |
| F-at-age | 100.00 | -100.00 | -0.016 | Approximate state sharing; scale and trend failure |
| Recruitment | 52790650 | 91888160 | 0.870 | Matched definition; scale and trend failure |
| SSB | 700540800 | 948653700 | 0.972 | Matched definition; severe scale failure |

Rounded diagnostics are retained here rather than used as tuning targets. Fixed M agrees by construction. The current evidence does not isolate whether omitted catch scaling, different F-state sharing, initial-state integration or their interaction causes the pathological fit. No observation corrections, new priors or core-model changes were introduced to rescue it. The dashboard at `results/cache/ices_cod_west_scotland_2022/dashboard.html` is explicitly labelled as exploratory and unsuitable as evidence of close replication.

Validation: source checks, fitting-disabled translation checks, conditional-start helper tests and multi-assessment runner tests pass. The complete canonical database validates; other stocks and assessments are unchanged, and the importer is byte-for-byte idempotent. The full package test suite and `R CMD check --no-manual --no-vignettes` pass (0 errors, 0 warnings, 0 notes). The runner test emits two package-serialization availability warnings in its temporary parallel workers, without a failed test. No package code or model mathematics changed in this curation/translation work.
