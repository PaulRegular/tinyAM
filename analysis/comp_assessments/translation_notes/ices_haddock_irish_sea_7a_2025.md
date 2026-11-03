# Irish Sea haddock: 2025 age-structured assessment: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_irish_sea_7a_2025.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_irish_sea_7a_2025.md).

## Expanded-library review

### tinyAM translation

The translation fits 1993–2024 at ages 0–5, keeps all four index identities and their reported ages, uses the annual weight and maturity tables, and supplies the reconstructed fixed M values. The model uses an IID N process with exponential initialization, random-walk recruitment, correlated random-walk F, one catch observation SD and survey-specific q and observation SDs. This is an approximation of the single-fleet ASAP fit; the report does not give enough detail to reproduce every likelihood, variance-sharing rule or initial-state treatment. The survey CV and effective-sample-size settings are recorded in the audit but are not fully represented by tinyAM's simple observation-error structure.

The report contains zero catches and zero index values, which tinyAM treats as missing under its lognormal observation likelihood. Those source values remain in the database. The three missing-index groups (Q4 2024, NI-MIK 2020/2022 and UK-FSPW 2014) are omitted only from the analysis-side observation table; no values are imputed. Fitting uses package starting values first. Since official N-at-age and F-at-age are unavailable, a source-aligned population-start retry cannot be constructed from the published outputs.

The package-default fit converged with optimizer code 0, maximum absolute gradient 0.00087 and a successful positive-definite sdreport. Agreement is still poor: mean absolute percent differences are 178.9% for recruitment and 523.5% for the approximate SSB comparison; terminal differences are +258.5% and +842.3%, respectively. Recruitment trend correlation is 0.695. Official Fbar is not directly compared because its age range matches but its weighting cannot be aligned without official N-at-age. The fixed M comparison checks only that the reconstructed 2016 input was passed through; it is not an official fitted M output.

### Alternative process checks

Turning the IID N process off, to approximate deterministic cohort survival,
returned optimizer code 0 and a small gradient, but the fixed covariance was
non-finite and the Hessian was not positive definite. Keeping IID N while
replacing correlated F increments with independent temporal random walks did
pass convergence checks (objective 903.142, gradient 0.001067 and a positive-
definite Hessian), but worsened both reported comparisons: median percent
differences rose from 178.6% to 208.7% for recruitment and from 562.9% to
620.9% for SSB. Recruitment correlation increased only slightly (0.695 to
0.701), while SSB correlation fell from 0.922 to 0.915. The original
correlated-F translation is retained.
