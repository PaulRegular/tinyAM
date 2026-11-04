# Western Baltic cod assessment history: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_ices_wgbfas_morhua_bs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_ices_wgbfas_morhua_bs.md).

## Expanded-library review

### tinyAM 2026 comparison

The restricted recipe uses 1985–2023 and ages 1–7+. This avoids missing age-0 maturity and the missing age-7+ stock weights in 2024–2025. It retains BITS Q1 and Q4 for 2001–2023 at ages 1–3; the reported age-4 row is pooled 4+ and is not treated as a single age. FEJUCS is omitted because it indexes age 0. Survey timing is carried forward from the 2025 native SAM fit (Q1 0.1653494, Q4 0.8477816) because the 2026 report does not provide within-quarter fractions. SSB is calculated at the start of the year by setting spawning time to zero. The fitted model uses IID N and F processes, IID age-1 recruitment, fixed source M, and the package's default initialization; no source-based starting values were used.

The default fit passed numerical checks: optimizer code 0, objective 455.812, maximum absolute gradient 0.00066, successful uncertainty estimation, and a positive-definite Hessian (23 fixed and 539 random effects). `check_tam()` also flagged residual spreads outside its reference range, large standardized residuals, and adjacent-year residual correlation in some groups.

Over the common 1985–2023 period, mean absolute percent differences were 17.5% for N-at-age, 13.1% for total abundance, 14.0% for recruitment, and 14.5% for approximate SSB. Their trend correlations were 0.952, 0.980, 0.970 and 0.980, respectively. At 2023, tinyAM was 62.3% lower for the mean age-specific N comparison, 14.0% lower for total abundance, 4.5% lower for recruitment and 6.7% lower for approximate SSB. Total biomass had a 14.0% mean absolute percent difference and 0.984 trend correlation, but was 54.7% lower in 2023. Fixed M-at-age matches by construction; F-at-age and Fbar are unavailable from the report. The SSB comparison is approximate because age 0 is excluded. These results are a restricted comparison to published trend estimates, not a reconstruction of the unverified WBCod2026 SAM fit. The dashboard is cached locally at `results/cache/ices_cod_western_baltic_2026/dashboard.html`.

### tinyAM fit result

The default-initialized tinyAM fit passed numerical checks: optimizer code 0,
objective 671.63, maximum absolute gradient 0.000356, successful uncertainty
report, and a positive-definite Hessian (29 fixed and 648 random effects). The
fit did not use source-based starting values.

Over common years, mean absolute percent differences were 16.9% for matched
SSB, 15.5% for derived total abundance, 16.1% for total biomass, and 18.4% for
N-at-age. Their annual trend correlations were 0.98, 0.98, 0.98, and 0.96,
respectively. However, at the terminal comparison year 2024, tinyAM was 70.5%
lower for SSB, 31.1% lower for total abundance, and 58.1% lower for total
biomass. Same-year age-0 recruitment is an approximate comparison and was
26.1% lower in 2024. Fixed M matches by construction; F and Fbar are not
like-for-like with the SAM age-state definitions. Convergence therefore does
not imply close terminal agreement. The exploratory comparison is available
in the local cached dashboard at
`results/cache/ices_cod_western_baltic_2025/dashboard.html`.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ck_ices_wgbfas_morhua_bs.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The translation uses an IID recruitment process and cohort deviations,
IID F states, fixed source M, source survey timings and survey-by-age q groups. It cannot reproduce SAM's joint survey residual correlation, native catch
observation weights, or the exact initial-state and process parameterization. These
are recorded as approximations; no observations are adjusted to improve
agreement.
