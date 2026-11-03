# Scup (US Atlantic): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_nefsc_saw_chrysops_usatl.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_nefsc_saw_chrysops_usatl.md).

## Expanded-library review

### Translation results and remaining limits

The translation fits 1984-2014 and ages 0-7+, using the four source catch-at-age series, fixed M, annual biology, and the NEFSC fall age-index table. It uses tinyAM's normal initialization with no explicit start values. The model uses deterministic exponential N survival, IID recruitment, and AR1 F. The age-index table is not likelihood-equivalent to the official age-composition plus aggregate-index likelihood; tinyAM represents its positive age-specific values as lognormal indices and omits zero cells. The 7+ weight is calculated by weighting source ages 7-10 in Table A24 by the four source catch-at-age series. The missing 1999 plus weight is interpolated between adjacent years for the fit only.

The fit passed numerical convergence checks: optimizer code 0, adjusted and raw maximum gradient 0.000222, positive-definite Hessian, and successful standard-error report. It has 24 fixed effects and 278 random effects, with objective 559.463. `check_tam()` also reports residual spread outside its advisory range for some groups and standardized residuals over three SDs. The dashboard is cached locally at `results/cache/nefsc_scup_2015/dashboard.html`.

The common-year comparison remains uneven. Mean absolute percent differences are 312.8% for F-at-age, 286.1% for Fbar, 51.1% for N-at-age, 43.3% for total abundance, and 38.6% for recruitment. Common-definition biomass and SSB comparisons are approximate and differ by 50.2% and 56.1%, respectively. Their terminal differences are -77.9% for biomass and -82.5% for SSB. Trend correlations range from 0.73 for recruitment to 0.95 for SSB; F-at-age correlation is 0.88. These results show a numerically usable but not close reproduction. Fixed M matches by construction because the same 0.20 input is supplied; it is not evidence of estimated agreement.

For the dashboard's common-definition comparison, source biomass-at-age is derived from official N-at-age and translated annual weights. Source SSB-at-age is derived from official N, F, and M with the same translated weights and maturity, using the June midpoint at 0.5 of the year. The directly reported aggregate SSB and its available uncertainty remain in the canonical source record. These derived comparison values are analysis-only and carry no invented uncertainty. The five omitted surveys, other active survey series and CV adjustments, and complete ASAP control settings remain unresolved, so the 2015 record stays partial and the fit remains a one-survey approximation to S60_BASE_18.
