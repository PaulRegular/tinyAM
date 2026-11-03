# Northeast Arctic saithe source review: tinyAM translation notes

Transferred from the pre-cleanup SAIO review without refitting. Check the registered stock script and review tier for current settings. [Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_ices_afwg_virens_nea_1_2.md).

## Outputs and tinyAM translation

The native run supplies annual N-at-age and F-at-age, fitted catch and survey
predictions, survey q estimates, and annual recruitment, SSB, Fbar (ages 4–7)
and total biomass. Age-specific N/F standard errors are obtained from the
native TMB report on the log scale and converted to natural-scale standard
errors with the delta method; displayed limits use a normal approximation on
the log scale. These are model-based approximations, not intervals transcribed
from a published age-surface table. The aggregate HTML table's Low/High
columns are retained as given; their interval level is not identified in the
cached table. The 2026 forecast recruitment and SSB are taken from ICES SAG;
age-specific forecast values and forecast intervals are not available.

The translation fits ages 3–12 over 1960–2025, uses free initial older-age
abundance, IID recruitment and survival residuals, correlated random-walk F,
and fixed source M. It preserves the two surveys, their q sharing, and the
native relative precision weights (unweighted rows use the source's unit
multiplier). tinyAM estimates a separate base observation SD for each survey.
It does not reproduce the source's cross-age FLT14 residual correlation or
catch prediction-variance link. Its F states for ages 8–12 are independent
age-specific states rather than one shared source state. The published Fbar is
an arithmetic summary; the dashboard labels its comparison with tinyAM's
N-weighted Fbar as approximate. M is checked as a reused fixed input, not
treated as an estimated comparison.

The source age surfaces allow total abundance to be derived over ages 3–12.
The common-definition SSB comparison uses the same official N, weight and
maturity surfaces and the start-of-year convention implied by the zero
pre-spawning mortality fractions. Dashboard comparisons retain native
aggregates and label these common-definition summaries separately.

The first tinyAM fit used package starting values and failed with a non-finite
gradient. One documented retry used the official N-at-age and F-at-age only as
starting values; no observations or model settings changed. That fit passed
`check_tam()`: optimizer code 0, objective 616.9131, maximum adjusted gradient
0.0009991, successful standard-error calculation and a positive-definite
Hessian. It used 29 fixed effects and 1,308 random effects. The diagnostics
also flag weak estimation of the Beverton–Holt curve, one strongly correlated
fixed-parameter pair, and residual outliers/spread differences.

Across common historical years, mean absolute percent differences were 9.9%
for N-at-age, 10.3% for F-at-age, 5.7% for recruitment, 6.0% for SSB, 3.4% for
total abundance and 3.6% for total biomass. Within-age temporal correlations
were 0.988 for N and 0.929 for F; aggregate trend correlations ranged from
0.986 to 0.991. Terminal-year differences were larger: N-at-age +30.7%,
F-at-age −24.1%, recruitment +3.1%, SSB +30.3%, total abundance +7.5%, and
total biomass +17.3%. These differences are retained for review rather than
used to alter source observations or tune the translation to the official
estimates. The dashboard is stored in the local ignored results cache.
