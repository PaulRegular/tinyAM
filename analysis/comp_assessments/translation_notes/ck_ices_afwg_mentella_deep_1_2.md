# Northeast Arctic beaked redfish (`Sebastes mentella`): tinyAM translation notes

Transferred from the pre-cleanup SAIO review without refitting. Check the registered stock script and review tier for current settings. [Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_ices_afwg_mentella_deep_1_2.md).

## Limited tinyAM translation

The maintained recipe is
[preliminary recipe](../scripts/translation/preliminary/jrn_redfish_beaked_2026.R). It fits 1992–2025
at ages 6–19+, using the non-overlapping total catch-at-age series once,
fixed M = 0.05, published age-6+ weight and maturity, and the ecosystem and
winter number indices for ages 6–15. The million-fish indices are converted
to thousand fish by the shared observation converter. It does not alter
observations using fitted q. The source's age-2 recruitment is not represented
or renamed as tinyAM recruitment; entry at age 6 is an age-range boundary.

The source framework says the ecosystem and winter survey selectivity is high
for young fish and declines at older ages. The tinyAM recipe estimates a
separate intercept and log-linear age effect for each survey. This allows a
decline and is a compact approximation to the source's survey-specific
selectivity curves. The source fleet structure, selectivity processes,
catch-weight likelihood, and observation weights are not represented.

The default-start fit converged and its Hessian was positive definite. It had
41 fixed effects, 509 random effects, an objective of 1031.742, and a maximum
absolute gradient of 0.00134. The optimizer returned code 0 with relative
convergence (4), and uncertainty estimation succeeded. Residual checks still
flagged some age/series groups with standardized residual spread outside
0.5–1.5 and some residuals above three standard deviations. Numerical
convergence is not evidence of close replication.

| Comparison | Status | Median absolute percent difference | Terminal mean percent difference | Trend correlation |
|---|---|---:|---:|---:|
| N at ages 6–19+ | Approximate | 31.6% | -44.9% | 0.926 |
| F at ages 6–19+ | Approximate | 55.8% | +85.9% | 0.819 |
| Fbar, ages 12–18 | Approximate | 38.4% | +90.5% | 0.790 |
| SSB, common definition at ages 6–19+ | Approximate | 46.4% | -50.7% | 0.983 |

The common-definition source SSB is calculated without uncertainty from
reported age-6+ N, fixed weight and maturity. The source's native SSB remains
available in the dashboard. Fixed M is shown as an input-consistency check;
its agreement is not evidence of estimated-M performance. Recruitment and
whole-stock abundance or biomass are non-equivalent because tinyAM starts at
age 6. Age-specific biomass is unavailable from the source outputs.

The dashboard and fitted object are cached locally under
`results/cache/jrn_redfish_beaked_2026/`; they are generated outputs and are
not committed.
