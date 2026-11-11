# NAFO translation readiness and fits

The original review used the local SAIO snapshot at commit `eee9dd7e2d683041e8f25d1d28e2db8d483c22b1`. Source tables and reviews remain in SAIO. The three fit candidate recipes use ages and inputs retained in the source tables; a fourth preliminary recipe records a source blocker. No age data were digitized from plots or inferred from fitted outputs.

| Assessment record | Readiness and fitted grid | Main limitation |
|---|---|---|
| `nafo_cod_3m_2026` | Ready for a limited fit, 1988–2025, ages 2–8+ | Age 1 is omitted; source zeros are omitted as specified; the incomplete 2018 age-8+ survey value is excluded. Source-estimated M is represented by fixed prior medians. |
| `nafo_cod_3no_2021` | Inputs support an approximate ADAPT translation, 1960–2019, ages 3–12 | Age 2 is omitted because stock weights begin at age 3. Source catch is exact and survey residuals have equal log-scale weights; tinyAM uses a tight finite catch SD and unit index SD. Age-12 F constraints are not reproduced. |
| `nafo_american_plaice_3lno_2021` | Inputs support an exploratory translation, 1960–2020, ages 5–15+ | The 2026 record is monitoring only. The historical 2021 run has no tabulated annual N/F/SSB series; its terminal survivor parameters are not substitutes. The source warns that absolute estimates are unreliable. |
| `nafo_redfish_3m_2024` | Explicitly blocked before fitting | The detailed 2024 assessment report is inaccessible; available sources provide total catch and projections, but no age-specific catch, index, weight or maturity matrices. |

All three fit candidate grids passed the no-fit `run_assessment(..., fit = FALSE)` checks after the final settings were applied. The redfish recipe returns a source-readiness blocker without constructing data or fitting. Catch and index observation SDs are supplied as fixed values, so the fit recipes use `sd_form = ~0`.

## Curated Flemish Cap cod

The curated 2026 recipe uses ages 1–8+ with structural-zero fishing mortality at age 1 and a shared F state at ages 7–8+. The EU survey's age-1 abundance exponent is estimated with `q_power()`. Other ages retain proportional catchability. The translation represents the survey interval by its midpoint, retains the baseline 2018 survey omission and fixes M at the source prior medians. It approximates the source fleet mortality process and omits the Bayesian priors.

Both fits use SAIO `cecb45ad` and cover 1988–2025. Default starts for the age-1 trial exceeded the gradient tolerance. The recipe now uses one `fit_tam()` call with the successful source N and positive F starting values. Its maximum gradient was 0.00554 and its Hessian was positive definite. The age-2 baseline passed from default starts.

| Metric | Baseline median absolute difference | Curated median absolute difference | Curated terminal difference | Curated trend correlation |
|---|---:|---:|---:|---:|
| SSB | 15.3% | 13.3% | +12.0% | 0.995 |
| Recruitment | 21.9% | 10.3% | +6.9% | 0.983 |
| Abundance | 12.6% | 11.3% | +7.6% | 0.986 |
| Fbar | 18.0% | 13.1% | −17.6% | 0.889 |

The curated recruitment comparison uses age 1, while the baseline uses age 2 as a proxy. SSB and abundance use all modeled ages in each fit, with both reference and tinyAM series aligned to the same ages within each comparison. Fbar remains approximate because the official average is unweighted.

## Historical 3M cod fit and comparison

The final fit used package-default starting values and converged. The supplied-SD formulas were changed from `~1` to `~0` so the supplied values are not accompanied by an additional estimated scale. This changed the fitted parameter count and objective; the values below and the dashboard use the final settings.

- Optimizer code 0, relative convergence; objective 932.982.
- Maximum absolute gradient 0.00213; SD report succeeded and Hessian was positive definite.
- 20 fixed effects, 265 random effects; fit time about 12 seconds.
- Dashboard: [NAFO 3M cod preliminary comparison](../results/preliminary/cache/eee9dd7e2d683041e8f25d1d28e2db8d483c22b1/nafo_cod_3m_2026/dashboard.html).

The comparison covers 38 common years (1988–2025). Percent differences and correlations are descriptive; age-specific mean percent differences can be inflated when official values are near zero.

| Metric | Comparison | Mean absolute % difference | Median absolute % difference | Terminal difference | Trend correlation |
|---|---|---:|---:|---:|---:|
| F at age | matched ages 2–8+ | 1,057.2 | 24.7 | −30.5% | 0.861 |
| Fbar | approximate: both ages 3–5, but tinyAM weights by abundance | 286.4 | 18.0 | −51.7% | 0.787 |
| N at age | matched ages 2–8+ | 28.2 | 20.2 | +47.7% | 0.969 |
| Total abundance | approximate common age range 2–8+ | 28.9 | 12.6 | +59.1% | 0.964 |
| Biomass | approximate common age range 2–8+ | 19.9 | 11.3 | +51.7% | 0.962 |
| Recruitment | approximate: tinyAM age-2 N vs source age-1 recruitment | 43.4 | 21.9 | −5.3% | 0.954 |
| SSB | approximate mature biomass over ages 2–8+ | 18.4 | 15.3 | +33.8% | 0.986 |

M is unavailable as an estimated-M comparison because tinyAM fixes the source prior-median age values. The source report's SSB intervals are central 80% intervals; age-specific N and F outputs have no tabulated intervals, so their uncertainty is not directly comparable with tinyAM's SD report.

## 3NO cod and 3LNO American plaice fit attempts

The final-settings 3NO default fit returned after 599.7 seconds but failed numerical checks: optimizer code 1 (`false convergence (8)`), objective 1621.154, maximum gradient 2162.5, and non-positive-definite Hessian. An SD report object was returned, but the fixed-parameter covariance was not finite, so uncertainty checks failed. The policy-permitted source-start retry was then attempted and interrupted after approximately 30 minutes without returning a fit or retry diagnostics. An earlier attempt with the initial SD formulas also failed with NA/NaN evaluations; those results are retained separately and are not substituted for the final-settings diagnostics.

The historical 3LNO 2021 default fit converged after 527.9 seconds (optimizer code 0, objective 1720.647, maximum gradient 0.00818, 43 fixed effects and 731 random effects; SD report succeeded and the Hessian was positive definite). `check_tam()` flagged residual spread outside 0.5–1.5 for some age/series groups and adjacent-year residual correlation above 0.5. No comparable annual reference rows are available: the accessible numerical outputs are terminal survivor estimates, not annual N/F/SSB. The fit is retained for structural exploration, but there are no annual agreement statistics or Official–tinyAM dashboard. The source has no aligned annual N/F surfaces for a source-start retry; none was needed because the default fit converged. The 2026 interim-monitoring record remains separate and is not treated as a model assessment.

The final fit outcomes are recorded in `results/preliminary/nafo_fit_diagnostics.csv` and their local initialization records. The targeted `test_expansion_01_translations.R` checks passed with 153 expectations, including the redfish blocker. No tinyAM package code or core model mathematics was changed.
