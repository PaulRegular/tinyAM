# NAFO translation readiness and preliminary fits

This review uses the selected local SAIO snapshot at commit `eee9dd7e2d683041e8f25d1d28e2db8d483c22b1`. Source tables and reviews remain in SAIO. The three fit recipes use ages and inputs retained in the source tables; no age data were digitized from plots or inferred from fitted outputs.

| Assessment record | Readiness and fitted grid | Main limitation |
|---|---|---|
| `nafo_cod_3m_2026` | Ready for a limited fit, 1988–2025, ages 2–8+ | Age 1 is omitted; source zeros are omitted as specified; the incomplete 2018 age-8+ survey value is excluded. Source-estimated M is represented by fixed prior medians. |
| `nafo_cod_3no_2021` | Inputs support an approximate ADAPT translation, 1960–2019, ages 3–12 | Age 2 is omitted because stock weights begin at age 3. Source catch is exact and survey residuals have equal log-scale weights; tinyAM uses a tight finite catch SD and unit index SD. Age-12 F constraints are not reproduced. |
| `nafo_american_plaice_3lno_2021` | Inputs support an exploratory translation, 1960–2020, ages 5–15+ | The 2026 record is monitoring only. The historical 2021 run has no tabulated annual N/F/SSB series; its terminal survivor parameters are not substitutes. The source warns that absolute estimates are unreliable. |
| `nafo_redfish_3m_2024` | Not fit-ready | The accessible sources provide total catch and short-term scenario outputs, but no age-specific catch, index, weight or maturity matrices. |

All three recipe grids passed the no-fit `run_assessment(..., fit = FALSE)` checks after the final settings were applied. Catch and index observation SDs are supplied as fixed values, so the recipes use `sd_form = ~0`.

## 3M cod fit and comparison

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

The final-settings 3NO default fit was stopped after about four minutes without returning an optimizer result. Its optimizer code, objective, gradient and uncertainty diagnostics are unavailable; this is an interrupted fit, not a confirmed numerical failure. An earlier attempt with the initial SD formulas failed after 40.5 seconds with NA/NaN evaluations, and its source-start retry was interrupted after more than eleven minutes. That earlier attempt used different settings and is not treated as a result for the final recipe. No final-recipe source retry was started because the final default attempt did not return convergence diagnostics.

The 3LNO 2021 default fit was stopped after about four minutes without returning a result. The accessible source output does not provide annual N/F surfaces for aligned starting values, so no source-based retry is available. The 2026 interim-monitoring record is not treated as a model assessment. No Official–tinyAM dashboard is produced because the 2021 output does not provide a comparable annual population series.

The interrupted attempts are marked in `results/preliminary/nafo_fit_diagnostics.csv` and their local initialization records. No tinyAM package code or core model mathematics was changed.
