# Validate sampled-state quantile residuals

In 80 correctly specified tinyAM fits, drawing the random effects once restored nearly the expected spread and frequency of extreme residuals. Evaluating the residuals at the fitted random-effect modes pulled them toward zero, especially for observations that strongly informed their own fitted states.

Run from the package root with R and the tinyAM development dependencies available.

```sh
Rscript --vanilla analysis/formula_expansion/validate_sampled_quantile_residuals.R
```

The script generates 40 independent 32-year histories with ages 1–4+, constant recruitment and mortality, and IID deviations from predicted survivors. Each history supplies the same catch and stock to two scenarios. The survey observes the population through either a lognormal or a Poisson-link delta lognormal distribution. The fits estimate fishing mortality, recruitment and the cohort-process SD. They use known natural mortality, survey catchability, observation SDs and the PLD detection scale to focus this check on the diagnostic. Each fit supplies one paired comparison of fitted-state (`mle-eb`) and sampled-state (`mle-mvn`) quantile residuals. Both use the same random numbers for PLD zeros.

All 40 pairs completed in each scenario without warnings. The table averages within-history SD and tail frequency across the 40 histories. For a calibrated standard normal diagnostic, SD is one and 5% of residuals fall outside ±1.96.

| Survey family | Observation | SD, fitted state | SD, sampled state | Outside ±1.96, fitted state | Outside ±1.96, sampled state |
|---|---|---:|---:|---:|---:|
| Lognormal | Catch | 0.859 | 0.987 | 2.66% | 4.65% |
| Lognormal | Index | 0.802 | 1.001 | 1.78% | 5.23% |
| PLD lognormal | Catch | 0.747 | 0.982 | 1.52% | 4.77% |
| PLD lognormal | Index | 0.944 | 0.980 | 3.91% | 4.43% |

Paired Monte Carlo SEs for the change in SD were 0.004–0.008. The average sampled-state residual was within 0.019 of zero in all four comparisons. The script writes all attempts, per-history scores, average calibration and paired differences to the ignored results directory.

The simulation gives a calibration check under these simple, correctly specified operating models. It leaves catchability, detection and observation dispersion fixed, and uses independent observation errors. In a fit with age-AR(1) observation errors, marginal quantile residuals can still be correlated within series and year. The conditional draw relies on a Gaussian approximation around the fitted state and uses the assessment’s own observations. Hindcasts and one-step-ahead residuals answer different predictive questions.
