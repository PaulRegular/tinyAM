# Information-criterion validation

Run `Rscript analysis/information_criteria/validate.R` from the repository root.
The shorter `--quick` run checks execution only. Generated CSVs and session
information stay in the ignored `results/` directory.

The Gaussian random-intercept experiment uses 100 replicates. Its sparse trace
is checked against an independent dense calculation. The two tinyAM experiments
use 30 replicates each: Beverton–Holt recruitment and active IID cohort survival.
These deliberately exercise nonlinear latent-state distributions, where the
Gaussian-random-effect assumption behind Equation 6 is not guaranteed.

Each experiment compares `2 * EDF` with prediction optimism: the expected
deviance of new observations sharing the *same realized latent population*,
minus the fitted training deviance. The expectation is analytic for Gaussian
log observations. This is not a hindcast or new-population forecast experiment.
The same generating parameter values initialize every assessment replicate;
there is no retry or rescue fitting. Every warning and error is retained.
Numerical convergence and criterion calculation are summarized separately.
Monte Carlo uncertainty is reported for the mean penalty-minus-optimism difference.

## Method conventions

The implementation follows Equation 6 of Zheng, Cadigan and Thorson (2026),
[doi:10.1080/00949655.2026.2736799](https://doi.org/10.1080/00949655.2026.2736799).
The accompanying arXiv source archive, version 2, contains
`R_TMB_codes_for_Gaussian_example.pdf`. Its `cAIC_calc()` uses `p - 1`, excluding
a prior-only correlation parameter from the conditional fixed-parameter penalty.
The example is noncentred: some variance parameters enter observation predictions.
tinyAM counts structural dependencies in its own fitted parameterization.

Reference implementations reviewed on 19 October 2026:

- [sdmTMB `R/caic.R`](https://github.com/sdmTMB/sdmTMB/blob/main/R/caic.R):
  process-only objective, sparse random Hessians and selected inverse.
- [WHAM `R/aic.R`](https://github.com/timjmiller/wham/blob/devel/R/aic.R):
  process-only objective and dense inverse for the trace.

Both currently count all optimized outer parameters. tinyAM distinguishes this
marginal count (`n_fixed`) from the conditional observation-likelihood count
(`edf_fixed`). Process-only parameters are excluded from the latter. EDF
contributions follow the fitted latent coordinates; changing to a parameterization
that moves variance parameters into the observation model can change this
approximation. These are not claimed to be exact finite-sample criteria.

## Validation results

The seeded run used RTMB 2.0 and R 4.6.1. All 160 attempts produced a criterion;
the table below reports numerical success separately and summarizes optimism
only among numerically converged fits. The cached summary also includes every
attempt, including failed fits.

| Experiment | Numerical success | Mean `2 * EDF` | Mean prediction optimism | Difference (Monte Carlo SE) |
|---|---:|---:|---:|---:|
| Gaussian random intercept | 100/100 | 25.03 | 23.37 | 1.66 (1.58) |
| Beverton–Holt recruitment | 28/30 | 91.78 | 119.54 | -27.76 (15.89) |
| IID cohort survival | 24/30 | 143.47 | 183.54 | -40.06 (28.67) |

Gaussian sparse-versus-dense trace contributions differed by at most
`2.8e-17`. The Gaussian mean penalty agrees with measured optimism within Monte
Carlo uncertainty. In the nonlinear experiments, the estimated penalties were
lower on average, but uncertainty was substantial. These runs do not establish
calibration for nonlinear latent states.

Failed fits had extreme prediction losses and dominated the all-attempt mean
optimism (about `4.7e9` for BH and `3.4e12` for cohort survival). They remain in
the raw results; they were not retried, tuned or silently discarded. All
nonlinear calls produced the intended applicability warning, and failed
numerical checks produced an additional warning. A computable criterion does
not make a failed fit suitable for model selection.

## Package verification

The targeted criterion tests passed 64 assertions. The final package check
passed all 2,922 assertions, with no failures, warnings or skips;
`R CMD check --no-manual` returned `Status: OK` on R 4.6.1.
Comparing the likelihood against `main` at `8b4730c` gave exactly zero
differences in objectives, gradients and N/F/M/SSB/observation reports for
representative N-IID/F-RW and M-AR1/F-IID configurations.
