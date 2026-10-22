# Cohort-covariate validation

`cohort_recovery.R` is a focused pilot, not validation for routine migration
inference. Run `cohort_gaussian_recovery()` and `cohort_recovery()` after loading
the development package from the repository root. Generated attempts, warnings,
errors and checks are saved under ignored `results/`.

The mean multiplies the survivor prediction by `exp(X beta)`, using the
**destination** year and age. The plus-group effect follows the incoming plus
retained survivor sum. Recruitment, initialization and biological-age shares
are unchanged. Coefficients are signed log-abundance adjustments. A temperature
level means a proportional annual gain/loss; it does not automatically mean a
temperature-driven change in distribution. Supply a change/lag explicitly if
that is the scientific hypothesis.

## Results

Seed 3141, 100 known-survivor Gaussian datasets: true coefficient 0.12,
mean estimate 0.1194; 96% Wald interval coverage.

Seed 2718, ten complete simulated assessments per scenario (30 years, six ages,
IID N/F, recruitment RW, two lognormal observation streams). The true cohort
coefficient is 0.12. Each model is fitted from package starts, with no rescue
fitting; all attempts are retained.

| Fitted specification | Numerical success | Median coefficient | Coverage, all attempts | Median log-N RMSE |
|---|---:|---:|---:|---:|
| Cohort covariate only | 9/10 | 0.1180 | 10/10 | 0.0899 |
| Same covariate also in M | 9/10 | 0.1745 | 8/10 | 0.1214 |
| Same covariate also in q | 8/10 | 0.1228 | 10/10 | 0.0986 |

Coverage includes numerically failed attempts with available uncertainty; it
is not a claim that those fits are usable. This small pilot does not estimate
success or coverage precisely. Mortality overlap visibly weakens attribution;
observation overlap also costs numerical stability. More information does not
guarantee that mechanisms can be separated. Do not interpret a cohort effect
as immigration, emigration or a causal environmental effect from this fit alone.

The first development run omitted the survey covariate from the index table,
so its ten survey-formula attempts correctly failed input validation. That run
was retained locally; the documented design includes the covariate in both
its population and survey locations. No scientific model was tuned in response.

## Regression checks

Default parameter lists, joint objectives and gradients are exactly unchanged
for N off, IID, AR1 and RW against the pre-change snapshot. Hand-calculated
cohort/plus equations, finite derivatives, generative simulation, unused-row
missing values, polynomial projection bases and intrinsic-RW redundancy are
covered by unit tests. A fitted cod example converged (maximum gradient 0.00109,
positive-definite Hessian); one warm-started hindcast fold and a two-year update,
process simulation and a locally rendered dashboard completed. Its hindcast
scores illustrate workflow execution, not a recommendation for this covariate.
