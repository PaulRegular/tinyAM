# Formula expansion

Experimental branch: `formula-expansion`, starting at `e801800` on
`comp-assessments`. Milestones are committed separately; this branch is not merged.

## Stages and review gates

1. Centralize numerical, structural and advisory diagnostics in `check_tam()`.
   Use them in fitting, retrospectives, print/summary methods and the dashboard.
2. Expand `q_form` with IID, random intercepts, RW, AR1 and logistic selectivity.
   Validate mathematics, simulation, uncertainty, projections, warm starts and
   dashboard displays before evaluating stock applications.
3. Review the Stage 2 evidence with the user before extending F/M mean formulas.
   Begin with shared temporal mean processes and IID residual processes.
   Multiple temporal processes need a separate identifiability review.
4. Review possible variance structures as a research question. Do not expose
   stochastic volatility or another variance process without a further decision.

## Agreed conventions

- Gaussian formula effects enter the selected log/logit predictor.
- Numeric `by` multiplies one shared process; categorical `by` gives separate
  group trajectories sharing a term-level SD and, for AR1, correlation.
- SDs are estimated by default. IID SD is marginal; RW/AR1 SD is innovation SD.
- RWs have a zero first effect; ordinary formula terms supply the baseline.
  Densities contain only unique increments. Calendar gaps retain missing steps.
- AR1 effects are stationary, with positive estimated correlation. A supplied
  correlation may be zero. Proper IID/AR1 effects are not arbitrarily centered.
- Logistic selectivity is a multiplicative rising curve with midpoint and
  positive slope. Categorical `by` gives separate curves. Numeric `by` is not
  supported initially. Ordinary terms supply baseline q; logit q remains below 1.
- Preserve ordinary formulas and mono plateaus. Avoid redundant parameters and
  reject exact variance aliases, rather than imposing hidden biological constraints.
- Future point predictions follow the process: RW holds its last effect, AR1
  returns toward zero, and a new IID level has zero effect. Process simulations
  include future variability. These are conditional link-scale predictions.
- Warm starts match stable term and state names. One registry controls random
  parameters, simulation and reporting.

## Reporting and dashboard

Show signed IID/random-intercept effects with intervals and RW/AR1 trajectories
with intervals. Show RW increments separately. Distinguish numeric-by raw
effects from their multiplied contributions. Report process SDs, correlation,
and logistic midpoint/slope with explicit SE scales. Catchability pages show
the combined response-scale curve and uncertainty, observed support and
projections. Retain existing visual style and framed ribbons; allow incomplete
assessment references alongside fitted models.

## Validation

Each milestone needs deterministic unit tests and relevant integration tests.
Compare densities with exact normal calculations (about 1e-10 tolerance),
check derivatives, and exercise simulation, projections, updates and retrospectives.
Run the full suite and R CMD check for substantial changes.

Recovery experiments use reproducible seeds: 100 standalone well-informed
replicates and 30 full-model replicates per representative specification.
Aim for at least 90% numerical success in well-informed experiments; interpret
bias and interval coverage against Monte Carlo uncertainty. Sparse/confounded
cases must be rejected or clearly flagged without convergence rescue hacks.

Evaluate candidates separately from the accepted models for GOA cod, EBS and
GOA pollock, Southern Gulf cod and spring herring. Change one assumption at a
time and preserve their data and accepted baseline settings. Agreement with
accepted trajectories alone is not evidence of statistical correctness.
