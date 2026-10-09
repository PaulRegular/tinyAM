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

## F/M mean validation

Stage 3 was authorized on 2026-10-09. IID/RW/AR1 effects now enter the F/M mean
formulas using the same unique-state densities as q. Initially use IID residuals
(or M residuals off), and one temporal mean term. Overlapping estimated IID
variances, saturated fixed/random terms and inconsistent M age blocks are rejected.
This does not add an F residual-off option or change projection/boundary rules.

Run from the repository root:

```r
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_mean_effects.R", "100", "30"))
rmarkdown::render("analysis/formula_expansion/mean_effects_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

Track these two sources and the small deterministic package tests. Generated
results, representative fits, figures and the HTML report are all in ignored
`results/`. The report retains seeds and the result object records session
information. Routine tests verify mathematics and integration; Monte Carlo
recovery remains an explicit analysis, not a slow or flaky CI test.

### Evidence and decision

The completed study contains 600 isolated fits and 360 full assessment fits.
All isolated fits passed the numerical criteria. Full-model success ranged
from 67% to 97% by specification; the 90% target was not met uniformly.
F SD recovery was generally reasonable, but gradient thresholds were often
missed. Small shared M variation was difficult to separate from IID M residuals;
a larger shared signal improved recovery. AR1 correlation remained imprecise.

Retain one structured mean with the existing restrictions. Prefer mean-only M
before attempting to estimate a second variance component, and require
stock-specific numerical and sensitivity checks. Do not enable overlapping
temporal components or variance processes on this evidence. The report shows
success counts, parameter recovery, intervals and example trajectories.

Final package checks: 1,564 passing assertions, one interactive-only skip,
and `R CMD check --no-manual` Status OK. PDF manual generation is blocked by
existing Unicode mathematical symbols in the simulation help page.
