# Formula expansion

Experimental branch: `formula-expansion`, starting at `e801800` on
`comp-assessments`. Milestones are committed separately; this branch is not merged.

## Stages and review gates

1. Centralize numerical, structural and advisory diagnostics in `check_tam()`.
   Use them in fitting, retrospectives, print/summary methods and the dashboard.
2. Expand `q_form` with IID, random intercepts, RW, AR1 and logistic selectivity.
   Validate mathematics, simulation, uncertainty, projections, warm starts and
   dashboard displays before evaluating stock applications.
3. Extend F/M mean formulas, as subsequently authorized by the user.
   Begin with shared temporal mean processes and IID residual processes.
   Multiple temporal processes need a separate identifiability review.
4. Extend catch/index SD formulas with age IID/RW/AR1 effects, as subsequently
   authorized by the user. Evaluate replicated and sparse SD recovery before
   considering temporal stochastic volatility or latent-process SD formulas.

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
with intervals. The Parameters menu has Fixed and Random pages; formula tabs
show effect trends only, with short readable names. Increments and numeric-by
contributions remain available in tidy data, but are not dashboard tabs.
Report process SDs, correlation,
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

## Catchability validation

Phase 2 recovery and stock experiments are separate from the F/M mean study.
Run from the repository root:

```r
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/test_q_validation.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_q_effects.R", "100", "30"))
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/validate_q_workflow.R")
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/evaluate_q_stocks.R")
rmarkdown::render("analysis/formula_expansion/q_effects_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

`q_validation.R` contains the shared designs and recovery calculations. The study
uses both q links, categorical and numeric multipliers, six effect specifications,
replicated and sparse observations, and estimated process/observation SDs.
Failures and warnings are retained, including structural rejections. The recovery
file records every attempt, seed, interval, and numerical check; checkpoint files
are marked incomplete until the entire study finishes.

Stock experiments fit the settled recipes afresh from a pinned committed database.
Candidate q formulas are in a separate experiment script; canonical stock scripts,
their observations and other model settings are not changed. GOA pollock candidates
replace its annual ADF&G fixed effects rather than duplicating them. Failed
candidates remain labelled in local dashboards. Aggregate comparative-assessment
CSV files are not overwritten. All fits, diagnostic tables, dashboards and report
figures go under ignored `results/`; retain only code, tests and report text in Git.

The recovery report distinguishes mathematical correctness, numerical success,
parameter/curve recovery, and agreement with accepted assessment outputs. Neither
good numerical checks nor closer accepted trajectories establish identifiability.
Retrospective tests report failed folds and never treat projection years as
assessment terminal years. Stock decisions need review before replacing baselines.

### Phase 2 evidence

Replicated full-model numerical success was 73–100%, below the 90% target in
several cases. Isolated Gaussian fits all passed with replicated observations;
logistic fits were less reliable. Sparse numeric-by effects could converge while
underestimating their SD, and AR1 correlation tended to be underestimated.
Keep these options experimental and inspect uncertainty, rather than treating
convergence as evidence that changing q is identifiable.

All five stock baselines converged; ten of fourteen candidates converged.
EBS RW/AR1, Southern Gulf RW and herring logistic did not pass. Baseline stock
recipes remain unchanged. One logistic retrospective fold also failed; its actual
retained fold count is recorded. The HTML report includes recovery figures,
intervals, diagnostic tables and definition-matched assessment comparisons.

Phase 2 checks: 1,603 passing package assertions, one interactive-only skip,
all 29 comparative-assessment test files passed, and `R CMD check --no-manual`
Status OK. Recovery comprises 2,400 isolated and 480 full-model attempts.

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

## Observation SD validation

Phase 4 adds Gaussian formula effects to `catch_settings$sd_form` and
`index_settings$sd_form`. Effects enter log SD, not the observation mean.
Supplied SDs remain offsets; process SDs quantify variation in log observation
SD. Begin with replicated age effects; a small number of ages may poorly
identify AR1 correlation even with many observations per age. One RW/AR1
term per SD surface is currently permitted. This does not add formulas for
latent N/F/M SDs or validate temporal stochastic volatility.

Run from the repository root:

```r
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/test_sd_validation.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_sd_effects.R", "100", "30", "3"))
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/validate_sd_workflow.R")
rmarkdown::render("analysis/formula_expansion/sd_effects_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

`sd_validation.R` contains the designs and recovery calculations. The study
records 900 isolated process-recovery attempts, 500 noisy-tail fits comparing
common/quadratic/IID/RW/AR1 curves, and 180 full assessments varying catch or
survey age SDs. Every attempt, seed, warning, interval and error is retained.
Numerical success and SD-curve recovery are reported separately; coverage is
conditional on numerical success and available intervals. No failed fit is
rescued or silently removed. The optional third argument controls independent
R workers (default one). Checkpoints resume missing attempts only; fixed seeds
inside each replicate make results independent of scheduling. Use the same
replicate counts when resuming. Workflow checks cover named warm starts,
projections, retrospective/hindcast folds, simulation, tidy effects and two
dashboards. Only the source, concise report text and deterministic tests are
tracked; results, representative fits, figures and dashboards are ignored.
Source fingerprints prevent resuming a checkpoint with changed model or study
code; move the old checkpoint aside to start a new study after source edits.

### Phase 4 evidence

All 1,580 attempts are complete. Full assessments passed numerical checks in
30/30 catch IID, 13/30 catch RW, 30/30 catch AR1, 30/30 survey IID, 28/30 survey
RW and 29/30 survey AR1 cases. All well-replicated isolated processes passed,
but sparse SD curves were poorly recovered despite frequent numerical success.
AR1 correlations were imprecise and underestimated. Catch process SDs also
tended to be underestimated. Retain the options experimentally, prefer IID
age effects as a starting point, and flag estimated catch-SD age RWs.

The correctly specified quadratic curve outperformed random effects for the
quadratic noisy-tail example; the random curves nevertheless captured the
tails. A recovery-curve synchronization bug after `sdreport()` was corrected
and covered by an intercept-score regression test. Standalone experiments
were recomputed with identical seeds; full assessment fits were unaffected.
Both workflow examples passed projections, two retrospective folds, two
hindcast folds, simulation and dashboard checks. The final full package suite
and `R CMD check --no-manual` passed.

## Recruitment validation

Recruitment formulas use the existing absolute log-recruitment states. They
support a default RW, IID/AR1 residuals, annual covariates, and Beverton–Holt or
Ricker curves with IID/AR1 residuals. The curve describes median recruitment.
Parent SSB is start-of-year; unavailable pre-model parents leave free fixed
boundary recruitment states rather than reconstructed SSB.

Run from the repository root:

```r
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/test_recruitment_validation.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_recruitment.R", "100", "30", "3"))
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/validate_recruitment_workflow.R")
rmarkdown::render("analysis/formula_expansion/recruitment_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

`recruitment_validation.R` holds the designs and recovery calculations. The
study records 1,400 isolated attempts and 180 full assessments across six
specifications, with narrow-SSB and correlated-covariate stress cases. Seeds,
warnings, failures, uncertainty, and source fingerprints are retained. The
optional worker count changes scheduling only; checkpoints resume missing
attempts with unchanged source and replicate counts. Numerical success,
parameter recovery, and coverage are separate outcomes.

The workflow checks projections, named starts, retrospective/hindcast folds,
simulation, tidy uncertainty, and three dashboards. The runnable package
example is `inst/examples/example_recruitment.R`. Keep these sources and the
concise report text; generated results, fits, figures, and HTML remain ignored.

A study-only survey-row assignment error invalidated an earlier set of full
attempts. Those attempts and the original study source remain in the ignored
`results/recruitment_invalid_observation_order/` archive. The full experiment
was restarted with corrected row alignment; a deterministic regression test
checks that simulated observations survive preparation without reordering.
Package simulation and the earlier formula studies were unaffected.

### Recruitment evidence

All 1,400 isolated attempts passed numerical checks. Full assessments passed
30/30 BH–IID, 29/30 BH–AR1, 28/30 Ricker–IID, 27/30 Ricker–AR1, 24/30
covariate–IID and 28/30 covariate–RW attempts. Four covariate–IID failures came
from the prespecified RW warm start; two had gradients just above 0.01. Keep
these options experimental: curve support and interval calibration remain
important limitations even when fitting succeeds. The full suite passed 1,818
assertions with one interactive-only skip; `R CMD check --no-manual` was OK.

## Recruitment self-tests, cross-tests and noise

The follow-up crosses RW, BH–IID, BH–AR1, Ricker–IID and Ricker–AR1 truths
with all five fitted recruitment specifications. Each dataset is shared across
its five fits. This explicitly tests fitting curves to RW-generated recruitment,
as well as recovering the correct curve. The original recovery study fits only
the generating specification.

```r
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/test_recruitment_cross_validation.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_recruitment_cross.R", "100", "30", "6"))
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_recruitment_noise.R", "100", "30", "3"))
rmarkdown::render("analysis/formula_expansion/recruitment_cross_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

The cross-test records 2,500 isolated and 750 full requested fits. The BH noise
study adds 600 isolated and 180 full fits at SDs 0.15, 0.35 and 0.60. It separates
narrow parent SSB from wider contrast; full wider-SSB cases use a declared
synthetic fishing-pressure covariate in both generating and fitted F means.
No assessment recipe is changed. These checkpoints also retain failures,
warnings, seeds and source fingerprints. Curve parameters have no recovery
truth when recruitment is a RW. Diagnostic flags are descriptive cautions,
not a formal test for a biological stock–recruit relationship.
