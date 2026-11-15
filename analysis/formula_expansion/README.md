# Statistical validation studies

These reproducible studies examine the package's formula and process options.
They complement deterministic package tests by separating numerical success,
parameter/curve recovery and interval coverage. They are not model-selection
rules or promises of agreement with an official assessment.

## Files and execution

Keep one generating/recovery implementation and runner per study, its meaningful
smoke checks, and the report source for its distinct results. Current stock
applications belong in `analysis/comp_assessments`; superseded stock experiments
and development plans remain in Git history. This directory is retained for
statistical validation, not as an archive of assessment tuning.

Run from the repository root. `simulate_*` scripts accept isolated and full
replicate counts; supported optional worker counts change scheduling only.
Full studies are explicit analyses, not routine CI tests.

| Study | Runner / generating code | Report source |
|---|---|---|
| Catchability effects | `simulate_q_effects.R`, `q_validation.R` | `q_effects_report.Rmd` |
| F/M mean effects | `simulate_mean_effects.R` | `mean_effects_report.Rmd` |
| Observation SD effects | `simulate_sd_effects.R`, `sd_validation.R` | `sd_effects_report.Rmd` |
| RV observation covariance | `benchmark_rv_covariance.R` | [RV covariance benchmark](rv_covariance_benchmark.md) |
| Poisson-link delta and mortality recovery | `simulate_pld.R`, `pld_validation.R`, `pld_validation_fit.R` | [PLD validation report](pld_validation_report.Rmd) |
| Recruitment recovery | `simulate_recruitment.R`, `recruitment_validation.R` | `recruitment_report.Rmd` |
| Recruitment self-/cross-tests | `simulate_recruitment_cross.R`, `recruitment_cross_validation.R` | `recruitment_cross_report.Rmd` |
| BH noise and parent-SSB contrast | `simulate_recruitment_noise.R` | `recruitment_cross_report.Rmd` |
| Rising F means | `f_mean_curve_validation.R` | `f_mean_curve_report.Rmd` |
| Correlated F and process-SD sharing | `validate_process_expansion.R` | `process_expansion_report.Rmd` |
| Fleet curves, processes and smooths | `multi_fleet_recovery.R`, `multi_fleet_validation.R` | `multi_fleet_validation.md` |

For example:

```r
system2(file.path(R.home("bin"), "Rscript"),
        "analysis/formula_expansion/test_q_validation.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("analysis/formula_expansion/simulate_q_effects.R", "100", "30"))
rmarkdown::render("analysis/formula_expansion/q_effects_report.Rmd",
                  output_dir = "analysis/formula_expansion/results")
```

The recruitment cross-report requires both the cross-test and noise studies.
Run `run_f_curve_validation()` from `f_mean_curve_validation.R` for the rising-F
study; its defaults are 100 isolated, 30 full and 20 stress replicates.
`validate_process_expansion.R` defaults to 100 isolated and 30 full replicates.

`test_*_validation.R` checks the generating mathematics and observation alignment.
`validate_q_workflow.R`, `validate_sd_workflow.R` and
`validate_recruitment_workflow.R` check simulation, starts, projections, folds,
tidying and dashboard integration against their study caches. They do not tune
stock recipes. Run small smoke cases after changes; rerun the long studies only
when their mathematics or generating designs change.

## Do sampled zeros improve stock estimates?

The [controlled PLD study](pld_simple_report.Rmd) uses known survey scale,
three fitted fishing phases, two recruitment levels and a true 14+ group.
It isolates encounter information with fixed detection under exact PLD,
then compares ordinary lognormal and fitted PLD with finite age samples.
The last comparison estimates early and late natural mortality in a rise
and a constant-rate control. Each stage starts with eight histories and
expands when at least seven paired fits pass the numerical checks.

```sh
Rscript analysis/formula_expansion/simulate_pld_simple.R run analysis/formula_expansion/results/pld-simple-refined newton
Rscript analysis/formula_expansion/simulate_pld_simple.R report analysis/formula_expansion/results/pld-simple-refined analysis/formula_expansion/results/pld-simple/summary.rds analysis/formula_expansion/results/pld-simple-tight/summary.rds
```

Supply a fresh directory as the second argument to run against current package
source. On Windows use one below-normal R worker. To continue or replay the
original package implementation, set
`TINYAM_PLD_SIMPLE_SOURCE` to an extracted checkout of `5a7f6a0d`. The run
records source hashes, all fit attempts and independent calibration.

The default optimizer pilot stopped the finite-sample and mortality
comparisons at eight histories. All numerical failures exceeded the gradient
threshold. A separate pilot with relative optimizer tolerance `1e-13` produced
singular-convergence codes. The final comparison uses one damped Newton step
after the usual optimizer, accepted only when the objective and maximum
gradient decrease. All runs used the same acceptance checks.

```sh
Rscript analysis/formula_expansion/simulate_pld_simple.R run analysis/formula_expansion/results/pld-simple default
Rscript analysis/formula_expansion/simulate_pld_simple.R run analysis/formula_expansion/results/pld-simple-tight tight
```

[Controlled study report](results/pld-simple-refined/pld_simple_report.html)
[Default optimizer report](results/pld-simple/pld_simple_report.html)

## Poisson-link delta and mortality recovery under spatial monitoring

Use the SimSurvey source at [29a1b03](https://github.com/PaulRegular/SimSurvey/commit/29a1b03350c274f88dac108fb95f9446eac79233)
or [bb4a06d](https://github.com/PaulRegular/SimSurvey/commit/bb4a06d8fef8f4eb74ecfa179d6dc831f589e5fe),
which has the same simulation code. Point `TINYAM_PLD_SIMSURVEY` at the checkout
or its extracted source archive. Run these commands from the tinyAM root.

```sh
Rscript analysis/formula_expansion/simulate_pld.R calibration 6 40 analysis/formula_expansion/results/pld-30yr-calibration
Rscript analysis/formula_expansion/simulate_pld.R sampling 24 40 analysis/formula_expansion/results/pld-30yr-simsurvey
Rscript analysis/formula_expansion/simulate_pld.R fit 24 40 analysis/formula_expansion/results/pld-30yr-simsurvey
Rscript analysis/formula_expansion/simulate_pld.R report 24 40 analysis/formula_expansion/results/pld-30yr-simsurvey
```

Use one replicate to check the workflow before launching the study. An
optional fourth argument selects a different ignored results directory.
Stages run sequentially and retain one checkpoint per generating and
fitting attempt. The data manifest records the simulation sources. The
fit manifest allows reuse of generated inputs with explicit provenance
when the fitter changes.

The 30-year histories combine an early fishing pulse with a later increase
in natural mortality. The three monitoring effort levels share population
histories and realized removals. Use the Northern cod comparison in the
report to check zero frequencies at ages 2–14+ under the retained monitoring
settings. Six separate standard-effort histories estimate each program's
half-detection scale from occurrence and known observation targets for the
fixed-scale comparison. They use population identifiers 101–106, separate
from the main study's 1–24, and keep an individual input hash in the fit
manifest. Set `TINYAM_PLD_CALIBRATION` to their directory if using a different
location. Keep these outputs separate from the original 20-year experiment.

The optional `sampling_fit` stage generates and fits one SimSurvey history
at a time. Use
`TINYAM_PLD_REPORT_DIR` to refresh a report after each history. An optional
`controlled` run exercises exact-PLD and multinomial generation. Use the
same `fit` stage in a separate directory to fit those simulations.
Set `TINYAM_PLD_FISHERY_CACHE` to reuse verified removals for the same
population and F schedule when changing only observation settings. On Windows,
`run_pld_study.ps1` can start each stage with one below-normal R worker and
records CPU time and working memory. The launcher inherits the caller's R
library and Pandoc settings.

The main run uses a private, pinned SimSurvey adapter that updates cell-by-age
counts after each fishing operation. It skips terminal allocation of remaining
deaths by cell and length because subsequent age–year abundance and mortality
come from the supplied population. Monitoring is drawn separately from the
same retained removals. On six calibration snapshots spanning 5–183 operations,
the adapter matched raw operations, age–length catches, mortality references,
biological samples, and monitored estimates exactly for the same seeds. The
calibration histories retain their original-engine manifests.

## Gillnet and RV observation error

[**SimSurvey observation properties**](observation_error_report.Rmd) examines
how repeated monitoring of a fixed cod history changes age-specific precision,
zero frequency, and correlations among ages. It varies spatial patterns and
monitoring effort. The seed 21 results also supply independent calibration
for the Cadigan-style survey SD tested in the pilot.
[Generated sampling report](results/observation_error_seed21/study-grid40/observation_error_report.html)

[**tinyAM population pilot**](observation_error_pilot_report.Rmd) compares
11 SD models across 24 stable, declining and rise-and-fall stock histories.
It checks fitted SD against independent monitoring draws and uses paired
comparisons to judge changes in spawning biomass error, interval coverage
and fit reliability.
[Generated pilot report](results/observation_error_population_pilot/study-grid40/extension/observation_error_pilot_report.html)

### Render from saved results

Run from the package root with PowerShell 7.4 or later, using R with
`rmarkdown` and Pandoc available. Set `$rscript` to the Rscript executable.

```powershell
$runner = './analysis/formula_expansion/run_observation_error.ps1'
& $runner -Phase SamplingReport -Rscript $rscript
& $runner -Phase Report -Rscript $rscript
```

### Generate study checkpoints

Set `$simSurveySnapshot` to the extracted
`SimSurvey-29a1b03350c274f88dac108fb95f9446eac79233` source. The launcher uses
one below-normal worker with RAM and CPU safeguards. Skip stages with complete
checkpoints.

```powershell
& $runner -Phase Pilot -Grid 40 -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase Study -Grid 40 -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase PopulationPilot -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase PopulationDiagnostic -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase PopulationExtension -SimSurveySource $simSurveySnapshot -Rscript $rscript
```

`Pilot` checks the sampling run and `Study` generates the full sampling
results. The population diagnostic tests observation means on three histories
and supplies some of the independent SD draws. The extension uses
`results/observation_error_seed21/study-grid40/summary.rds` for the Cadigan
calibration.

The earlier single-population assessment experiments are recorded on the
`observation-error-study` branch through `bee37cb3`.

### Mortality and sampling covariance

The [Z covariance study](observation_error_z_report.Rmd) crosses the three
spatial designs with independent, age-correlated, and shared mortality
shocks. It includes a temporal-only contrast under default spatial fields.
Repeat monitoring in years 5 and 16 of four paired histories separates
sampling covariance from covariance between stock states. A small tinyAM
comparison checks whether age correlations in residuals track the known
sampling errors.

```powershell
& $runner -Phase ZPilot -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase ZStudy -SimSurveySource $simSurveySnapshot -Rscript $rscript
& $runner -Phase ZReport -Rscript $rscript
```

The 40 × 40 experiment uses standard effort. The smoke run, full checkpoints,
and [generated report](results/observation_error_z_covariance/study-grid40/observation_error_z_report.html)
are under `results/observation_error_z_covariance`.

The [RV covariance benchmark](rv_covariance_benchmark.md) compares independent,
positive age-AR1 and exchangeable log errors on held-out draws from this study.

### Study to-do

- [ ] Test whether an assessment that estimates both observation and process
  correlation can recover their strengths. Use simulations with known
  observation covariance and repeat surveys to assess when the model can
  separate the two sources.

## Outputs and reproducibility

All results, fits, figures, HTML, warnings and failures stay under ignored
`results/`. Preserve every attempt and seed; do not drop failures, relax tolerances
or retry fits to improve recovery rates. Checkpoints retain provenance and,
where implemented, reject changed source fingerprints. Start a fresh checkpoint
after substantive source/design changes rather than mixing revisions. Report
interval coverage separately from numerical success and unavailable uncertainty.
Generated reports read cached experiments; rendering does not refit models.

The current fleet recovery runner uses `results/multi_fleet/recovery_supported/`.
Its seven isolated and seven full configurations retain their original seeds.
The earlier mixed benchmark and retired exploration caches remain historical;
their provenance is recorded in `multi_fleet_validation.md` and Git history.

Population mean/SD and recruitment covariates now belong in `pop_covariates`.
Recruitment uses `rec_settings`; parent SSB follows supplied maturity timing,
with zero timing in the original start-of-year experiments. These are current
interfaces, not additional changes to the generating models.

## Findings and limits

- **Catchability:** the original study used 2,400 isolated and 480 full attempts.
  Full numerical success ranged from 73–100%. Sparse numeric-by effects could
  converge with underestimated SD, and AR1 persistence was imprecise. Inspect
  support and uncertainty rather than equating convergence with identifiability.
- **F/M means:** 600 isolated fits passed; success among 360 full fits ranged
  from 67–97%. Small shared M effects were hard to separate from residual
  variation. One supported temporal mean is useful; overlapping processes need
  separate evidence. Prefer mean-only M before adding another variance component.
- **Observation SD:** 1,580 attempts included process recovery, noisy tails and
  full assessments. Full passes were 30/30 catch IID, 13/30 catch RW, 30/30 catch
  AR1, 30/30 survey IID, 28/30 survey RW and 29/30 survey AR1. Replication helps;
  sparse age curves and correlation can remain weakly estimated. A correctly
  specified quadratic curve outperformed random curves for the quadratic truth.
- **Recruitment:** all 1,400 isolated fits passed. Full passes were 30/30 BH–IID,
  29/30 BH–AR1, 28/30 Ricker–IID, 27/30 Ricker–AR1, 24/30 covariate–IID and
  28/30 covariate–RW. A corrected study observation-row alignment is protected
  by a regression test; superseded invalid experiments remain historical.
- **Recruitment cross-tests:** 2,498/2,500 isolated and 698/750 full fits passed.
  Ricker–IID self-tests passed 25/30, with beta coverage only 76% among passing
  full fits. RW truth has no curve-parameter recovery target. Weak curve support
  is a scientific limitation even when a curve appears plausible.
- **BH contrast:** lower recruitment noise alone did not reliably identify a
  curve with narrow parent SSB. Wider SSB contrast improved BH recovery. The
  noise study uses SDs 0.15, 0.35 and 0.60; its declared synthetic F covariate is
  held identical between generating and fitted models.
- **Rising F means:** isolated/full passes were mono–IID 96/100 and 22/30,
  mono–AR1 94/100 and 28/30, logistic–IID 100/100 and 30/30, and logistic–AR1
  100/100 and 25/30. Full mean-F RMSEs were 0.0099, 0.0194, 0.0044 and 0.0104;
  mean-curve coverage was 94.6%, 92.5%, 95.0% and 92.9%. All 15 primary full
  failures missed the fixed gradient criterion despite optimizer/Hessian success.
  Exact plateaus and youngest-age logistic midpoints were also tested, for 680
  total attempts. Static mean curves cancel from RW increments and are rejected.
- **Correlated F / latent SDs:** all 600 isolated and 90 informative full fits
  passed. Full median older/young SDs were F 0.199/0.101, N 0.203/0.097 and
  M 0.199/0.085 against 0.20/0.10. F correlation coverage was 26/30. Some young-M
  SDs approached zero with wide intervals; numerical success did not imply
  precise mortality-variance inference. Process-SD coefficients and derived SDs
  are distinct scales, and M groups must respect fitted state blocks.

Recovery reports retain equations, figures, support warnings and interval
calculations. The package help documents the current interfaces and safeguards.
No validation finding automatically promotes an assessment candidate or changes
its retained scientific assumptions.

## Cohort covariates

See [cohort_report.md](cohort_report.md) and `cohort_recovery.R` for the focused
cohort-mean pilot, including mortality/catchability overlap. Generated attempts
remain under ignored `results/`.

## Formula-driven F

`deterministic_f_recovery.R` and [deterministic_f_report.md](deterministic_f_report.md)
cover F means without a separate residual process, including projection and
uncertainty checks. Settled assessment recipes are unchanged.
