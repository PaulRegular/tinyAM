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
| Recruitment recovery | `simulate_recruitment.R`, `recruitment_validation.R` | `recruitment_report.Rmd` |
| Recruitment self-/cross-tests | `simulate_recruitment_cross.R`, `recruitment_cross_validation.R` | `recruitment_cross_report.Rmd` |
| BH noise and parent-SSB contrast | `simulate_recruitment_noise.R` | `recruitment_cross_report.Rmd` |
| Rising F means | `f_mean_curve_validation.R` | `f_mean_curve_report.Rmd` |
| Correlated F and process-SD sharing | `validate_process_expansion.R` | `process_expansion_report.Rmd` |

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

## Outputs and reproducibility

All results, fits, figures, HTML, warnings and failures stay under ignored
`results/`. Preserve every attempt and seed; do not drop failures, relax tolerances
or retry fits to improve recovery rates. Checkpoints retain provenance and,
where implemented, reject changed source fingerprints. Start a fresh checkpoint
after substantive source/design changes rather than mixing revisions. Report
interval coverage separately from numerical success and unavailable uncertainty.
Generated reports read cached experiments; rendering does not refit models.

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
