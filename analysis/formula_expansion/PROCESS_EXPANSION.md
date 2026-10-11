# Mortality processes, SD sharing and spawning-time SSB

Study branch: `process-expansion`, from `comp-assessments` at `a048fb9`.
Stock recipes and the canonical database are unchanged. Generated fits,
figures, detailed tables and dashboards remain in ignored
`results/process_expansion/`.

## Implemented options

- F `cor_rw`: annual log-F residual increments have AR1 age correlation,
  signed correlation in (-1, 1), marginal age-specific SDs, and no penalty on
  the first row. Zero correlation gives independent RW increments.
- N/F/M `sd_form`: ordinary age-based log-SD formulas evaluated from weight
  covariates. N uses destination ages; M uses unique state blocks. Temporal
  SDs and random-effect SD terms are rejected. Default `~ 1` retains scalar
  parameters and existing likelihood/simulation paths.
- `ssb_settings$spawn_time`: a stock-wide fraction in [0, 1]. Survival through
  that fraction of F + M affects reported SSB and stock–recruit parent SSB in
  the same chronological population path. Default zero preserves previous SSB.

Coefficients and derived SDs have distinct reporting scales. SD profiles and
SSB retain RTMB uncertainty; correlation and weak-SD advisories are separate
from numerical convergence. The dashboard shows SD age profiles and SSB timing.
Accepted N and biology alone no longer generate a purported spawning-time
reference SSB without accepted mortality.

## Recovery evidence

| Experiment | Attempts | Numerical passes |
|---|---:|---:|
| Isolated F correlated increments, common/grouped SD | 200 | 200 |
| Isolated N IID residuals, common/grouped SD | 200 | 200 |
| Isolated M IID residuals, common/grouped SD | 200 | 200 |
| Full assessment, correlated F and grouped SD | 30 | 30 |
| Full assessment, grouped N SD | 30 | 30 |
| Full assessment, grouped M SD | 30 | 30 |

All Hessians passed. Full assessments had 30 years, eight ages, two surveys
and known observation log-SD 0.08: deliberately informative conditions.
Full-assessment median older/young SDs were F 0.199/0.101, N 0.203/0.097 and
M 0.199/0.085, against truth 0.20/0.10. Young-age M SD sometimes approached
zero; four M attempts emitted the new wide-SD-interval advisory. Numerical
success did not guarantee precise variance recovery.

Isolated coefficient interval coverage ranged from 91% to 98%. Full F/N
SD-coefficient coverage was 90–93%; F correlation coverage was 26/30.
M coverage was 93–100%, with some very wide intervals. Thirty replicates
do not establish general interval calibration or performance with sparse data.

## Source-supported trials

All optimizer codes were zero and all Hessians passed. Maximum gradients:

| Assessment | Baseline | Individual changes | Combined |
|---|---:|---|---:|
| EBS pollock | 0.00194 | Spawning fraction 0.25: **0.03695 (fails 0.01 criterion)** | Not applicable |
| Northern Shelf haddock | 0.00098 | Correlated F: 0.00109; separate plus-group N SD: 0.00137 | 0.00080 |
| North Sea plaice | 0.00045 | Correlated F: 0.00109; four source F-SD groups: 0.00038 | 0.00048 |

EBS's pinned source spawning month supports 0.25. Its structured selectivity
does not justify claiming an exact AR1 age-correlated F-RW mapping. Spawning
time leaves this recipe's likelihood bit-for-bit unchanged at the same
parameters because recruitment has no SSB curve. Its warm-start restart
missed the gradient criterion; this attempt was retained without rescue.
Approximate SSB mean absolute percent difference fell from 39.6% to 27.3%;
trajectory correlation remained about 0.896. Biological definitions remain
labelled approximate. Population states and residuals were essentially unchanged.

Haddock's baseline F is stationary AR1; the correlated-RW trial replaces that
process rather than only adding correlation. Alone it did not improve the
overall comparison. Source-supported N sharing improved SSB scale difference
from 24.2% to 3.0% and trajectory correlation from 0.966 to 0.999; N and F
scale agreement also improved. Recruitment scale difference worsened from
17.2% to 23.3%, and the largest survey residual rose from 3.57 to 6.73 SDs.
The combined trial did not resolve that trade-off.

Plaice's correlated F improved F scale difference from 22.9% to 20.3%, but
worsened SSB (11.2% to 15.6%), N and recruitment scale agreement. Source F-SD
groups improved recruitment slightly, with little SSB change and worse F
agreement. Combinations also traded metrics. Correlation estimates had finite
intervals, but existing outliers, residual dependence and strongly correlated
fixed parameters persisted. No candidate is automatically promoted.

## Reproduce and inspect

From the repository root, run `validate_process_expansion.R` (defaults: 100
isolated and 30 full replicates), then `trial_process_assessments.R`, both
under this directory. Render `process_expansion_report.Rmd` for the figures,
coverage tables, comparison metrics and residual summaries. Every attempt,
warning and error is retained locally. The three assessment dashboards are
in `results/process_expansion/assessments/`.

Validation passed 1,965 package assertions, with one interactive-only skip;
all 32 comparative-assessment test scripts and general database validation
passed. `R CMD check --no-manual` finished with zero errors, warnings or notes.
A preliminary stricter `--as-cran` check reported the pre-existing undeclared
`withr` test dependency, top-level `AGENTS.md`, and inability to verify current
time. Those unrelated packaging/environment findings were left unchanged.

Retain these capabilities, with uncertainty checks and source-definition
review. F/N recovery is encouraging; M variance inference remains more
fragile. The application evidence supports user review of trade-offs, not
automatic recipe changes or promises of closer assessment agreement.
