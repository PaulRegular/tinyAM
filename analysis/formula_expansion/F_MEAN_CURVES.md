# Recruitment settings and rising F means

Recruitment now has its own `rec_settings` list. Process and mean are separate;
IID/AR1 require an explicit mean, including `~ 1` for a constant mean. The
default RW, Gaussian densities, stock–recruit equations, lags and boundary
states are unchanged. Six pre-refactor cases retained parameter layouts,
designs, likelihoods, gradients and reports within `1e-10`.

F mean formulas now permit one `mono()` or `logistic()` curve with IID or AR1
residuals. The former uses non-negative log-F increments; the latter multiplies
baseline F by logistic selectivity. These constrain the typical age pattern,
rather than every annual F realization. Static curves cancel from RW increments
and are rejected there. M curves and overlapping temporal means with AR1
residuals remain outside the enabled scope. No stock recipe was retuned.

## Recovery evidence

Run `f_mean_curve_validation.R` from the repository root. It retains every
attempt, seed, warning and failure under ignored `results/f_mean_curves/`.
Render `f_mean_curve_report.Rmd` for tables and figures. Code and report text
are tracked; fits, detailed results and HTML are local.

The study used 100 isolated Gaussian and 30 full-population replicates per
primary combination. Full models had 30 years, eight ages, two surveys, known
M and observation log-SDs of 0.08. Numerical success requires optimizer code
zero, a positive Hessian and adjusted maximum gradient below `0.01`.

| Curve / residual process | Isolated success | Full success | Full mean-F RMSE | Mean-F interval coverage |
|---|---:|---:|---:|---:|
| Mono / IID | 96/100 | 22/30 | 0.0099 | 94.6% |
| Mono / AR1 | 94/100 | 28/30 | 0.0194 | 92.5% |
| Logistic / IID | 100/100 | 30/30 | 0.0044 | 95.0% |
| Logistic / AR1 | 100/100 | 25/30 | 0.0104 | 92.9% |

Recovery summaries include attempts failing numerical checks when estimates
and uncertainty were available. Pooled curve-coefficient coverage was
93.3–98.9% in full models; local boundary intervals are not boundary-adjusted
inference.

Mono/IID and logistic/AR1 missed the 90% full-fit target. Every primary full
fit had optimizer code zero and a positive Hessian. All 15 failures were
gradient misses, between 0.0107 and 0.0184. Repeating two failed seeds solely
to audit gradients gave AD/central-difference discrepancies below `4e-7`.
This supports numerical stopping precision as the issue; it does not make
the failed checks pass. No tolerance, optimizer or model was adjusted.

Twenty replicates per stress case were run at both levels, for 680 total
attempts. Exact mono plateaus passed 20/20 and 19/20 isolated checks, and
18/20 and 19/20 full checks for IID and AR1. Logistic curves with their
midpoint at the youngest observed age passed all 80 stress checks. This
particular stress retained enough contrast; it does not establish reliable
recovery when all observed ages lie on a nearly flat asymptote.

## Recommendation

Retain both curves as source-supported options with the existing diagnostics.
The tested equations and simulations are consistent, and curve recovery was
useful in these informative settings. Numerical success was imperfect, so
inspect `check_tam()` and uncertainty rather than assuming every specification
will fit reliably. Do not use either curve to force agreement with an accepted
assessment or infer an unobserved selectivity limb without supporting data.

## Validation

All 33 comparative-assessment test scripts passed, including interactive
previews and parallel execution. The final installed-package suite passed
2,230 assertions with one interactive-only skip. `R CMD check --no-manual`
finished with status OK. Recruitment-study smoke tests and the six numerical
parity checks also passed. Mean-F dashboards retain their own uncertainty
ribbons and omit unavailable reference means without changing annual-F plots.
