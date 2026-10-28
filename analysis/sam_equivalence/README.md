# Northern Shelf haddock: SAM likelihood equivalence

This focused comparison tests whether tinyAM can reproduce the **2026 official
Northern Shelf haddock SAM run**, `NShaddock_WGNSSK2026_Run1`, over 1972–2026,
ages 0–8+. It does not alter comparative-assessment recipes or result tables.

## Reproduce

Use the native model object retained by the source review and stockassessment
0.12.0 at **`1cc464b80f6f2a4a38b35cce276f13e521761da3`**. The script verifies
the installed revision and input MD5 (`3ae172debf40dbc59d202335d9d08371`).
The local pinned installation was built with R 4.4.1; use a compatible runtime.

```r
# From the repository root:
# Rscript analysis/sam_equivalence/compare_haddock.R path/to/run_model.RData path/to/SAM-library
```

The default output directory is `results/haddock/` beneath this document and is
ignored. It retains the fitted object, numerical check tables, warnings and a
self-contained HTML dashboard. A failed run records its error and warnings.
The script uses source parameters as starts for a numerical equivalence test;
it does not claim a test of default-start convergence or tune failed models.

Sources:

- [Official 2026 assessment chapter](https://ndownloader.figshare.com/files/66104291).
- [Official SAM run](https://stockassessment.org/datadisk/stockassessment/userdirs/user259/NShaddock_WGNSSK2026_Run1/).
- [Pinned SAM implementation](https://github.com/fishfollower/SAM/tree/1cc464b80f6f2a4a38b35cce276f13e521761da3/stockassessment).

## Matched specification

Both models integrate all first-year N states without an initial density, use
random-walk recruitment, independent cohort-survival innovations with a shared
age-1–7 SD and a separate 8+ SD, and age-correlated F random-walk increments.
There is one F state per age, one F innovation SD, and no first-year F density.
SAM's correlation transform is `tanh(itrans_rho)`, identical to tinyAM's
`tanh(atanh_rho_f)`. The parameter mapping is otherwise a permutation.

Both use the same supplied M, stock weights, maturity, zero spawning fractions,
observation values, survey times, relative precision, 17 baseline-q coefficients,
three positive q powers, and 12 observation-SD groups. Power affects abundance
at survey time: `log I = log q0 + p * (log N - s * Z) + error`.
No reference estimates constrain the fitted parameters.

SAM has no catch-weight input for the survey-only year 2026. The script repeats
2025 catch weights solely to satisfy tinyAM's complete input grid. This enters
neither the numbers likelihood nor SSB; **2026 yield is not compared**. Native
stock weights, M and maturity for 2026 are retained unchanged.

The initial source preparation also excludes the unused recruitment-SD category
from the older-cohort SD design. Recruitment retains its own source SD; no
active variance parameter or density is changed.

## Numerical result

Checks compare the independently compiled, pinned SAM likelihood with RTMB at
the native fitted parameter vector and a reproducible small perturbation.

| Quantity | Result |
|---|---:|
| Included observations | 1,153 |
| Estimated outer parameters | 37 |
| Integrated N/F states | 990 |
| Maximum absolute joint-objective difference | 7.96e-13 |
| Maximum absolute gradient difference | 8.42e-12 |
| Maximum absolute full Hessian difference | 1.17e-10 |
| Maximum absolute log-observation-prediction difference | 3.56e-15 |
| Laplace objective, both implementations | 650.04377631258 |
| Optimized objective difference | -4.55e-13 |
| Maximum absolute outer-parameter difference | 5.11e-11 |
| Maximum relative N difference | 1.10e-10 |
| Maximum relative F difference | 9.26e-11 |

The full joint Hessian check includes the random-state block and its cross
derivatives. The tinyAM fit passes numerical checks: optimizer code 0, maximum
absolute gradient 9.82e-7 and a positive-definite Hessian.

Scientific advisories remain: F-increment correlation is within 0.05 of a
boundary, three fixed-parameter pairs exceed 0.95 absolute correlation, some
standardized residuals exceed three SDs, and some groups have substantial
adjacent-year residual correlation. Equivalence does not validate these
assumptions or establish that the model is suitable for operational use.

The source's arithmetic Fbar and tinyAM's abundance-weighted Fbar are different
reported summaries. Native SAM interval tables also use different conventions
from tinyAM's default 95% Wald limits, and conditional state SEs must not be
compared with SEs that propagate fixed-parameter uncertainty. Align these
definitions before comparing summaries; they do not change the likelihood.

## Package validation

On R 4.6.1, the complete package suite passed all 3,807 expectations with no
failures, warnings or skips. `R CMD check --no-manual` returned **Status: OK**,
including examples and documentation checks. Focused regression tests cover
latent recruitment boundaries, power-formula validation and derivatives,
simulation, warm starts, projections, hindcasts, dashboard rendering and OSA
row identities. These checks complement the pinned SAM comparison above;
they do not remove its statistical advisories or default-start limitation.
