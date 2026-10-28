# SAM likelihood equivalence checks

## SAMjr fixture comparisons

Run from the repository root with R 4.4 and the pinned SAM library:

```sh
Rscript analysis/sam_equivalence/compare_samjr.R SAM-library all
```

Use `baseline` instead of `all` for only `nscodbiopro`'s supplied-biology `fit0`.
Sources are pinned to SAMjr
`0384774a3f2f8ec3dc3330375998bd73ef6bbe7c` and native SAM
`1cc464b80f6f2a4a38b35cce276f13e521761da3` (stockassessment 0.12.0).
`fixtures.csv` records Git blob identifiers and content hashes; the runner
downloads missing files and checks hashes, normalizing line endings and trailing
empty lines for portability. Fits, warnings, input provenance and check tables
remain under the ignored `results/` directory. Every attempted fit is retained;
failures and failed tolerances give a nonzero exit status.

All four fitted specifications pass the 33 checks each (132 in total). No scientific settings
or failed configurations were tuned. Native fitted values initialize the main
comparison; one reproducible perturbation of 0.01 SD on the unconstrained
coordinates and two calendar-year retrospective peels are also fitted by both
packages. These are not tests of tinyAM's default starting values.

| Specification | Maximum joint-objective difference | Optimized-objective difference | Largest scaled SE difference | Largest relative N difference across perturbed/peel fits |
|---|---:|---:|---:|---:|
| `nscodbiopro`: supplied-biology baseline | 2.28e-13 | 6.83e-13 | 4.50e-12 | 8.33e-6 |
| `nscodXtraSd`: standard | 2.28e-13 | 6.83e-13 | 1.31e-12 | 2.50e-6 |
| `nscodXtraSd`: extra SD | 5.69e-13 | 1.60e-12 | 2.29e-12 | 2.31e-6 |
| `codIN3Ben`: unfished young ages | 6.83e-13 | 1.71e-12 | 1.15e-9 | 3.26e-5 |

Checks cover joint objectives, gradients and Hessians at native and perturbed
coordinates; Laplace objectives; estimated fixed/common random-coordinate
counts; observed-row identities/counts; fitted parameters and their SEs; N/F
surfaces; observation predictions; marginal AIC; and aligned recruitment, SSB,
biomass, catch, catchability and abundance-weighted Fbar. Normal 0.975 quantiles align interval
probabilities. Fbar is calculated using native N and F with tinyAM's weighting;
SAM's native arithmetic Fbar intervals are not treated as comparable intervals.

The preset limits are absolute joint-objective difference 1e-7, scaled derivative
differences 1e-6, Laplace-objective difference 1e-5, and optimized objectives and
appropriately scaled summaries 1e-4. Scaled differences divide by
`1 + abs(native value)`; relative N differences divide by native N. Count
differences must be exactly zero. These are numerical-tolerance checks, not
claims of identical floating-point bits.

Across the four fixtures, the largest scaled derivative difference was
1.60e-12 and the largest Laplace-objective difference was 1.71e-12.
The largest scaled difference among the aligned estimate/interval tables was
6.51e-10 (recruitment intervals); catchability-table differences were at most
1.34e-13. All fitted and perturbed-start attempts and both peels per fixture
passed the numerical-success checks.

`nscodXtraSd` binds actual F at ages 5--6. Its additional observation SD is an
estimated multiplicative factor for catch observations in 1993--2005, exactly
represented by a numeric indicator in `catch_settings$sd_form`. `codIN3Ben`
fixes F at ages 1--2 to zero and shares F at ages 9--10. It also contains 115
missing observations integrated as auxiliary Normal variables in native SAM.
Their independent densities integrate to one. Joint comparisons subtract those
auxiliary densities and their derivatives, comparing only the common states;
the original native fit and full Laplace objective are unchanged. Native and
tinyAM random counts therefore differ by those 115 auxiliary coordinates.

Survey-only final years without native catch weights are padded solely for
tinyAM's complete input grid; those years' yield is not a source comparison.
Native retrospective calls use common calendar cutoffs across fleets, rather
than removing each fleet's own final observations. Source scripts are evaluated
only through the selected fits; later biological processes are not enabled.

Render a saved two-model dashboard without refitting:

```sh
Rscript analysis/sam_equivalence/render_samjr.R analysis/sam_equivalence/results/samjr/nscodXtraSd/extra_sd.rds
```

It writes a standalone comparison HTML and native `tam_ref` beside the RDS.
Neither includes newly calculated OSA residuals. OSA and seeded process-draw
equivalence are deferred; matching Gaussian samples requires matching parameter
ordering and factorization, not only matching seeds.

An exploratory package-data trial combining unfished youngest ages, shared
older F, latent initial abundance and a freely estimated catch SD drove catch
SD to about 4e-5. Projection updates retained identical objectives/parameters
but failed positive-definite uncertainty checks. That trial is not an official
SAM fixture and was not rescued by tuning; its warnings remain with the ignored
validation output. The routine integration test uses the established `exp`
initialization. Likelihood equivalence does not establish identifiability, and
intrinsic-process AIC/cAIC normalization warnings remain applicable.

## Northern Shelf haddock

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

To render both saved fits together without refitting, run
`Rscript analysis/sam_equivalence/render_haddock.R`. This creates
`results/haddock/haddock-comparison.html` and `SAM-reference.rds`. The SAM
`tam_ref` uses native states, predictions, coefficients and available uncertainty;
its Background page identifies aligned summary definitions and unavailable
intervals. The curves overlap to numerical precision.

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

On R 4.6.1, the complete package suite passed all 3,867 expectations with no
failures, warnings or skips. `R CMD check --no-manual` returned **Status: OK**,
including examples and documentation checks. The focused F-block file passed
60 expectations. Regression coverage includes singleton/default density and
gradient parity, explicit versus legacy M boundaries, structural zeros, unique
state densities, finite derivatives, single-block correlation removal, warm
starts, simulations, catch targets, information criteria and dashboard rendering.
The existing suite also exercises catch weights, plus-age biology, projections,
retrospectives, hindcasts and OSA extraction.

The 69 standalone comparative-assessment scripts were run without changing
recipes: 59 passed and 10 failed. The failures were reproduced on unchanged
`main` (`96c670d`):

- `test_atka_translation.R`, `test_bc_canary_translation.R`,
  `test_ebs_cod_translation.R`, `test_ebs_translation.R`,
  `test_goa_cod_translation.R`, `test_goa_translation.R`,
  `test_pacific_M_priors.R` and `test_recipe_blocks.R` have expectations that
  differ from the current recipes or observation structure.
- `test_coastal_hake_translation.R` fails with the local RTMB/R error
  `attempt to set an attribute on NULL` while setting `environment(partial)`.
- `test_age_subsets.R` fails its direct reduced-age fit's positive-definite
  uncertainty requirement on both branches. Earlier reduced-age checks pass.
  The original default F calculation was retained after an equivalent
  simplification changed rounding enough to fail an earlier fit in this file.

The dashboard preview script passed on direct rerun after the isolated runner
hit its 180-second limit during concurrent validation. All automated renders use
`open_file = FALSE`. The saved SAM comparison HTML contains both model traces
and 56 widgets; static widget and reference-table checks passed. Interactive
browser inspection could not be completed because the available browser tool
rejects local-file URLs. Rotation and control interaction remain a manual check.

These results complement the pinned SAM comparison; they do not remove its
statistical advisories, exploratory uncertainty failure or default-start limitation.
