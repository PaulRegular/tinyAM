# SAM likelihood equivalence checks

These tests compare tinyAM with native [SAM](https://github.com/fishfollower/SAM) using pinned SAMjr fixtures and the official 2026 Northern Shelf haddock assessment. **All supported comparisons pass.** This establishes numerical equivalence for the tested configurations, not general model validity or convergence from default starts.

## SAMjr printed-output tests

```sh
Rscript analysis/sam_equivalence/compare_samjr.R SAM-library testmore
Rscript analysis/sam_equivalence/test_testmore.R
```

Using unchanged SAMjr inputs and expected outputs (revision `0384774a3f2f8ec3dc3330375998bd73ef6bbe7c`), all three implementations reproduce both scripts' printed results exactly:

| Script | Matched output |
|---|---|
| `codIN3Ben` | Objective `751.6052` |
| `nscodXtraSd` | Log-likelihood, parameter count, AIC and LR test (*p* = `0.3031836`) |

Implementations are ordinary `fit_tam()`, tinyAM's likelihood with SAMjr optimizer controls, and native SAM (**six passing comparisons**). Each starts from the source's **unoptimized `defpar()` values**, not a fitted SAM solution. The protocol optimizer uses `nlminb` (`eval.max = iter.max = 2000`, `rel.tol = 1e-10`) and three Newton steps for every fit. All fits pass convergence/Hessian checks. No scientific settings were tuned.

The runner verifies fixture hashes and exact printed output (ignoring trailing carriage returns only), saves fits/checks/warnings under ignored `results/testmore/`, and exits nonzero on failure. This is **printed-output equivalence**, not bitwise identity. Harness tests pass on R 4.4.1 and 4.6.1; the complete comparison was rerun on R 4.4.1.

`samjr_testmore.R` adapts attributed GPL-2 source workflows; [SOURCES.md](SOURCES.md) documents changes. The likelihood is implemented by tinyAM; native SAM provides input/configuration utilities. SAMjr readers are not tested. The full `nscodbiopro` output script is excluded because its later biology processes and stochastic forecasts are unsupported here.

## Native SAM likelihood comparisons

Run from the repository root with R 4.4 and the pinned SAM library:

```sh
Rscript analysis/sam_equivalence/compare_samjr.R SAM-library all
```

Use `baseline` instead of `all` for only `nscodbiopro`'s supplied-biology `fit0`; `all` includes the printed-output tests. Sources are pinned to SAMjr `0384774a3f2f8ec3dc3330375998bd73ef6bbe7c` and SAM `1cc464b80f6f2a4a38b35cce276f13e521761da3` (stockassessment 0.12.0). `fixtures.csv` records provenance and hashes. Results and failures are retained under ignored `results/`.

**All four specifications pass 33 checks each (132 total)**, including objectives, derivatives, Laplace likelihoods, parameter estimates/SEs, N/F, predictions, AIC and aligned stock summaries. Both packages pass a 0.01-SD perturbed-start fit and two retrospective peels per specification.

| Specification | Max. joint objective difference | Optimized difference | Max. scaled SE difference |
|---|---:|---:|---:|
| `nscodbiopro`: supplied biology | 2.28e-13 | 6.83e-13 | 4.50e-12 |
| `nscodXtraSd`: standard | 2.28e-13 | 6.83e-13 | 1.31e-12 |
| `nscodXtraSd`: extra SD | 5.69e-13 | 1.60e-12 | 2.29e-12 |
| `codIN3Ben`: unfished young ages | 6.83e-13 | 1.71e-12 | 1.15e-9 |

Tolerance limits are `1e-7` for joint objectives, `1e-6` for scaled derivatives, `1e-5` for Laplace objectives and `1e-4` for optimized objectives/scaled summaries. Scaled differences divide by `1 + abs(native value)`; relative N differences use native N. State and observation counts must match exactly. The largest observed scaled derivative and Laplace differences are `1.60e-12` and `1.71e-12`, respectively; maximum relative N difference across perturbed/peel fits is `3.26e-5`. Comparisons initialize from **native fitted parameters**, so do not test tinyAM's default starts.

Model-specific handling is unchanged: `nscodXtraSd` binds F at ages 5–6 and estimates a multiplicative 1993–2005 catch SD using `catch_settings$sd_form`; `codIN3Ben` fixes F to zero at ages 1–2 and shares F at ages 9–10. The latter has 115 missing observations represented as auxiliary Gaussian states in SAM. Their independent densities integrate to one and are subtracted, with derivatives, only for common-state joint comparisons. The original native fit and full Laplace objective remain intact; random-state counts differ by 115. Survey-only years may require padded catch weights; their yield is not compared. Fbar is aligned using tinyAM's abundance-weighted definition, not SAM's arithmetic Fbar intervals.

Render a saved comparison without refitting:

```sh
Rscript analysis/sam_equivalence/render_samjr.R analysis/sam_equivalence/results/samjr/nscodXtraSd/extra_sd.rds
```

OSA residual and process-draw equivalence remain untested. An exploratory, non-fixture configuration estimated catch SD near `4e-5` and failed positive-definite uncertainty checks despite matching objectives. **Likelihood equivalence does not establish identifiability**; AIC/cAIC normalization cautions still apply.

## Northern Shelf haddock (2026)

The official `NShaddock_WGNSSK2026_Run1` SAM model (1972–2026, ages 0–8+) is reproduced using the source's fitted parameters and a small perturbation. Run with a compatible R 4.4 installation and pinned SAM `1cc464b80f6f2a4a38b35cce276f13e521761da3`:

```sh
Rscript analysis/sam_equivalence/compare_haddock.R path/to/run_model.RData path/to/SAM-library
Rscript analysis/sam_equivalence/render_haddock.R
```

The script checks the source input MD5 (`3ae172debf40dbc59d202335d9d08371`). Results, warnings and an HTML dashboard are saved under ignored `results/haddock/`.

Both implementations integrate first-year N without an initial density, use random-walk recruitment, independent survival innovations (one SD for ages 1–7, another for 8+), and age-correlated F random walks with one state per age. They retain the same M, maturity, weights, survey timing, zero spawning fractions and relative precision, including 17 baseline catchability coefficients, three positive survey powers and 12 observation-SD groups. The F-correlation transformations are equivalent; no source estimates are fixed as parameter constraints. The source's 2026 survey-only year lacks catch weights: 2025 weights are repeated for tinyAM's input grid, but **2026 yield is excluded**.

| Measure | Difference or value |
|---|---:|
| Observations / fixed parameters / N-F states | 1,153 / 37 / 990 |
| Maximum joint-objective difference | 7.96e-13 |
| Maximum gradient / full-Hessian difference | 8.42e-12 / 1.17e-10 |
| Laplace objective (both models) | 650.04377631258 |
| Optimized-objective difference | -4.55e-13 |
| Maximum parameter / relative N / relative F difference | 5.11e-11 / 1.10e-10 / 9.26e-11 |

The tinyAM fit converges (optimizer code 0, maximum gradient `9.82e-7`) with a positive-definite Hessian. The joint Hessian comparison includes the random-state block and cross derivatives. Scientific cautions remain: F-correlation approaches a boundary, three parameter pairs are highly correlated, and some residuals are large or serially correlated. SAM's arithmetic Fbar differs from tinyAM's abundance-weighted Fbar. SAM intervals also differ from default 95% Wald intervals; conditional state SEs should not be compared with SEs that propagate fixed-parameter uncertainty. **Numerical agreement does not validate model assumptions.**

Sources: [assessment chapter](https://ndownloader.figshare.com/files/66104291), [official SAM run](https://stockassessment.org/datadisk/stockassessment/userdirs/user259/NShaddock_WGNSSK2026_Run1/), [pinned SAM code](https://github.com/fishfollower/SAM/tree/1cc464b80f6f2a4a38b35cce276f13e521761da3/stockassessment).

