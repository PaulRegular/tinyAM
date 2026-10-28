# SAMjr-derived test code

Upstream: [SAMjr](https://github.com/fishfollower/samjr), revision
`0384774a3f2f8ec3dc3330375998bd73ef6bbe7c`. Its package DESCRIPTION names
Anders Nielsen and Christoffer Moesgaard Albertsen as authors and specifies
GPL-2. The adapted `samjr_testmore.R` file retains that GPL-2 license; the
package's license and model implementation are unchanged. The entire analysis
directory is excluded from the tinyAM package build. The GPL-2 terms are
included in [COPYING.SAMJR](COPYING.SAMJR).

Sources, all at the pinned revision:

| Source | Use and changes |
|---|---|
| [`testmore/codIN3Ben/script.R`](https://github.com/fishfollower/samjr/blob/0384774a3f2f8ec3dc3330375998bd73ef6bbe7c/testmore/codIN3Ben/script.R) | Original ICES input preparation, configuration and `defpar()` values are evaluated. The fitting call is intercepted, then replaced with the corresponding tinyAM specification. The original objective-only output expression is retained. |
| [`testmore/nscodXtraSd/script.R`](https://github.com/fishfollower/samjr/blob/0384774a3f2f8ec3dc3330375998bd73ef6bbe7c/testmore/nscodXtraSd/script.R) | Same substitution for both Standard and XtraSd models; the original comparison-table output is retained. |
| [`samjr/R/tables.R`](https://github.com/fishfollower/samjr/blob/0384774a3f2f8ec3dc3330375998bd73ef6bbe7c/samjr/R/tables.R), `modeltable.default()` | Adapted table assembly, sorting, likelihood-ratio calculation and column labels. Extraction uses the fitted objective and active outer-parameter count, shared by tinyAM and native SAM. |
| [`samjr/R/fit.R`](https://github.com/fishfollower/samjr/blob/0384774a3f2f8ec3dc3330375998bd73ef6bbe7c/samjr/R/fit.R), `sam.fit()` | A separately labelled protocol run adapts its `nlminb()` settings and three subsequent Newton steps, applied to a fresh tinyAM RTMB objective. Ordinary `fit_tam()` results are also tested, without these additional steps. |
| [`Makefile`](https://github.com/fishfollower/samjr/blob/0384774a3f2f8ec3dc3330375998bd73ef6bbe7c/Makefile), `testmore-%` | Exact generated-output comparison against `res.EXP`, ignoring only trailing carriage returns as in `diff --strip-trailing-cr`. Implemented in R for Windows portability. |

The source `.dat`, `.cfg`, `script.R` and `res.EXP` files are downloaded into
ignored local caches. Their blob identifiers and canonical content hashes are
tracked in `fixtures.csv`; expected outputs are never regenerated from tinyAM.
The existing source fixtures and the new expected files are checked before use.

The ICES readers, configuration loader and `defpar()` come from pinned native
SAM (`stockassessment` 0.12.0, revision
`1cc464b80f6f2a4a38b35cce276f13e521761da3`). This deliberately tests tinyAM's
models, not its ability to reproduce SAMjr's input-reader implementations.
No SAMjr likelihood code is copied or used to calculate tinyAM predictions.
`samjr_helpers.R` is tinyAM-specific translation and numerical comparison code.

Only the two complete scripts above count as upstream output tests.
`nscodbiopro` is retained as a separate supplied-biology likelihood comparison;
its full script's biology processes and stochastic forecasts are not supported
and are not reported as a passing upstream test.

An exact output match establishes agreement at the upstream printing precision,
not bitwise equality of every fitted double. Differences, warnings and failed
attempts remain in the ignored output directory, with the two optimizer paths
reported separately. No scientific configuration or reference output is tuned.
