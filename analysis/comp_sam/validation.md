# Validation record

Completed 2026-09-30 on Windows 11, R 4.4.1.

| Check | Result |
|---|---|
| Offline SAM tests (`testthat::test_local(filter = "sam")`) | 93 assertions passed; no failures, warnings or skips |
| Full suite in final built-package check | 1,154 assertions passed; no failures or warnings; one interactive-only dashboard check skipped |
| `R CMD build` | Passed |
| `R CMD check --no-manual` | **Status: OK**; no errors, warnings or notes |
| Pinned public input checksums | All 14 files verified |
| `check_obs()` on translated North Sea cod inputs | Passed |
| Simplified tinyAM model | Finite initial joint objective and gradient; no fit claimed |
| Installed SAM input/configuration comparison | Passed |
| Saved/current reference tables versus SAM table functions | Passed |
| Fresh SAM fit | Optimizer code 0; finite SDs; positive-definite Hessian; maximum absolute gradient 4.07 × 10⁻¹¹ |
| Final diff whitespace check | Passed |

The final package check includes the fixed-q formula regressions, rejection of
unfitted SAM reference objects and the initial-state integration audit. The
interactive skip checks browser launching, while the noninteractive HTML-render
tests passed. PDF manual generation was excluded. Pandoc printed deprecation
notices about its highlighting option; these did not produce test/check warnings.

SAM integration used installed stockassessment 0.12.0 and TMB 1.9.16, with SAM
source commit `c6cfd035c7de59f7b3421dde31901efbca4cb0e8`. The workflow is in
`001_download.R`, `002_nscod.R` and optional `003_validate_sam.R`. Exact fit
diagnostics are in `results/SAM_current_diagnostics.csv`.

Work is on `comp-sam`, based on `main` commit
`a06e309f999cfc62d08e42cd9358f974996ad7f7`. No existing core model R files,
defaults or dependencies were changed. No merge was performed.

The audit is specific to the selected public repository example. Scientific
choices remain about F innovation correlation, F variance groups, catch scaling,
initial-state integration and comparable biomass/Fbar definitions; see README.md.
