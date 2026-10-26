# Assessment scripts

Run these scripts from the repository root in RStudio.

## Database

`database/` contains reproducible source imports, extraction helpers and
validators. Retain scripts that recover canonical values or assumptions even
when their names contain `review`. They are curation code, not model trials.
Run `002_validate_database.R` to validate the canonical database;
`003_observation_readiness.R` regenerates the observation-readiness table.
Source-dependent validators require the corresponding ignored source cache.

## Translation and fitting

Open `translation/run_stock.R`, choose `assessment_id`, and source it. Each
stock script prepares observations and shows its literal `fit_tam()` call.
The workspace retains `obs`, years, ages, `fit`, background and `assessment`.
Rerun the Preview line after editing a fit; previewing never refits it.

`run_assessment(..., fit = FALSE)` prepares observations without final or
preliminary fits. Model settings are stored in a fitted object's `dat`.

`translation/run_all.R` reads the committed database, fits with four workers,
and writes aggregate results only in the parent process. It caches returned
fits before optional dashboards, then runs the two overview-figure scripts.
The controls at the top select parallel execution, worker count, dashboards
(off by default), and a 20-minute per-assessment deadline. A finite deadline
requires `callr` and runs each assessment in a supervised R subprocess, including
when `parallel = FALSE`. A timed-out process is stopped and the batch continues
with an explicit `timed_out` diagnostic; settings and starting-value rules are
unchanged. Set `timeout <- Inf` for unbounded fitting.

Stock starts and completions are printed as they happen. With `cache = TRUE`,
each completed stock is saved immediately as `results/cache/<id>/result.rds`,
alongside its fit (if available), `run.log` with stage timestamps, and
`run.log.output` with subprocess output. These checkpoints do not wait for the
whole batch. Existing old fit files are retained; `result.rds` records the current
run's outcome and database revision. `cache` saves results; it does not skip fits.
Checkpoint replacement happens only after serialization completes, so an
interrupted save does not overwrite a previous complete checkpoint.
Reduce worker count when memory is limited; use serial execution for debugging.

Current summary CSVs and overview figures are tracked. Fits, dashboards,
source files, detailed trial outputs and rendering failures remain local and
ignored. The complete batch alone refreshes aggregate results.
