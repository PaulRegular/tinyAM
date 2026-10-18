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
fits before rendering dashboards, then runs the two overview-figure scripts.
Reduce worker count when memory is limited; use serial execution for debugging.

Current summary CSVs and overview figures are tracked. Fits, dashboards,
source files, detailed trial outputs and rendering failures remain local and
ignored. The complete batch alone refreshes aggregate results.
