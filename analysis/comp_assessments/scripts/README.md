# Assessment scripts

Run these scripts from the repository root in RStudio.

## Database

Scripts in `database/` seed, import, validate, and report the assessment database.
For GOA pollock, run `032_import_goa_pollock_catch_numbers.R` after the base
import to extract the detailed age-1–15 catch table from the cached final SAFE
and refresh its canonical input rows.

## Translation and fitting

Open and source `translation/run_stock.R` to fit one assessment interactively.
Change `assessment_id` in that file first. It reads the working-tree database
and leaves `obs`, `settings`, `fit`, `ref`, `audit`, `diagnostics`, and
`comparison` in the workspace.

Source `translation/run_all.R` for the reproducible all-assessment run. It reads
the committed database snapshot, fits with four multisession workers, and writes
the aggregate diagnostics and comparison summaries from the main R process.

Functions can also be called directly:

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")

x <- run_assessment("dfo_cod_2j3kl_2025")
runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE
)
```

Use `run_assessment(..., fit = FALSE)` to inspect the translated observations
and settings before fitting. Use `run_assessments(parallel = FALSE)` when
debugging a multi-assessment run.

Run `translation/check_northern_cod_indices.R` to compare the RV-only,
RV/Smith Sound and juvenile-index trial fits. It writes compact diagnostics
to `results/northern_cod_indices.csv` without saving fitted objects.
