# Assessment scripts

Run these scripts from the repository root in RStudio.

## Database

Scripts in `database/` seed, import, validate, and report the assessment database.

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
runs <- run_assessments(parallel = TRUE, workers = 4)
```

Use `run_assessment(..., fit = FALSE)` to inspect the translated observations
and settings before fitting. Use `run_assessments(parallel = FALSE)` when
debugging a multi-assessment run.
