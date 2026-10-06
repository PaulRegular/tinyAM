# Comparing assessment fits

This analysis translates selected production age-structured assessments into tinyAM models and compares their reported population estimates. The database preserves source values and assumptions; stock scripts make the tinyAM-specific choices. A close fit is an empirical comparison, not a claim that the models use identical likelihoods or assumptions.

## Scope and provenance

The assessment database was constructed with substantial agentic-AI assistance.
Source records have not all been exhaustively human-verified; check values,
definitions and provenance against the authoritative assessment sources before
using them in research. The tinyAM translations are illustrative abstractions
and proof-of-concept models. They are not exact reproductions, official stock
assessments or management advice.

## Where to look

| Material | Location |
|---|---|
| Canonical stock, assessment, assumption, input, and output records | `database/` |
| Database fields and units | `DATABASE_STRUCTURE.md` |
| Source selection, extraction, and validation rules | `PROTOCOL.md` |
| Translation, fitting, and comparison workflow | `TINYAM_TRANSLATION.md` |
| Stock-specific source checks and unresolved gaps | `source_reviews/` |
| Shared readers, converters, audits, and runners | `R/` |
| Database construction and stock model choices scripts | `scripts/translation/` |
| Readiness and aggregate run results | `results/` |

The five CSVs in `database/` are canonical. Source PDFs, data files, and model objects are cached locally under the gitignored `source_cache/` directory when available. Do not place tinyAM settings or derived run outputs in the canonical tables.

The assessment marked `is_current` is the latest accepted assessment with
detailed inputs, assumptions, and outputs recoverable for the database. It may
predate newer summary-only advice or an FSAR.

## Validate the database

From the repository root, run:

```r
source("analysis/comp_assessments/scripts/database/002_validate_database.R")
source("analysis/comp_assessments/scripts/database/003_observation_readiness.R")
```

`observation_readiness.csv` is an observation/data readiness screen: it reports input
coverage, whether observation translation succeeds, and whether
`tinyAM::check_obs()` passes. It does not load each stock recipe or call
`make_dat()`; recipe-level model readiness is checked during translation. A
ready observation object does not guarantee model convergence.

## Fit one assessment

Open `scripts/translation/run_stock.R`, set `assessment_id`, then source it. It uses the working-tree database and leaves the source record, observations, settings, fit, reference, audit, diagnostics, and comparison summary in the workspace. Use `fit = FALSE` through `run_assessment()` to inspect a translation first.

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
x <- run_assessment("dfo_cod_2j3kl_2025")
```

## Fit all applicable assessments

Source `scripts/translation/run_all.R` for the reproducible batch. It reads the committed database snapshot, uses four `future::multisession` workers through `furrr`, and writes aggregate results once in the parent process. Assessments without stock recipes and individual fit failures remain visible in the returned diagnostics. Use `parallel = FALSE` when debugging.

The batch includes every assessment marked `is_current`, including a current
accepted benchmark that has not yet been applied for advice.

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE
)
```

## Committed results and local cache

`results/` retains `observation_readiness.csv`, `fit_diagnostics.csv`, `comparison_summary.csv`, a compact `sensitivity_summary.csv`, and `northern_cod_indices.csv` for the Northern cod integration trials. Fitted objects and dashboards are generated on request and may be saved under the gitignored `results/cache/<assessment_id>/` directory with `cache = TRUE`. They are not committed by the normal workflows.

Comparison summaries separate scale agreement (absolute and percent differences)
from trajectory agreement (`trend_correlation`). Each row records its unit,
definition and matching status. Native aggregates are retained for dashboard
context; common-age calculations used in the summary are labelled separately.
See `TINYAM_TRANSLATION.md` for the matching rules and their limitations.
