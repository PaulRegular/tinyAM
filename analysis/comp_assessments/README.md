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
| Source imports and validators | `scripts/database/` |
| Visible stock model calls and run drivers | `scripts/translation/` |
| Remaining model capabilities and priorities | `CAPABILITY_GAPS.md` |
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
`prepare_tam()`; recipe-level model readiness is checked during translation. A
ready observation object does not guarantee model convergence.

## Fit one assessment

Open `scripts/translation/run_stock.R`, set `assessment_id`, then source it. It loads the working-tree database and sources the selected stock script. The observations, years, ages, literal `fit_tam()` call, fit, background and comparison metadata are visible in your workspace. You can then rerun sections of that stock script interactively.

For diagnostics, references and an optional dashboard, use the shared runner.
It sources the same script in an isolated environment. `fit = FALSE` prepares
observations without running final or preliminary fits. Settings are available
in `x$fit$dat`; the result does not duplicate them in a settings list.

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
x <- run_assessment("dfo_cod_2j3kl_2025", dashboard = TRUE, cache = TRUE)
check_tam(x$fit)
```

## Fit all applicable assessments

Source `scripts/translation/run_all.R` for the reproducible batch. It reads the committed database snapshot, uses `future::multisession` through `furrr`, and writes aggregate results once in the parent process. The script currently uses one worker; increase `workers` when memory permits. Assessments without stock recipes and individual fit failures remain visible in the returned diagnostics. Use `parallel = FALSE` when debugging.

The batch includes every assessment marked `is_current`, including a current
accepted benchmark that has not yet been applied for advice.

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 1,
  save_results = TRUE,
  cache = TRUE,
  dashboard = TRUE
)
```

## Committed results and local cache

`results/` retains `observation_readiness.csv`, `fit_diagnostics.csv`, `comparison_summary.csv`, a compact `sensitivity_summary.csv`, and `northern_cod_indices.csv` for the Northern cod integration trials. Fitted objects and dashboards are generated on request and may be saved under the gitignored `results/cache/<assessment_id>/` directory with `cache = TRUE`. They are not committed by the normal workflows.

Comparison summaries separate scale agreement (absolute and percent differences)
from trajectory agreement (`trend_correlation`). Each row records its unit,
definition and matching status. Native aggregates are retained for dashboard
context; common-age calculations used in the summary are labelled separately.
See `TINYAM_TRANSLATION.md` for the matching rules and their limitations.

## Reviewing model choices

`scripts/translation/review_translations.R` records source-led candidates and
the pinned baseline revision. Source it after the runner, then call
`review_translations()`. Every attempt, including failures and advisories, is
cached locally; `results/translation_review.csv` is the compact tracked record.
The review includes all current detailed assessments and records blockers for
those without recipes. Baselines remain when changes weaken diagnostics or
trade one comparison improvement for another loss. Numerical convergence and
statistical support are reviewed separately; fixed-M agreement is not evidence
of improvement. See `CAPABILITY_GAPS.md` for the remaining limitations.

The review table has one row per attempted model. `*_percent` is the mean
absolute percent difference; `*_terminal_percent` is the signed terminal-year
difference; `*_trend` is trajectory correlation. The matching status is retained
alongside each metric. Detailed definitions and units are in the aggregate
comparison summary and cached reference objects. Diagnostic failures remain
visible; their comparison values are not grounds for selecting that model.
Pinned baseline models and candidates use the same current comparison
definitions, including corrections to reported age-group interpretation.
