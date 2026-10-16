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

The assessment marked `is_current` is the latest official assessment with
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
`tinyAM::check_tam_data()` passes. It does not load each stock recipe or call
`prepare_tam()`; recipe-level model readiness is checked during translation. A
ready observation object does not guarantee model convergence.

## Fit one assessment

Open `scripts/translation/run_stock.R`, set `assessment_id`, then source it. It
loads the working-tree database, fits the selected stock script and opens a
temporary dashboard comparing tinyAM with the Official assessment. Observations,
years, ages, fit, background and the `assessment` comparison context remain in
your workspace. Each stock script retains its visible `fit_tam()` call. Observations are a
`tam_data` list; population formula columns belong in `obs$pop_covariates`.
Catch/index covariates describe observation SD and q, while weight, maturity
and spawning timing stay in their biological tables.

Edit and rerun the stock script's Model section, or update the current fit. Then
rerun the Preview line in `run_stock.R`:

```r
fit <- update(fit, catch_settings = list(sd_form = ~ age + I(age^2)))
preview <- view_assessment(fit, assessment)
preview$summary
```

The preview uses the supplied fit without fitting again, rebuilds the Official
reference on its grid, and saves nothing into the repository. Its temporary HTML
opens in the browser. Use `open_file = FALSE` to render without opening, or
`output_file = "experiment.html"` to retain a dashboard. A non-converged fit with
usable predictions can be previewed with a warning; missing uncertainty remains
missing. Inspect `check_tam(fit)` before interpreting such a fit. If you change
the age range, revise incompatible mappings in `assessment$comparison` explicitly.

For diagnostics, references and an optional dashboard, use the shared runner.
It refits the stock script rather than previewing an edited fit. `fit = FALSE` prepares
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
official benchmark that has not yet been applied for advice.

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

Pacific expansion source reviews: [BSAI Pacific ocean perch](source_reviews/afsc_pop_bsai.md),
[BSAI Atka mackerel](source_reviews/afsc_atka_mackerel_bsai.md).
[Aleutian Islands pollock](source_reviews/afsc_pollock_ai.md).
[Alaska plaice](source_reviews/afsc_alaska_plaice_bsai.md).
[Eastern Bering Sea yellowfin sole](source_reviews/afsc_yellowfin_sole_ebs.md).
[Eastern Bering Sea northern rock sole](source_reviews/afsc_northern_rock_sole_ebs.md).
[Joint coastal Pacific hake](source_reviews/joint_hake_coastal.md).
[Eastern Bering Sea Pacific cod](source_reviews/afsc_cod_ebs.md).
[Strait of Georgia Pacific herring](source_reviews/dfo_herring_sog.md) is curated,
but translation is blocked by its aggregate spawning-biomass survey without
survey ages; it has not received a tinyAM fitting attempt.
[British Columbia Pacific ocean perch](source_reviews/dfo_pop_bc.md) has a
pooled female translation, with weakly identified M and an unresolved native
spawning-output definition flagged for review.
[Aleutian Islands Pacific cod](source_reviews/afsc_cod_ai.md) remains blocked
pending native inputs matching the final production assessment.

`results/` retains `observation_readiness.csv`, `fit_diagnostics.csv`, `comparison_summary.csv`, a compact `sensitivity_summary.csv`, and `northern_cod_indices.csv` for the Northern cod integration trials. Fitted objects and dashboards are generated on request and may be saved under the gitignored `results/cache/<assessment_id>/` directory with `cache = TRUE`. They are not committed by the normal workflows.

Comparison summaries separate scale agreement (absolute and percent differences)
from trajectory agreement (`trend_correlation`). Each row records its unit,
definition and matching status. Native aggregates are retained for dashboard
context; common-age calculations used in the summary are labelled separately.
See `TINYAM_TRANSLATION.md` for the matching rules and their limitations.

## Reviewing model choices

`results/translation_review.csv` retains the historical candidate review and
its pinned baseline revision. The retired review script remains in Git history;
current model choices are visible in each stock script's `fit_tam()` call.
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
