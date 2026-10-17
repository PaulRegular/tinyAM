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
[British Columbia Pacific ocean perch](source_reviews/dfo_pop_bc.md) remains
curated but is blocked from fitting and comparison: matching numerical historical
coastwide population outputs have not been recovered. Partial regional biomass
draws are now documented in its source review. Earlier exploratory fits are archived;
constant M and boundary biomass alone do not support the intended comparison.
[British Columbia petrale sole](source_reviews/dfo_petrale_bc.md) has a partial
official 2024 summary/output record, including annual SSB posterior summaries.
Original observations and biology remain missing; it is not a current detailed
record or a completed translation.
[British Columbia Canary Rockfish](source_reviews/dfo_canary_bc.md) has verified
native 2022 inputs, fixed biological functions and annual female SSB posterior
summaries; official N/F-at-age and annual recruitment remain unavailable.
[Aleutian Islands Pacific cod](source_reviews/afsc_cod_ai.md) remains blocked
pending native inputs matching the final production assessment.
[BSAI arrowtooth flounder](source_reviews/afsc_arrowtooth_flounder_bsai.md)
has recovered sex-specific states and biology, but missing likelihood
compositions, an age–length key and timing block a defensible fit.

`results/` retains `observation_readiness.csv`, `fit_diagnostics.csv`, `comparison_summary.csv`, a compact `sensitivity_summary.csv`, and `northern_cod_indices.csv` for the Northern cod integration trials. Fitted objects and dashboards are generated on request and may be saved under the gitignored `results/cache/<assessment_id>/` directory with `cache = TRUE`. They are not committed by the normal workflows.

Comparison summaries separate scale agreement (absolute and percent differences)
from trajectory agreement (`trend_correlation`). Each row records its unit,
definition and matching status. Native aggregates are retained for dashboard
context; common-age calculations used in the summary are labelled separately.
See `TINYAM_TRANSLATION.md` for the matching rules and their limitations.

## Reviewing model choices

The historical candidate review and its pinned baseline revision remain in Git
history (the parent of cleanup commit `4830fb1`, at
`analysis/comp_assessments/results/translation_review.csv`). The retired review script also remains in Git history;
current model choices are visible in each stock script's `fit_tam()` call.
That historical review covered the 23 then-current detailed assessments,
including blockers for those without recipes. The expanded-library review is
in [PACIFIC_REVIEW.md](PACIFIC_REVIEW.md), including the source/model findings
and a link to the archived one-off trials. Experimental trial runners and
duplicate review tables are retained in Git history rather than the routine workflow.
Delivery and validation details are in [PACIFIC_EXPANSION.md](PACIFIC_EXPANSION.md).
Baselines remain when changes weaken diagnostics or
trade one comparison improvement for another loss. Numerical convergence and
statistical support are reviewed separately; fixed-M agreement is not evidence
of improvement. See `CAPABILITY_GAPS.md` for the remaining limitations.

The agreement heatmap includes exploratory, non-converged fits where reporting
is available. Check the numerical outcomes and comparison definitions in
`PACIFIC_REVIEW.md` before interpreting its colours; blank or grey cells do not
indicate agreement. Its median differences complement the mean and terminal
differences used for review flags.

Both figure scripts share `R/plot_stock_order.R`: ocean first, then
phylogenetic groups, scientific species names and stock labels. Higher groups
follow [Betancur-R. et al. (2017)](https://doi.org/10.1186/s12862-017-0958-3);
alphabetical ties within groups do not imply a fully resolved species tree.
Trend facets follow the heatmap from top to bottom, reading facets left to
right and top to bottom. Add new genera to the shared ordering before plotting.
Trend figures use a seven-column, four-row layout. The stock name appears
above a smaller grey ocean label, rendered with `ggtext`.

Current comparison tables retain matching status, definitions, units and
mean/median/terminal differences. Annual trajectories in
`results/population_timeseries.csv` must agree with those summaries. Diagnostic
failures remain visible; their comparison values are not grounds for selecting
a model. Cached exploratory fits are retained without inventing missing
official estimates.
