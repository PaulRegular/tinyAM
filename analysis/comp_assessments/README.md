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
| Charbonneau–Keith Master catalogue crosswalk | `source_reviews/charbonneau_keith_master_catalogue.md` |
| Gulf of Maine haddock source review | [source_reviews/nefsc_haddock_gulf_of_maine.md](source_reviews/nefsc_haddock_gulf_of_maine.md) |
| Georges Bank winter flounder | [source review](source_reviews/nefsc_winter_flounder_georges_bank.md) · [tinyAM translation](scripts/translation/stocks/nefsc_winter_flounder_gb_2022_vpa.R) |
| Eastern English Channel plaice | [source review](source_reviews/ices_plaice_english_channel_7d_2026.md) · [tinyAM translation](scripts/translation/stocks/ices_plaice_english_channel_7d_2026.R) |
| Irish Sea plaice | [source review](source_reviews/ices_plaice_irish_sea_7a_2026.md) · [tinyAM translation](scripts/translation/stocks/ices_plaice_irish_sea_7a_2026.R) |
| Scup source review | [source_reviews/ck_nefsc_saw_chrysops_usatl.md](source_reviews/ck_nefsc_saw_chrysops_usatl.md) |
| Central Baltic herring | [source review](source_reviews/ck_ices_wgbfas_harengus_bs_25_29_32.md) · [tinyAM translation](scripts/translation/stocks/ices_herring_central_baltic_2026.R) |
| Bay of Biscay sardine | [source review](source_reviews/ck_ices_wghansa_pilchardus_sp8abd.md) |
| Southern New England/Mid-Atlantic yellowtail flounder | [source review](source_reviews/ck_nefsc_garmiii_ferruginea_sne_ma.md) |
| Georges Bank yellowtail flounder | [source review and production-data gap](source_reviews/ck_dfo_nefsc_ferruginea_gb.md) |
| Northeast Arctic beaked redfish | [source review and remaining input gaps](source_reviews/ck_ices_afwg_mentella_deep_1_2.md) |
| 4VWX Atlantic herring | [component indicators and production-input gap](source_reviews/ck_dfo_mar_harengus_4vwx.md) |
| Atlantic striped bass | [base-run outputs and missing production inputs](source_reviews/ck_nefsc_garmiii_saxatilis_gom_chatt.md) |
| Eastern Georges Bank cod | [selected framework, later summaries and missing age inputs](source_reviews/ck_nefsc_trac_cert_morhua_gb.md) |
| North Sea sandeel | [catalogue identity](source_reviews/ck_ices_hawg_dubius_ns.md); [Area 1r](source_reviews/ices_sandeel_area_1r.md); [Area 2r](source_reviews/ices_sandeel_area_2r.md); [Area 3r](source_reviews/ices_sandeel_area_3r.md); [Area 4](source_reviews/ices_sandeel_area_4.md) production assessments |
| Celtic Sea and Bristol Channel plaice | [current indicator assessment and age-based translation gap](source_reviews/ices_plaice_celtic_sea_7fg.md) |
| Northeast Arctic golden redfish | [2026 SPiCT states and age-based translation gap](source_reviews/ices_golden_redfish_northeast_arctic.md) |
| Western English Channel plaice | [2025 chr rule, surviving discards and translation gap](source_reviews/ices_plaice_western_channel_7e.md) |
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

Source `scripts/translation/run_all.R` for the reproducible batch. It reads the committed database snapshot, uses `future::multisession` through `furrr`, and writes aggregate results once in the parent process. The script currently uses four workers; reduce `workers` when memory is limited. Assessments without stock recipes and individual fit failures remain visible in the returned diagnostics. Use `parallel = FALSE` when debugging.

The batch includes every assessment marked `is_current`, including a current
official benchmark that has not yet been applied for advice.

```r
pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
runs <- run_assessments(
  database = read_committed_database(),
  parallel = TRUE,
  workers = 4,
  save_results = TRUE,
  cache = TRUE,
  dashboard = FALSE
)
```

## Committed results and local cache

`results/` retains observation readiness, fit diagnostics, comparison summaries,
annual population comparisons and five overview figures. These are compact
snapshots of the complete batch; individual stock trials do not overwrite them.
Fitted objects, dashboards, warnings and failures stay in the gitignored
`results/cache/` directory. Authoritative downloaded sources stay in
`source_cache/`. Neither directory is pruned by routine repository tidying.

Current source availability and translation blockers are recorded in each
stock's `source_reviews/` note. A stock recipe may deliberately reject fitting
when the necessary source data or comparison definitions remain unavailable.

Translations and comparison displays use thousand fish, average weight in kg,
and biomass in tonnes. Count-rate surveys retain their effort denominator;
relative indices retain their source scale. See the units section of
[TINYAM_TRANSLATION.md](TINYAM_TRANSLATION.md#translation-units).
The complete `run_all.R` workflow saves returned fits before making dashboards,
including exploratory previews of usable non-converged fits. Rendering failures
are recorded locally. Figure fallbacks use this complete batch rather than
older stock caches; aggregate comparison eligibility is unchanged.

Comparison summaries separate scale agreement (absolute and percent differences)
from trajectory agreement (`trend_correlation`). Each row records its unit,
definition and matching status. Native aggregates are retained for dashboard
context; common-age calculations used in the summary are labelled separately.
See `TINYAM_TRANSLATION.md` for the matching rules and their limitations.

## Reviewing model choices

Current model choices are visible in each stock script's `fit_tam()` call.
Superseded trials and review reports remain in Git history. Source notes retain
current rationale, extraction definitions and limitations rather than a model
development diary. See `CAPABILITY_GAPS.md` for supported approximations and
remaining capabilities.

Retain the baseline when a candidate weakens diagnostics or produces unclear
comparison trade-offs. Numerical convergence and statistical support are
separate; fixed-M agreement is not evidence of successful estimation.

The agreement heatmap includes exploratory, non-converged fits where reporting
is available. Check `results/fit_diagnostics.csv` and the comparison definitions
before interpreting its colours. Blank or grey cells do not indicate agreement.

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
