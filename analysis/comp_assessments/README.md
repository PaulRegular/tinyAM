# Comparing stock assessments with tinyAM

Source data are maintained in the independent **SAIO — Stock Assessment Inputs and Outputs Database**. This folder contains tinyAM translations and comparisons.

## Curated and preliminary translations

- scripts/curated/: the paper working set, built from the 35 scripts retained at 78f3d347331b0e770d55a855de4c5b90307fd400 and subsequent reviews. Membership does not imply convergence or scientific equivalence.
- scripts/preliminary/: additional translations and changed alternatives from comp-assessments. They have not been promoted into the paper working set.
- translations.csv: the explicit registry, including each script's originating commit. An assessment can have a curated version and a preliminary alternative.

Curated runs are the default. Preliminary runs require an explicit tier. Extraction completeness, translation review and numerical convergence are separate. Move a script and update the registry only after reviewing its sources, assumptions and comparisons. If both tiers contain the same assessment, selecting its ID without a tier uses the curated version.

The first 20-stock expansion beyond Charbonneau–Keith is tracked in
[expansion_batch_01.csv](expansion_batch_01.csv), using the selected stock units
in SAIO. The registry shows which recipes from that batch have been promoted.
The registry's `source_commit`
records the origin of migrated scripts; new recipes leave it blank and use Git
history and run diagnostics for their model revision. The SAIO data revision
is recorded separately.

The batch register also records translation limitations. `source_readiness_blocked`
means the reviewed source record does not yet support a defensible fit; it does
not establish that data can never be recovered. `reconstruction_pending` means
native observations are present but their age/length, sex, area or removal
definitions still need a reviewed conversion. Neither status counts as a fit.

## 1. Read local SAIO data

Both tiers read the five database tables directly from the current commit of your local SAIO checkout. Update or commit data in SAIO before rerunning an analysis. [SAIO source curation and validation](https://github.com/PaulRegular/saio)

```r
source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
```

Run from the tinyAM repository root. By default the reader looks for `saio` next to the primary tinyAM checkout, including when you work from a linked worktree. For another checkout, set `SAIO_REPO` or pass its path as `saio_repo` to `read_database()`. The reader uses committed Git data and locally available Git LFS objects. It ignores uncommitted edits and performs no network requests or checkout changes. Parsed tables are reused within the R session for the same repository and commit.

To reproduce an earlier analysis, pass the recorded SHA as `revision` to `read_database()`. That commit and its table data must be available locally. If Git LFS objects are missing, run `git lfs fetch origin <sha>` from the SAIO checkout. Pass the returned `database` to the fitting and overview scripts to use the same revision throughout.

One preliminary English Channel plaice recipe retains its documented native-SAM starting-value retry. That optional retry needs `source_cache/ices_plaice_7d_2026/fit_local_locked.RData`, copied manually from the reviewed source artifact.

## 2. Fit assessments

For one stock, open `scripts/run_stock.R`, select `translation_tier` and `assessment_id`, and source its Fit section. The literal model call remains in the registered recipe. The workspace exposes its inputs, fit and comparison context. Rerun the Preview line after editing the model; previewing never refits it.

For a batch, open `scripts/run_all.R`. Its stages are deliberately separate:

```r
future::plan(future::multisession, workers = 10)
runs <- run_assessments(database = database, tier = "curated",
  cache = TRUE, save_results = TRUE, timeout = 20 * 60)
```

Use `future::plan(future::sequential)` to troubleshoot. `run_assessments()` honors the caller's plan and never changes it. Its default `timeout = Inf` is unbounded; the batch script visibly selects a 20-minute deadline. A finite limit uses a supervised R subprocess so native optimizer code can be stopped. Every failed or timed-out attempt remains in diagnostics; model settings and starting-value policies are unchanged.

Each completed assessment is checkpointed immediately when caching is enabled. Checkpoints are replaced only after serialization completes. `cache` saves results; it never skips fitting. Aggregate CSVs are written only with `save_results = TRUE`, after workers finish. Targeted trials must not overwrite the canonical batch summaries.

## 3. Render and inspect dashboards

```r
dashboards <- render_assessment_dashboards(runs)
```

This stage uses returned fits and their original source context. It does not reread SAIO, refit models or calculate OSA residuals. Nonconverged fits with usable predictions may produce clearly identified **troubleshooting** dashboards; they do not become numerically eligible comparisons. Missing fits and render errors remain in `dashboard_manifest.csv`. A rendering failure never loses the saved fit.

Open `results/<tier>/index.html`. Inspect Fisheries/Surveys residuals, Tables → Diagnostics and Model structure before interpreting trends. Rerun dashboard generation independently after fixing a plotting problem. For an edited individual fit, use `view_assessment(fit, assessment)`; explicit mismatches between fit and source revisions remain errors.

The overview scripts `scripts/stats_plots.R` and `scripts/trend_plots.R` can also run separately. They reject result tables from another data revision.

## Retrospectives and hindcasts

The same caller-owned future plan controls package evaluation functions:

```r
future::plan(future::multisession, workers = 10)
retro <- fit_retro(fit, folds = 5)
hind <- fit_hindcast(fit, folds = 5, start_from_fit = TRUE)
print(hind)
vis_tam(hind)
future::plan(future::sequential)
```

The hindcast default conditions each forecast on actual aggregate catch, while holding out survey observations and catch composition. Catch-total agreement is therefore not forecast skill. `proj_settings = list(F_mult = 1)` is the status-quo-F alternative. Compare validation coverage before comparing scores; failed folds remain in `fold_diagnostics`. See `?fit_retro` for Mohn's rho, log-RMSE and MASE definitions and multi-model calls.

## Results and translation guidance

Generated results remain local and ignored under `results/curated/` and `results/preliminary/`. Fits, checkpoints and dashboards are separated by `cache/<SAIO commit>/<assessment_id>/`. Batch checkpoints share that revision folder. Cached fits record both the SAIO commit and the tinyAM revision. Existing historical results keep their original provenance.

The scientific translation instructions remain in [TINYAM_TRANSLATION.md](TINYAM_TRANSLATION.md). Stock-specific assumptions and inherited reviews are in [translation_notes/](translation_notes/README.md); official source evidence is maintained in SAIO. The migration record is in [migration/](migration/README.md). No operational assessment or study hypothesis is changed by this workflow.
