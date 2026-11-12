# Assessment scripts

Run these scripts from the tinyAM repository root in RStudio.

The scripts read committed tables directly from a local SAIO checkout. They use its current commit by default. Configure the checkout path with `SAIO_REPO`, or pass `saio_repo` and optionally `revision` to `read_database()` when supplying the database explicitly. [Analysis workflow](../README.md#1-read-local-saio-data)

1. `run_stock.R` chooses a stock and tier, fits its recipe, and previews the supplied fit independently.
2. `run_all.R` reads the database once, configures `future::plan(future::multisession, workers = 10)`, fits and checkpoints the batch, then renders dashboards and overview figures in separate sections.
3. `observation_readiness.R` audits the supplied `database` or reads the local SAIO checkout when run alone.
4. `stats_plots.R` and `trend_plots.R` rebuild overview figures from results matching the selected source revision.

`run_assessment(..., fit = FALSE)` prepares inputs without final or preliminary fits. The batch functions honor your existing future plan; choose `future::sequential` for troubleshooting. `run_all.R` selects a 20-minute per-assessment deadline; use `timeout = Inf` for unrestricted fitting. A finite deadline requires `callr` and records a timeout without retuning the model.

With `cache = TRUE`, each completed attempt is saved immediately under `results/<tier>/cache/<SAIO commit>/<assessment_id>/`, along with the fit when available and stage logs. Interrupted serialization preserves the previous complete checkpoint. `cache` does not reuse or skip fits. Only an intentional full batch should set `save_results = TRUE` to refresh aggregate CSVs.

Run `render_assessment_dashboards(runs)` after fitting. It records rendered, troubleshooting, skipped and failed outcomes in `dashboard_manifest.csv` and rebuilds the tier's index. It never fits or changes numerical eligibility. Troubleshooting dashboards may lack reliable estimates or uncertainty; inspect `check_tam()` and the diagnostic table.

For retrospectives or hindcasts, configure the same future plan and call `fit_retro(fit, folds = 5)` or `fit_hindcast(fit, folds = 5)`, then `vis_tam(result)`. Close workers with `future::plan(future::sequential)` when finished. See the [workflow](../README.md) for catch conditioning and score interpretation.

Generated fits, dashboards and figures stay local and ignored. [SAIO source curation, extraction and validation](https://github.com/PaulRegular/saio)
