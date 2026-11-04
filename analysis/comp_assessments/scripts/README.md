# Assessment scripts

Run these scripts from the tinyAM repository root in RStudio.

1. `update_database_cache.R`: manually copy the five clean, committed SAIO CSVs into `cache/saio/<commit>/`. Set `saio_repo` first. The CSVs are parsed once into `database.rds`; the tracked manifest is updated only after both steps finish; SAIO is never changed or fetched.
2. `run_stock.R`: choose a stock and tier, fit its literal recipe, then preview the supplied fit independently.
3. `run_all.R`: read the snapshot once, configure `future::plan(future::multisession, workers = 10)`, fit/checkpoint the batch, then render dashboards and overview figures in separate sections.
4. `observation_readiness.R`: audit the same supplied `database`; when run alone, read the selected local snapshot.
5. `stats_plots.R` and `trend_plots.R`: rebuild overview figures from results matching the selected source revision.

`run_assessment(..., fit = FALSE)` prepares inputs without final or preliminary fits. The batch functions honor your existing future plan; choose `future::sequential` for troubleshooting. `run_all.R` selects a 20-minute per-assessment deadline; use `timeout = Inf` for unrestricted fitting. A finite deadline requires `callr` and records a timeout without retuning the model.

With `cache = TRUE`, each completed attempt is saved immediately under `results/<tier>/cache/<SAIO commit>/<assessment_id>/`, along with the fit when available and stage logs. Interrupted serialization preserves the previous complete checkpoint. `cache` does not reuse or skip fits. Only an intentional full batch should set `save_results = TRUE` to refresh aggregate CSVs.

Run `render_assessment_dashboards(runs)` after fitting. It records rendered, troubleshooting, skipped and failed outcomes in `dashboard_manifest.csv` and rebuilds the tier's index. It never fits or changes numerical eligibility. Troubleshooting dashboards may lack reliable estimates or uncertainty; inspect `check_tam()` and the diagnostic table.

For retrospectives or hindcasts, configure the same future plan and call `fit_retro(fit, folds = 5)` or `fit_hindcast(fit, folds = 5)`, then `vis_tam(result)`. Close workers with `future::plan(future::sequential)` when finished. See the [workflow](../README.md) for catch conditioning and score interpretation.

All generated snapshots, fits, dashboards and figures stay local and ignored. Canonical source curation, extraction and validation are maintained in [SAIO](https://github.com/PaulRegular/saio).
