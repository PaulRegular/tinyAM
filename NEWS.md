# tinyAM 0.1.1.9000

## Evaluation workflows

* `fit_retro()` and `fit_hindcast()` accept multiple fits and return `tam_retro`
  objects with a concise print method, detailed fold diagnostics and retained failures.
* Hindcasts support multiple forecast horizons, biological averaging settings,
  persistence-relative MASE, survey age/aggregate scores and catch-composition scores.
* **Breaking:** forecast controls now belong in `proj_settings`; numeric vectors
  are ordered by horizon. Fitted folds are nested by model and terminal year.
  Scores are in `$scores`, replacing the former single `$hindcast_rmse` value.
* `start_par` accepts a fitted model, including when uncertainty estimation failed.

The case study and cheatsheet need user-led migration to the new hindcast
interface and nested results. They are intentionally unchanged in this work.
