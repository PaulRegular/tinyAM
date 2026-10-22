library(tinyAM)

models <- readRDS("analysis/comp_retro/outputs/001_models.rds")
future::plan(future::multisession, workers = 2)

retros <- lapply(models, function(fit) {
  fit_hindcast(fit, folds = 20, start_from_fit = TRUE, grad_tol = 1e-2)
})
saveRDS(retros, "analysis/comp_retro/outputs/002_retros.rds")
for (model in names(retros)) {
  cat("\n", model, "\n")
  print(retros[[model]])
}
cli::cli_inform("Before comparing these separate runs, check forecast_pairs for matching validation cases and units.")

for (model in names(retros)) {
  fitted <- Filter(Negate(is.null), retros[[model]]$fits[[1L]])
  if (length(fitted)) vis_tam(
    model_list = fitted,
    output_file = paste0("analysis/comp_retro/outputs/002_retro_", model, ".html"),
    open_file = FALSE
  )
}

collect <- function(table) {
  do.call(rbind, Map(function(result, model) {
    d <- result[[table]]
    if (nrow(d)) d$model <- model
    d
  }, retros, names(retros)))
}
scores <- collect("scores")
subset(scores, type == "index" & level == "overall" & is.na(horizon))
rhos <- collect("mohns_rho")
subset(rhos, metric == "ssb")
cols <- tam_pal(length(retros))
age_metrics <- unique(rhos$metric[is.finite(rhos$age)])

rhos |>
  subset(metric %in% age_metrics) |>
  plotly::plot_ly(x = ~age, y = ~rho, color = ~model, colors = cols, frame = ~metric) |>
  plotly::add_bars()

rhos |>
  subset(!metric %in% age_metrics) |>
  plotly::plot_ly(x = ~metric, y = ~rho, color = ~model, colors = cols) |>
  plotly::add_bars()

subset(rhos, model == "N_iid_F_rw" & metric %in% c("N", "F"))
