
library(tinyAM)
library(furrr)

# source("analysis/comp_retro/001_models.R")
models <- readRDS("analysis/comp_retro/outputs/001_models.rds")

future::plan(future::multisession, workers = 10)

retros <- lapply(names(models), function(nm) {
  try(fit_retro(models[[nm]], folds = 20, hindcast = TRUE, grad_tol = 1e-2))
})
names(retros) <- names(models)

saveRDS(retros, file = "analysis/comp_retro/outputs/002_retros.rds")

lapply(names(retros), function(nm) {
  if (!inherits(retros[[nm]], "try-error")) {
    vis_tam(
      retros[[nm]]$fits,
      interval = 0.9,
      output_file = paste0("analysis/comp_retro/outputs/002_retro_", nm, ".html"),
      open_file = FALSE
    )
  }
})

sapply(names(retros), function(nm) retros[[nm]]$hindcast_rmse) |>
  sort()

sapply(names(retros), function(nm) (retros[[nm]]$mohns_rho |> subset(metric == "ssb"))$rho) |>
  sort()

rhos <- lapply(names(retros), function(nm) (retros[[nm]]$mohns_rho))
names(rhos) <- names(retros)
rhos <- stack_list(rhos, label_type = "factor")
cols <- tam_pal(length(retros))

counts <- table(rhos$metric)
age_metrics <- names(which(counts == max(counts)))

rhos |>
  subset(metric %in% age_metrics) |>
  plot_ly(x = ~age, y = ~rho, color = ~model, colors = cols,
          frame = ~metric) |>
  add_bars()

rhos |>
  subset(!metric %in% age_metrics) |>
  plot_ly(x = ~metric, y = ~rho, color = ~model, colors = cols) |>
  add_bars()

retros$N_iid_F_rw$mohns_rho |>
  subset(metric %in% c("N", "F"))
