
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
sapply(names(retros), function(nm) (retros[[nm]]$mohns_rho |> subset(metric == "recruitment"))$rho) |>
  sort()

N_rho <- lapply(names(retros), function(nm) (retros[[nm]]$mohns_rho |> subset(metric == "N")))
names(N_rho) <- names(retros)
stack_list(N_rho)


