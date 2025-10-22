
library(tinyAM)

# source("analysis/comp_retro/001_models.R")
models <- readRDS("analysis/comp_retro/outputs/001_models.rds")

future::plan(future::multisession, workers = 10)

retros <- lapply(names(models), function(nm) {
  try(fit_retro(models[[nm]], folds = 20, hindcast = TRUE, grad_tol = 1e-2))
})
names(retros) <- names(models)

saveRDS(retros, file = "analysis/comp_retro/outputs/002_retros.rds")

for (nm in names(retros)) {
  if (!inherits(retros[[nm]], "try-error")) {
    vis_tam(
      retros[[nm]]$fits,
      interval = 0.9,
      output_file = paste0("analysis/comp_retro/outputs/002_retro_", nm, ".html"),
      open_file = FALSE
    )
  }
}



