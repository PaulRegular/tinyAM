pkgload::load_all("C:/Users/paulr/Documents/tinyAM", quiet = TRUE)
library(tinyAM)
repo <- "C:/Users/paulr/Documents/tinyAM"
cache <- file.path(repo, "cheatsheet", "cache")
attempt <- function(label, code) {
  warnings <- character()
  error <- NULL
  message(label, " started")
  elapsed <- system.time(fit <- withCallingHandlers(tryCatch(code(), error = function(e) {
      error <<- conditionMessage(e); NULL
    }), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }))[["elapsed"]]
  saveRDS(list(value = fit, warnings = warnings, error = error, elapsed_seconds = elapsed),
          file.path(cache, paste0(label, ".rds")))
  if (is.null(fit)) { cat(label, "ERROR", error, "\n"); return(NULL) }
  cat(label, "converged", fit$is_converged, "optimizer", fit$opt$convergence,
      "pdHess", fit$sdrep$pdHess, "seconds", elapsed, "\n")
  cat(warnings, sep = "\n")
  fit
}
control <- attempt("control_package_default", function() {
  source(file.path(repo, "inst", "examples", "example_fit_default.R"), local = TRUE)
  fit
})
attempt("control_no_N_rw", function() fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  F_settings = list(process = "rw", mu_form = NULL),
  M_settings = list(process = "off", mu_supplied = ~I(.3)),
  catch_settings = list(sd_form = ~1),
  index_settings = list(sd_form = ~1, q_form = ~q_block), silent = TRUE))
if (isTRUE(control$is_converged)) attempt("control_baseline_F_ar1", function() fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "ar1", mu_form = ~1),
  M_settings = list(process = "off", mu_supplied = ~I(.3)),
  catch_settings = list(sd_form = ~1),
  index_settings = list(sd_form = ~1, q_form = ~q_block), silent = TRUE))
