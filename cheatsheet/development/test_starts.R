pkgload::load_all("C:/Users/paulr/Documents/tinyAM", quiet = TRUE)
library(tinyAM)
cache <- "C:/Users/paulr/Documents/tinyAM/cheatsheet/cache"
dat <- readRDS(file.path(cache, "data.rds"))
settings <- list(years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "iid", mu_form = ~1),
  M_settings = list(process = "off", mu_supplied = ~I(.3)),
  catch_settings = list(sd_form = ~1),
  index_settings = list(sd_form = ~1, q_form = ~q_block), silent = TRUE)
attempt <- function(label, args) {
  warnings <- character()
  message(label, " started")
  elapsed <- system.time(fit <- withCallingHandlers(do.call(fit_tam, args),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }))[["elapsed"]]
  saveRDS(list(value = fit, warnings = warnings, elapsed_seconds = elapsed), file.path(cache, paste0(label, ".rds")))
  cat(label, "converged", fit$is_converged, "optimizer", fit$opt$convergence,
      "pdHess", fit$sdrep$pdHess, "seconds", elapsed, "\n")
  cat(warnings, sep = "\n")
}
attempt("original_cod_baseline", c(list(data = cod_obs), settings))
anchor <- log(max(tapply(dat$catch$obs, dat$catch$year, sum, na.rm = TRUE)))
cat("Catch-derived starting log_r0", anchor, "\n")
attempt("scaled_start_baseline", c(list(data = dat, start_par = list(log_r0 = anchor)), settings))
