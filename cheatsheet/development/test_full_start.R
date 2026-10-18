pkgload::load_all("C:/Users/paulr/Documents/tinyAM", quiet = TRUE)
library(tinyAM)
cache <- "C:/Users/paulr/Documents/tinyAM/cheatsheet/cache"
dat <- readRDS(file.path(cache, "data.rds"))
anchor <- log(max(tapply(dat$catch$obs, dat$catch$year, sum, na.rm = TRUE)))
starts <- list(log_r0 = anchor, log_r = setNames(rep(anchor, 41), 1984:2024),
                rec_beta = c("(Intercept)" = anchor))
warnings <- character()
elapsed <- system.time(fit <- withCallingHandlers(fit_tam(
  dat, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "iid", mu_form = ~1),
  F_settings = list(process = "iid", mu_form = ~1),
  M_settings = list(process = "off", mu_supplied = ~I(.3)),
  catch_settings = list(sd_form = ~1),
  index_settings = list(sd_form = ~1, q_form = ~q_block),
  start_par = starts, silent = TRUE),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }))[["elapsed"]]
saveRDS(list(value = fit, warnings = warnings, elapsed_seconds = elapsed, starts = starts),
        file.path(cache, "full_rec_start_baseline.rds"))
cat("full_rec_start_baseline converged", fit$is_converged, "optimizer", fit$opt$convergence,
      "pdHess", fit$sdrep$pdHess, "seconds", elapsed, "\n")
cat(warnings, sep = "\n")
