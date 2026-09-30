# Run from the repository root after installing/loading this branch of tinyAM.
# The fitted SAM reference is loaded directly, never reconstructed or refitted.
library(tinyAM)
case <- "WKCOD_combined_99"
dir <- file.path("analysis", "comp_sam")
dir.create(file.path(dir, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dir, "results"), recursive = TRUE, showWarnings = FALSE)
cache <- file.path(dir, "cache", paste0(case, ".rds"))
if (!file.exists(cache)) {
  sam_fit <- stockassessment::fitfromweb(case, character.only = TRUE)
  saveRDS(sam_fit, cache)
  provenance <- data.frame(model = case,
    url = paste0("https://stockassessment.org/datadisk/stockassessment/userdirs/user3/", case, "/run/model.RData"),
    downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stockassessment_version = as.character(packageVersion("stockassessment")),
    stockassessment_sha = packageDescription("stockassessment")$RemoteSha)
  write.csv(provenance, file.path(dir, "reference_provenance.csv"), row.names = FALSE)
}
sam_fit <- readRDS(cache)
write.csv(data.frame(file = basename(cache), md5 = unname(tools::md5sum(cache))),
          file.path(dir, "reference_checksum.csv"), row.names = FALSE)
# Original maturity is absent before 1983. This common period was selected
# explicitly; neither observations nor biology are filled from fitted SAM values.
years <- 1983:2022
tam_obs <- sam_to_tam_obs(sam_fit)
tam_obs <- lapply(tam_obs, function(d) d[d$year %in% years, , drop = FALSE])
settings <- sam_to_tam_settings(sam_fit, overrides = list(years = years))
audit <- sam_to_tam_audit(sam_fit, settings)
write.csv(audit, file.path(dir, "results", "audit.csv"), row.names = FALSE)
settings_table <- data.frame(argument = names(settings),
  value = vapply(settings, function(x) paste(deparse(x), collapse = " "), character(1)))
write.csv(settings_table, file.path(dir, "results", "settings.csv"), row.names = FALSE)
for (nm in names(tam_obs)) write.csv(tam_obs[[nm]], file.path(dir, "results", paste0("input_", nm, ".csv")), row.names = FALSE)

# SAM estimates provide a reproducible starting point only. Every tinyAM
# coefficient/state remains free under the translated model.
dat <- do.call(make_dat, c(list(obs = tam_obs), settings))
start <- make_par(dat)
y <- match(years, sam_fit$data$years)
start$log_r0 <- sam_fit$pl$logN[1, y[1]]
start$log_n0[] <- sam_fit$pl$logN[-1, y[1]]
start$log_r[] <- sam_fit$pl$logN[1, y[-1]]
start$log_n[,] <- t(sam_fit$pl$logN[-1, y[-1], drop = FALSE])
fkeys <- sam_fit$conf$keyLogFsta[which(sam_fit$data$fleetTypes == 0), ]
start$log_f[,] <- t(sam_fit$pl$logF[fkeys + 1L, y, drop = FALSE])
start$log_q[] <- sam_fit$pl$logFpar[sort(unique(tam_obs$index$q_key)) + 1L]
start$log_sd_catch[] <- sam_fit$pl$logSdLogObs[sort(unique(tam_obs$catch$sd_key)) + 1L]
start$log_sd_index[] <- sam_fit$pl$logSdLogObs[sort(unique(tam_obs$index$sd_key)) + 1L]
start$log_sd_f <- sam_fit$pl$logSdLogFsta[1]
start$log_sd_r <- sam_fit$pl$logSdLogN[1]
start$log_sd_n <- sam_fit$pl$logSdLogN[2]
fit_file <- file.path(dir, "results", "tinyAM_fit.rds")
timing <- system.time(tam_fit <- do.call(fit_tam, c(list(obs = tam_obs, start_par = start, silent = TRUE), settings)))
saveRDS(tam_fit, fit_file)
# Explicit common-period reporting view; the complete original SAM fit is retained.
sam_comparison <- sam_to_tam_comparison(sam_fit)
sam_comparison$dat$years <- years
sam_comparison$dat$is_proj <- rep(FALSE, length(years))
sam_comparison$dat$obs <- tam_obs
models <- list(SAM = sam_comparison, tinyAM = tam_fit)

grad <- tam_fit$obj$gr(tam_fit$opt$par)
diagnostics <- data.frame(model = "tinyAM", is_converged = tam_fit$is_converged,
  convergence_code = tam_fit$opt$convergence, optimizer_message = tam_fit$opt$message,
  objective = tam_fit$opt$objective, max_abs_gradient = max(abs(grad)),
  sdreport_success = all(is.finite(tam_fit$sdrep$sd)),
  positive_definite_hessian = tam_fit$sdrep$pdHess, elapsed_seconds = timing[["elapsed"]],
  fixed_parameters = length(tam_fit$obj$par), random_parameters = length(tam_fit$obj$env$random))
native_diagnostics <- data.frame(model = "SAM reference", is_converged = NA,
  convergence_code = sam_fit$opt$convergence, optimizer_message = sam_fit$opt$message,
  objective = as.numeric(sam_fit$opt$objective),
  max_abs_gradient = max(abs(sam_fit$sdrep$gradient.fixed)),
  sdreport_success = all(is.finite(sam_fit$sdrep$sd)),
  positive_definite_hessian = sam_fit$sdrep$pdHess, elapsed_seconds = NA_real_,
  fixed_parameters = length(sam_fit$sdrep$par.fixed), random_parameters = length(sam_fit$sdrep$par.random))
diagnostics <- rbind(native_diagnostics, diagnostics)
write.csv(diagnostics, file.path(dir, "results", "diagnostics.csv"), row.names = FALSE)
tabs <- tidy_tam(model_list = models)
for (nm in names(tabs$pop)) write.csv(tabs$pop[[nm]], file.path(dir, "results", paste0("native_", nm, ".csv")), row.names = FALSE)
for (nm in names(tabs$obs_pred)) write.csv(tabs$obs_pred[[nm]], file.path(dir, "results", paste0("predictions_", nm, ".csv")), row.names = FALSE)
write.csv(tabs$fixed_par, file.path(dir, "results", "parameters.csv"), row.names = FALSE)
vis_tam(model_list = models, output_file = file.path(dir, "results", "SAM_tinyAM_dashboard.html"),
        open_file = FALSE, render_args = list(quiet = TRUE))
source(file.path(dir, "002_summarize.R"))
print(diagnostics)
