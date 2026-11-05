# Exploratory model formulations for the fixed simulated-cod example.
# Run sequentially; original simulations and fits are never overwritten.
# Rscript analysis/simulated_cod/run_trials.R

pkgload::load_all(".", quiet = TRUE, helpers = FALSE)
directory <- "data-raw/cache/sim_cod_exploration/original-data"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
dat <- sim_cod_obs
dat$index$q_age_group <- cut(dat$index$age, c(0, 3, 5, 7, 9, Inf))
catchability <- ~ 0 + for_survey("RV", ~ q_age_group) +
  for_survey("Sentinel", ~ q_age_group)
separable <- ~ 0 +
  for_fleet("Trawl", ~ factor(year) + logistic(age)) +
  for_fleet("Longline", ~ factor(year) + logistic(age)) +
  for_fleet("Gillnet", ~ factor(year) + double_logistic(age))
flexible <- ~ 0 + for_fleet("Trawl", ~ logistic(age)) +
  for_fleet("Longline", ~ logistic(age)) +
  for_fleet("Gillnet", ~ double_logistic(age))

errors <- function(fit, truth = sim_cod_truth$annual) {
  tables <- list(SSB = fit$pop$ssb, abundance = fit$pop$abundance,
    recruitment = fit$pop$recruitment)
  columns <- c(SSB = "SSB", abundance = "N", recruitment = "recruitment")
  do.call(rbind, lapply(names(tables), function(quantity) {
    estimated <- tables[[quantity]]
    estimated <- estimated[!estimated$is_proj, ]
    position <- match(estimated$year, truth$year)
    relative <- estimated$est / truth[[columns[[quantity]]]][position] - 1
    if (any(!is.finite(relative))) stop("Incomplete or nonfinite truth comparison.")
    data.frame(quantity = quantity, n_years = length(relative),
      mean_relative_error = mean(relative), mean_absolute_relative_error = mean(abs(relative)),
      terminal_relative_error = relative[which.max(estimated$year)])
  }))
}

run_trial <- function(name, F_settings, M_settings, q_form = catchability) {
  path <- file.path(directory, paste0(name, ".rds"))
  if (file.exists(path)) stop("Preserve the existing attempt: ", path)
  warnings <- character()
  error <- NULL
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(fit_tam(dat, ages = 1:14,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "rw"),
    F_settings = F_settings, M_settings = M_settings,
    catch_settings = list(sd_form = ~ fleet),
    index_settings = list(q_form = q_form, sd_form = ~ survey), silent = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
  record <- list(name = name, fit = fit, warnings = warnings, error = error,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    source = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    input_hashes = tools::md5sum(c("data/sim_cod_obs.rda", "data/sim_cod_truth.rda")),
    criterion = "Absolute mean annual relative SSB error < 0.10; single-realization criterion.")
  saveRDS(record, path)
  if (!is.null(fit)) {
    record$diagnostics <- check_tam(fit)
    record$errors <- errors(fit)
  }
  saveRDS(record, path)
  cat("\n", name, "\n")
  print(record[c("elapsed_seconds", "warnings", "error")])
  if (!is.null(fit)) {
    print(record$diagnostics)
    print(record$errors)
  }
  invisible(record)
}

mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "M-fixed"
M <- list(process = "ar1", mu_supplied = ~ I(0.3), mean_ages = 5:10)
if (mode == "M-fixed") {
  run_trial("M-ar1-separable-F", list(process = "off", mu_form = separable, mean_ages = 5:10), M)
  run_trial("M-ar1-flexible-F", list(process = "ar1", mu_form = flexible,
    sd_form = ~ fleet, mean_ages = 5:10), M)
}
if (mode == "q-curves") {
  curves <- ~ 0 + for_survey("RV", ~ logistic(age)) +
    for_survey("Sentinel", ~ double_logistic(age))
  run_trial("M-ar1-separable-F-q-curves", list(process = "off", mu_form = separable,
    mean_ages = 5:10), M, curves)
  run_trial("M-ar1-flexible-F-q-curves", list(process = "ar1", mu_form = flexible,
    sd_form = ~ fleet, mean_ages = 5:10), M, curves)
}
