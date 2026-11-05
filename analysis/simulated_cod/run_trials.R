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
  if (!is.data.frame(truth) || !nrow(truth) ||
      !all(c("year", "SSB", "N", "recruitment") %in% names(truth)) ||
      !is.numeric(truth$year) ||
      any(!is.finite(truth$year)) || anyDuplicated(truth$year)) {
    stop("Truth must have one finite numeric year per annual row.")
  }
  tables <- list(SSB = fit$pop$ssb, abundance = fit$pop$abundance,
    recruitment = fit$pop$recruitment)
  columns <- c(SSB = "SSB", abundance = "N", recruitment = "recruitment")
  do.call(rbind, lapply(names(tables), function(quantity) {
    estimated <- tables[[quantity]]
    if (!is.data.frame(estimated) ||
        !all(c("year", "est", "is_proj") %in% names(estimated)) ||
        anyNA(estimated$is_proj)) {
      stop("Missing annual comparison columns for ", quantity, ".")
    }
    estimated <- estimated[!estimated$is_proj, , drop = FALSE]
    if (!is.numeric(estimated$year) || any(!is.finite(estimated$year)) ||
        anyDuplicated(estimated$year) || nrow(estimated) != nrow(truth) ||
        !setequal(estimated$year, truth$year)) {
      stop("The ", quantity, " comparison requires exactly one estimate for every truth year.")
    }
    position <- match(estimated$year, truth$year)
    relative <- estimated$est / truth[[columns[[quantity]]]][position] - 1
    if (any(!is.finite(relative))) stop("Incomplete or nonfinite truth comparison.")
    data.frame(quantity = quantity, n_years = length(relative),
      mean_relative_error = mean(relative), mean_absolute_relative_error = mean(abs(relative)),
      terminal_relative_error = relative[which.max(estimated$year)])
  }))
}

run_trial <- function(name, F_settings, M_settings, q_form = catchability, ages = 1:14) {
  path <- file.path(directory, paste0(name, ".rds"))
  if (file.exists(path)) stop("Preserve the existing attempt: ", path)
  warnings <- character()
  error <- NULL
  input_file <- tempfile()
  saveRDS(list(data = dat, truth = sim_cod_truth$annual), input_file)
  input_hash <- unname(tools::md5sum(input_file))
  unlink(input_file)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(fit_tam(dat, ages = ages,
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
    input_hash = input_hash,
    packaged_baseline_hashes = tools::md5sum(c("data/sim_cod_obs.rda", "data/sim_cod_truth.rda")),
    criterion = "Numerical checks pass and absolute mean annual relative SSB error < 0.10; single-realization criterion.",
    accepted = FALSE)
  saveRDS(record, path)
  if (!is.null(fit)) {
    record$diagnostics <- check_tam(fit)
    record$errors <- errors(fit)
    ssb_error <- record$errors$mean_relative_error[record$errors$quantity == "SSB"]
    record$accepted <- isTRUE(record$diagnostics$is_converged) &&
      length(ssb_error) == 1L && is.finite(ssb_error) && abs(ssb_error) < 0.10
  }
  saveRDS(record, path)
  cat("\n", name, "\n")
  print(record[c("elapsed_seconds", "warnings", "error", "accepted")])
  if (!is.null(fit)) {
    print(record$diagnostics)
    print(record$errors)
  }
  invisible(record)
}

if (sys.nframe() == 0L) {
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
if (mode == "q-age") {
  dat$index$q_age <- factor(pmin(dat$index$age, 14L))
  q_age <- ~ 0 + for_survey("RV", ~ q_age) + for_survey("Sentinel", ~ q_age)
  run_trial("M-ar1-separable-F-q-age", list(process = "off", mu_form = separable,
    mean_ages = 5:10), M, q_age)
  M_estimated <- list(process = "ar1", mu_form = ~ 1, mu_supplied = NULL, mean_ages = 5:10)
  run_trial("M-ar1-estimated-mean-separable-F-q-age", list(process = "off", mu_form = separable,
    mean_ages = 5:10), M_estimated, q_age)
}
if (mode == "q-year") {
  dat$index$q_age <- factor(pmin(dat$index$age, 14L))
  dat$index$availability_group <- factor(pmin(dat$index$age, 5L))
  q_year <- ~ 0 +
    for_survey("RV", ~ q_age + ar1(year, by = availability_group)) +
    for_survey("Sentinel", ~ q_age + ar1(year, by = availability_group))
  run_trial("M-ar1-separable-F-q-age-year", list(process = "off", mu_form = separable,
    mean_ages = 5:10), M, q_year)
  M_point_four <- list(process = "ar1", mu_supplied = ~ I(0.4), mean_ages = 5:10)
  run_trial("M-ar1-point-four-separable-F-q-age-year", list(process = "off", mu_form = separable,
    mean_ages = 5:10), M_point_four, q_year)
}
if (mode == "known-M-diagnostic") {
  dat$index$q_age <- factor(pmin(dat$index$age, 14L))
  dat$index$availability_group <- factor(pmin(dat$index$age, 5L))
  q_year <- ~ 0 +
    for_survey("RV", ~ q_age + ar1(year, by = availability_group)) +
    for_survey("Sentinel", ~ q_age + ar1(year, by = availability_group))
  p <- sim_cod_truth$population
  dat$pop_covariates <- data.frame(year = p$year, age = p$age,
    M_known = ifelse(is.finite(p$M), p$M, p$Z))
  known <- list(process = "off", mu_supplied = ~ M_known, mean_ages = 5:10)
  run_trial("known-M-diagnostic-separable-F", list(process = "off", mu_form = separable,
    mean_ages = 5:10), known, q_year, ages = 1:20)
  run_trial("known-M-diagnostic-flexible-F", list(process = "ar1", mu_form = flexible,
    sd_form = ~ fleet, mean_ages = 5:10), known, q_year, ages = 1:20)
}
}
