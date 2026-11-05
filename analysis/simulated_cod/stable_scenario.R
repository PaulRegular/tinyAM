# Prospective stable-availability benchmark, seed held fixed.
# No fitted values or SSB errors enter fleet-share calibration.
# Original 80 x 80 data and caches remain in the preceding worktree/history.
source("analysis/simulated_cod/run_trials.R")
source("data-raw/sim_cod_obs.R")

settings <- sim_cod_settings(c(7, 7))
settings$spatial$group_ages <- settings$ages
settings$spatial$group_years <- settings$years
# No age/year transitions occur in this static common field; do not inflate its
# initial marginal SD using unused persistence parameters.
settings$spatial$phi_age <- 0
settings$spatial$phi_year <- 0
settings$reference_M <- 0.3
directory <- "data-raw/cache/sim_cod_exploration/stable-common-40-reference-M03"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(settings, file.path(directory, "prospective-settings.rds"))
cat("Generating one predeclared scenario; PID", Sys.getpid(), "\n")
result <- build_sim_cod_obs("C:/Users/paulr/Documents/Codex/simsurvey-data-source",
  settings = settings, write_data = FALSE)
saveRDS(result, file.path(directory, "scenario.rds"))
dat <- result$sim_cod_obs
sim_cod_truth <- result$sim_cod_truth
cat("\nREFERENCE CALIBRATION\n")
print(sim_cod_truth$provenance$reference_M)
print(sim_cod_truth$checks)
fishing <- ~ 0 +
  for_fleet("Trawl", ~ logistic(mean_length)) +
  for_fleet("Longline", ~ logistic(mean_length)) +
  for_fleet("Gillnet", ~ double_logistic(mean_length))
q <- ~ 0 + for_survey("RV", ~ logistic(mean_length)) +
  for_survey("Sentinel", ~ double_logistic(mean_length))
F <- list(process = "ar1", mu_form = fishing, sd_form = ~ fleet, mean_ages = 5:10)
run_trial("M-ar1-scalar-baseline", F,
  list(process = "ar1", mu_supplied = ~ I(0.3), mean_ages = 5:10), q, ages = 1:20)
run_trial("M-ar1-reference-baseline", F,
  list(process = "ar1", mu_supplied = ~ M0, mean_ages = 5:10), q, ages = 1:20)
