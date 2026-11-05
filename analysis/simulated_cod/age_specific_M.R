# Same completed scenario; explicit age-specific M states replace the default
# shared actual-M block. No new simulation, starts or optimizer controls.
source("analysis/simulated_cod/run_trials.R")
directory <- "data-raw/cache/sim_cod_exploration/stable-common-40-reference-M03"
result <- readRDS(file.path(directory, "scenario.rds"))
dat <- result$sim_cod_obs
sim_cod_truth <- result$sim_cod_truth
fishing <- ~ 0 +
  for_fleet("Trawl", ~ logistic(mean_length)) +
  for_fleet("Longline", ~ logistic(mean_length)) +
  for_fleet("Gillnet", ~ double_logistic(mean_length))
q <- ~ 0 + for_survey("RV", ~ logistic(mean_length)) +
  for_survey("Sentinel", ~ double_logistic(mean_length))
F <- list(process = "ar1", mu_form = fishing, sd_form = ~ fleet, mean_ages = 5:10)
for (baseline in c("scalar", "reference")) {
  M <- list(process = "ar1", age_breaks = 1:21,
    mu_supplied = if (baseline == "scalar") ~ I(.3) else ~ M0, mean_ages = 5:10)
  run_trial(paste0("M-ar1-", baseline, "-baseline-age-specific"), F, M, q, ages = 1:20)
}
