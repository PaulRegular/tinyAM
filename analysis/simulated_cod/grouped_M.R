# Reduce M-state dimension; keep the same simulated data, seed and controls.
source("analysis/simulated_cod/run_trials.R")
directory <- "data-raw/cache/sim_cod_exploration/stable-common-40-reference-M03"
result <- readRDS(file.path(directory, "scenario.rds"))
dat <- result$sim_cod_obs
sim_cod_truth <- result$sim_cod_truth
breaks <- c(1, 4, 6, 11, 21)
dat$pop_covariates$M_group <- cut(dat$pop_covariates$age, breaks,
  right = FALSE, include.lowest = TRUE)
dat$pop_covariates$M0_block <- with(dat$pop_covariates,
  exp(ave(log(M0), year, M_group, FUN = mean)))
# Both surveys have at least twelve positive years through age 13; older
# observations share the final q coefficient instead of identifying a tail curve.
dat$index$q_age <- factor(pmin(dat$index$age, 13))
q <- ~ 0 + for_survey("RV", ~ q_age) + for_survey("Sentinel", ~ q_age)
fishing <- ~ 0 + for_fleet("Trawl", ~ logistic(mean_length)) +
  for_fleet("Longline", ~ logistic(mean_length)) +
  for_fleet("Gillnet", ~ double_logistic(mean_length))
F <- list(process = "ar1", mu_form = fishing, sd_form = ~ fleet, mean_ages = 5:10)
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "means"
if (mode == "means") {
for (baseline in c("scalar", "reference")) {
  M <- list(process = "ar1", age_breaks = breaks,
    mu_supplied = if (baseline == "scalar") ~ I(.3) else ~ M0_block, mean_ages = 5:10)
  run_trial(paste0("M-ar1-", baseline, "-four-blocks-free-q"), F, M, q, ages = 1:20)
}
}
if (mode == "F-flexibility") {
  M <- list(process = "ar1", age_breaks = breaks,
    mu_supplied = ~ M0_block, mean_ages = 5:10)
  dat$catch$F_age <- factor(pmin(dat$catch$age, 13))
  F$mu_form <- ~ 0 + fleet:F_age
  run_trial("M-ar1-four-blocks-free-F-q", F, M, q, ages = 1:20)
  F <- list(process = "cor_rw", sd_form = ~ fleet, mean_ages = 5:10)
  run_trial("M-ar1-four-blocks-cor-rw-F-free-q", F, M, q, ages = 1:20)
}
