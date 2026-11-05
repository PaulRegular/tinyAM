# A deliberately simpler teaching scenario, declared before generation/fitting.
# Same seed and recruitment/growth centres; smaller Z innovations and length
# spread, static common availability and a 20 x 20 grid. No seed search/retries.
source("analysis/simulated_cod/run_trials.R")
source("data-raw/sim_cod_obs.R")
settings <- sim_cod_settings(c(14, 14))
settings$spatial$group_ages <- settings$ages
settings$spatial$group_years <- settings$years
settings$spatial$phi_age <- settings$spatial$phi_year <- 0
settings$mortality$log_sd <- .1
settings$growth$log_sd <- .03
settings$reference_M <- .3
directory <- "data-raw/cache/sim_cod_exploration/teaching-common-20-M03"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(settings, file.path(directory, "prospective-settings.rds"))
cat("Simplified scenario; PID", Sys.getpid(), "\n")
result <- build_sim_cod_obs("C:/Users/paulr/Documents/Codex/simsurvey-data-source",
  settings = settings, write_data = FALSE)
saveRDS(result, file.path(directory, "scenario.rds"))
dat <- result$sim_cod_obs
sim_cod_truth <- result$sim_cod_truth
breaks <- c(1, 4, 6, 11, 21)
dat$pop_covariates$M_group <- cut(dat$pop_covariates$age, breaks,
  right = FALSE, include.lowest = TRUE)
dat$pop_covariates$M0_block <- with(dat$pop_covariates,
  exp(ave(log(M0), year, M_group, FUN = mean)))
fishing <- ~ 0 + for_fleet("Trawl", ~ logistic(mean_length)) +
  for_fleet("Longline", ~ logistic(mean_length)) +
  for_fleet("Gillnet", ~ double_logistic(mean_length))
q <- ~ 0 + for_survey("RV", ~ logistic(mean_length)) +
  for_survey("Sentinel", ~ double_logistic(mean_length))
F <- list(process = "ar1", mu_form = fishing, sd_form = ~ fleet, mean_ages = 5:10)
for (process in c("off", "ar1")) {
  M <- list(process = process, age_breaks = breaks,
    mu_supplied = ~ M0_block, mean_ages = 5:10)
  run_trial(paste0("M-", process, "-four-blocks-curved-q"), F, M, q, ages = 1:20)
}
