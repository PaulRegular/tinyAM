# Explicit survey-only revision. Preserve the failed allocation and all caches.
# Reuse complete population/fishery stages only after verifying that their
# sources and every setting except survey density are identical.
source("analysis/simulated_cod/run_trials.R")
source("data-raw/sim_cod_obs.R")
helpers <- new.env(parent = globalenv())
sys.source("data-raw/sim_cod_helpers.R", helpers)
old_cache <- "data-raw/cache/sim_cod/6584a885a4efffc92f50e55db9d88df5"
report <- readRDS(file.path(old_cache, "generation-report.rds"))
base <- report$provenance
base$reference_M <- NULL
stopifnot(identical(helpers$.sim_cod_hash(base), basename(old_cache)),
  identical(unname(tools::md5sum(file.path("data-raw", names(base$sources)))), unname(base$sources)))
settings <- base$settings
settings$survey_density[["Sentinel"]] <- 4 / 1000
stopifnot(identical(settings[names(settings) != "survey_density"],
  base$settings[names(base$settings) != "survey_density"]))
next_base <- base
next_base$settings <- settings
fingerprint <- helpers$.sim_cod_hash(next_base)
cache <- file.path("data-raw/cache/sim_cod", fingerprint)
if (dir.exists(cache)) stop("Preserve the existing cache: ", cache)
dir.create(cache, recursive = TRUE)
reuse <- list(fingerprint = basename(old_cache), source_files = base$sources,
  change = "Only Sentinel set density: 6 to 4 per 1000 km2; no population/fishery input changed.")
for (stage in c("population", "fisheries")) {
  saved <- readRDS(file.path(old_cache, paste0(stage, ".rds")))
  stopifnot(identical(saved$fingerprint, basename(old_cache)), is.null(saved$error))
  saved$fingerprint <- fingerprint
  saved$reused_from <- reuse
  saveRDS(saved, file.path(cache, paste0(stage, ".rds")))
}
directory <- "data-raw/cache/sim_cod_exploration/teaching-common-20-M03-survey4"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(settings, file.path(directory, "prospective-settings.rds"))
result <- build_sim_cod_obs("C:/Users/paulr/Documents/Codex/simsurvey-data-source",
  settings = settings, cache_dir = cache, write_data = FALSE)
for (stage in c("population", "fisheries")) {
  result$sim_cod_truth$provenance$records[[stage]]$reused_from <- reuse
}
saveRDS(result, file.path(directory, "scenario.rds"))
dat <- result$sim_cod_obs
sim_cod_truth <- result$sim_cod_truth
breaks <- c(1, 4, 6, 11, 21)
dat$pop_covariates$M0_block <- with(dat$pop_covariates,
  exp(ave(log(M0), year, cut(age, breaks, right = FALSE), FUN = mean)))
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
