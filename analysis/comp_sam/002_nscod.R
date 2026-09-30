# Run from the tinyAM repository root after 001_download.R.
# For development, load the local package with pkgload::load_all() first.
library(tinyAM)
here <- file.path("analysis", "comp_sam")
out <- file.path(here, "results")
dir.create(out, showWarnings = FALSE)
source <- read_sam_files(file.path(here, "source"), conf = file.path(here, "nscod.cfg"))
obs <- sam_to_tam_obs(source)
stopifnot(check_obs(obs))
audit <- sam_tam_assumptions(source)
utils::write.csv(audit, file.path(out, "assumptions.csv"), row.names = FALSE)
for (nm in names(obs)) utils::write.csv(obs[[nm]], file.path(out, paste0("input_", nm, ".csv")), row.names = FALSE)

# Load only public serialized values. Never call fit$obj or a SAM optimizer.
reference_env <- new.env(parent = baseenv())
load(file.path(here, "source", "fit.expected.Rdata"), reference_env)
fit <- reference_env[["fit.exp"]]
reference <- sam_reference(fit,
  "fishfollower/SAM c6cfd035c7de59f7b3421dde31901efbca4cb0e8: stockassessment/tests/nscod/fit.expected.Rdata")

# Verify the saved reference uses the same observations, fleet IDs, years and
# biological inputs. The only shape differences are legacy single-fleet arrays.
saved_obs <- data.frame(year = fit$data$aux[, 1], fleet_id = fit$data$aux[, 2], age = fit$data$aux[, 3], obs = exp(fit$data$logobs))
translated <- rbind(obs$catch[, c("year", "fleet_id", "age", "obs")], obs$index[, c("year", "fleet_id", "age", "obs")])
key <- function(d) paste(d$year, d$fleet_id, d$age, sep = ":")
compare <- translated$obs[match(key(saved_obs), key(translated))]
compare[compare <= 0] <- NA_real_
stopifnot(isTRUE(all.equal(compare, saved_obs$obs, check.attributes = FALSE)))
for (pair in list(c("sw", "stockMeanWeight"), c("mo", "propMat"), c("nm", "natMor"),
                 c("cw", "catchMeanWeight"), c("pf", "propF"), c("pm", "propM"))) {
  stopifnot(isTRUE(all.equal(as.numeric(source$data[[pair[1]]]), as.numeric(fit$data[[pair[2]]]))))
}
stopifnot(identical(as.numeric(source$fleets$samp_time), as.numeric(fit$data$sampleTimes)))
for (nm in intersect(names(source$conf), names(fit$conf))) {
  if (nm == "maxAgePlusGroup") next # legacy scalar versus modern per-fleet flags
  current <- source$conf[[nm]]
  saved <- fit$conf[[nm]]
  if (nm == "fixVarToWeight") saved <- rep(saved, length.out = length(current))
  stopifnot(isTRUE(all.equal(current, saved, check.attributes = FALSE)))
}

# Independently verify the stored reference summaries/predictions. This uses
# SAM's equations, not tinyAM's different summary conventions.
N <- exp(t(fit$pl$logN))
F <- exp(t(fit$pl$logF))
M <- fit$data$natMor
Z <- F + M
stopifnot(isTRUE(all.equal(reference$tables$SSB$est,
                         rowSums(N * fit$data$stockMeanWeight * fit$data$propMat), check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(reference$tables$recruitment$est, N[, 1], check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(reference$tables$Fbar$est, rowMeans(F[, 2:4]), check.attributes = FALSE)))
for (i in seq_len(nrow(saved_obs))) {
  y <- match(saved_obs$year[i], fit$data$years)
  a <- saved_obs$age[i] - fit$conf$minAge + 1L
  f <- saved_obs$fleet_id[i]
  if (fit$data$fleetTypes[f] == 0) {
    prediction <- log(N[y, a] * F[y, a] / Z[y, a] * (-expm1(-Z[y, a])))
    sy <- match(saved_obs$year[i], fit$conf$keyScaledYears)
    if (!is.na(sy)) prediction <- prediction - fit$pl$logScale[fit$conf$keyParScaledYA[sy, a] + 1L]
  } else {
    prediction <- log(N[y, a]) - Z[y, a] * fit$data$sampleTimes[f] + fit$pl$logFpar[fit$conf$keyLogFpar[f, a] + 1L]
  }
  stopifnot(abs(prediction - fit$rep$predObs[i]) < 1e-8)
}
for (nm in names(reference$tables)) utils::write.csv(reference$tables[[nm]], file.path(out, paste0("SAM_", nm, ".csv")), row.names = FALSE)
utils::write.csv(reference$availability, file.path(out, "reference_availability.csv"), row.names = FALSE)

# An explicitly simplified structural model, not an exact SAM replication.
# F increment correlation and age-specific F variance are removed; SAM catch
# scaling is omitted. All raw observations, M, q and observation-SD groups stay.
# No optimization is attempted in this first translation/audit phase.
dat <- make_dat(obs, years = 1963:2015, ages = 1:6,
  N_settings = list(process = "iid", init = "free"),
  F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:4),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
  catch_settings = list(sd_form = ~ 0 + sd_block, fill_missing = FALSE),
  index_settings = list(q_form = ~ 0 + q_block, sd_form = ~ 0 + sd_block, fill_missing = FALSE),
  proj_settings = list(n_proj = 0))
par <- make_par(dat)
obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
stopifnot(is.finite(obj$fn(obj$par)), all(is.finite(obj$gr(obj$par))))
summary <- data.frame(year_min = min(dat$years), year_max = max(dat$years),
  min_age = min(dat$ages), max_age = max(dat$ages),
  catch_rows = nrow(obs$catch), survey_rows = nrow(obs$index),
  q_blocks = ncol(dat$q_modmat), catch_sd_blocks = ncol(dat$sd_catch_modmat),
  index_sd_blocks = ncol(dat$sd_index_modmat), initial_joint_nll = obj$fn(obj$par),
  source_reference_objective = as.numeric(fit$opt$objective))
utils::write.csv(summary, file.path(out, "translation_diagnostics.csv"), row.names = FALSE)
saveRDS(list(source = source, obs = obs, audit = audit, reference = reference,
             dat = dat, par = par), file.path(out, "nscod_bridge.rds"))
cat("Inputs validated; public reference checked; simplified model has finite objective and gradient.\n")
print(audit[audit$tam_status %in% c("unsupported", "partially_supported"), c("component", "sam_setting", "tam_status")])
