# Render the saved equivalence comparison without refitting either model.
# Run from the repository root with R 4.4:
# Rscript analysis/sam_equivalence/render_haddock.R [comparison.rds] [dashboard.html]
args <- commandArgs(TRUE)
pkgload::load_all(".", quiet = TRUE)
input <- if (length(args)) args[1L] else "analysis/sam_equivalence/results/haddock/comparison.rds"
output <- if (length(args) > 1L) args[2L] else file.path(dirname(input), "haddock-comparison.html")
x <- readRDS(input)
fit <- x$fit; sam <- x$native
stopifnot(x$SAM_revision == "1cc464b80f6f2a4a38b35cce276f13e521761da3")
years <- as.integer(sam$data$years); ages <- 0:8
N <- t(exp(sam$pl$logN)); F <- t(exp(sam$pl$logF)); M <- sam$data$natMor
stopifnot(identical(dim(N), dim(F)), identical(dim(N), dim(M)))
Z <- F + M
grid <- expand.grid(year = years, age = ages)
z <- qnorm(.975)
table_at_age <- function(values, unit, log_se = NULL) {
  d <- transform(grid, est = as.vector(values), unit = unit, is_proj = FALSE)
  d$se <- d$lwr <- d$upr <- NA_real_; d$se_scale <- NA_character_
  if (!is.null(log_se)) {
    d$se <- as.vector(log_se); d$se_scale <- "log"
    d$lwr <- exp(log(d$est) - z * d$se)
    d$upr <- exp(log(d$est) + z * d$se)
  }
  d
}
annual <- function(values, unit) data.frame(year = years, est = values,
  se = NA_real_, se_scale = NA_character_, lwr = NA_real_, upr = NA_real_,
  unit = unit, is_proj = FALSE)
report <- summary(sam$sdrep, "report")
reported <- function(name, unit) {
  a <- report[rownames(report) == name, , drop = FALSE]
  stopifnot(nrow(a) <= length(years))
  d <- annual(rep(NA_real_, length(years)), unit)
  j <- seq_len(nrow(a))
  d$est[j] <- exp(a[, 1L]); d$se[j] <- a[, 2L]; d$se_scale[j] <- "log"
  d$lwr[j] <- exp(a[, 1L] - z * a[, 2L]); d$upr[j] <- exp(a[, 1L] + z * a[, 2L])
  d
}
W <- sam$data$stockMeanWeight
WC <- matrix(sam$data$catchMeanWeight, nrow = dim(sam$data$catchMeanWeight)[1L])
# No native catch weights in the final, survey-only year: do not invent them.
if (nrow(WC) < length(years)) WC <- rbind(WC, rep(NA_real_, length(ages)))
ssb_age <- N * W * sam$data$propMat
mean_rate <- function(rate, selected) {
  j <- match(selected, ages)
  rowSums(rate[, j, drop = FALSE] * N[, j, drop = FALSE]) / rowSums(N[, j, drop = FALSE])
}
pop <- list(N = table_at_age(N, "thousand fish", t(sam$plsd$logN)),
  F = table_at_age(F, "per year", t(sam$plsd$logF)),
  M = table_at_age(M, "per year"), Z = table_at_age(Z, "per year"),
  W = table_at_age(W, "kg"), W_catch = table_at_age(WC, "kg"),
  biomass_mat = table_at_age(N * W, "tonnes"), ssb_mat = table_at_age(ssb_age, "tonnes"),
  abundance = annual(rowSums(N), "thousand fish"),
  biomass = reported("logtsb", "tonnes"), ssb = reported("logssb", "tonnes"),
  recruitment = reported("logR", "thousand fish"),
  F_bar = annual(mean_rate(F, fit$dat$F_settings$mean_ages), "per year"),
  M_bar = annual(mean_rate(M, fit$dat$M_settings$mean_ages), "per year"),
  total_yield_pred = reported("logCatch", "tonnes"))
pop$ssb$definition <- "Start-of-year mature biomass; native spawning fractions are zero."
pop$F_bar$definition <- paste("Abundance-weighted F over ages", paste(range(fit$dat$F_settings$mean_ages), collapse = "–"))
pop$M_bar$definition <- "Abundance-weighted M over the tinyAM selected mean ages."
stopifnot(max(abs(pop$ssb$est - rowSums(ssb_age))) < 1e-6,
  max(abs(pop$recruitment$est - N[, 1L])) < 1e-5)

# Copy observation keys/values only. Every prediction is the native SAM result.
aux <- as.data.frame(sam$data$aux); names(aux)[1:3] <- c("year", "fleet", "age")
obs_pred <- lapply(c("catch", "index"), function(type) {
  d <- fit$obs_pred[[type]]
  keep <- intersect(c("year", "age", "survey", "samp_time", "q_key", "q_power_group",
    "sd_key", "relative_sd", "unit", "obs", "is_proj"), names(d))
  d <- d[keep]
  fleet <- if (type == "catch") rep(1L, nrow(d)) else match(d$survey, c("IBTS Q1", "IBTS Q3+Q4")) + 1L
  i <- match(paste(d$year, fleet, d$age), paste(aux$year, aux$fleet, aux$age))
  log_pred <- sam$rep$predObs[i]
  sd_key <- sam$conf$keyVarObs[cbind(fleet, d$age + 1L)] + 1L
  d$sd <- exp(sam$pl$logSdLogObs[sd_key])
  weight <- sam$data$weight[i]
  relative_sd <- ifelse(is.finite(weight), 1 / sqrt(weight), 1)
  d$sd <- d$sd * relative_sd
  d$pred <- exp(log_pred)
  d$std_res <- ifelse(is.finite(d$obs) & d$obs > 0, (log(d$obs) - log_pred) / d$sd, NA_real_)
  d$osa_res <- NA_real_
  if (type == "index") {
    coord <- cbind(match(d$year, years), match(d$age, ages))
    A <- N[coord] * exp(-d$samp_time * Z[coord])
    q_key <- sam$conf$keyLogFpar[cbind(fleet, d$age + 1L)] + 1L
    power_key <- sam$conf$keyQpow[cbind(fleet, d$age + 1L)] + 1L
    d$q_base <- exp(sam$pl$logFpar[q_key]); d$q_power <- 1
    selected <- power_key > 0L
    d$q_power[selected] <- exp(sam$pl$logQpow[power_key[selected]])
    d$q <- d$pred / A
    stopifnot(max(abs(log(d$pred) - log(d$q_base) - d$q_power * log(A))) < 1e-8)
    d$q_unit <- d$q_base_unit <- NA_character_
    # The saved fit has baseline-q uncertainty, not effective-q uncertainty.
    d$q_lwr <- d$q_upr <- d$q_se <- NA_real_; d$q_se_scale <- NA_character_
    d$q_base_se <- sam$plsd$logFpar[q_key]; d$q_base_se_scale <- "log"
    d$q_base_lwr <- exp(log(d$q_base) - z * d$q_base_se)
    d$q_base_upr <- exp(log(d$q_base) + z * d$q_base_se)
  }
  d
})
names(obs_pred) <- c("catch", "index")
catch <- obs_pred$catch
sum_catch <- function(column, weight = rep(1, nrow(catch))) vapply(years, function(y) {
  v <- catch[[column]][catch$year == y] * weight[catch$year == y]
  if (!length(v) || anyNA(v)) NA_real_ else sum(v)
}, numeric(1))
pop$total_catch_pred <- annual(sum_catch("pred"), "thousand fish")
pop$total_catch <- annual(sum_catch("obs"), "thousand fish")
catch_w <- WC[cbind(match(catch$year, years), match(catch$age, ages))]
pop$total_yield <- annual(sum_catch("obs", catch_w), "tonnes")

# Fixed-effect labels follow the matched specification; values/SEs come from SAM.
fixed <- fit$fixed_par
key <- function(coef) as.integer(sub(".*K", "", coef)) + 1L
for (i in seq_len(nrow(fixed))) {
  p <- fixed$par[i]; k <- key(fixed$coef[i])
  source <- switch(p, sd_r = "logSdLogN", sd_f = "logSdLogFsta",
    rho_f = "itrans_rho", sd_beta_n = "logSdLogN", sd_catch = "logSdLogObs",
    sd_index = "logSdLogObs", q = "logFpar", q_power = "logQpow")
  stopifnot(!is.null(source))
  if (p %in% c("sd_r", "sd_f", "rho_f")) k <- 1L
  value <- sam$pl[[source]][k]; se <- sam$plsd[[source]][k]
  transform <- if (p == "rho_f") tanh else if (p == "sd_beta_n") identity else exp
  fixed$est[i] <- transform(value); fixed$se[i] <- se
  fixed$lwr[i] <- transform(value - z * se); fixed$upr[i] <- transform(value + z * se)
}
ref <- structure(list(call = quote(SAM(NShaddock_WGNSSK2026_Run1)),
  dat = list(obs = fit$dat$obs, years = years, ages = ages, is_proj = rep(FALSE, length(years))),
  pop = pop, obs_pred = obs_pred, fixed_par = fixed, random_par = list()), class = c("tam_ref", "list"))
stopifnot(inherits(ref, "tam_ref"), !inherits(ref, "tam_fit"))
for (metric in c("N", "F", "ssb", "recruitment", "F_bar")) {
  stopifnot(max(abs(ref$pop[[metric]]$est / fit$pop[[metric]]$est - 1), na.rm = TRUE) < 1e-7)
}
for (type in names(obs_pred)) stopifnot(max(abs(obs_pred[[type]]$pred / fit$obs_pred[[type]]$pred - 1), na.rm = TRUE) < 1e-7)
saveRDS(ref, file.path(dirname(output), "SAM-reference.rds"))
background <- c("### Northern Shelf haddock: tinyAM and SAM",
  "Both traces use the matched 2026 official SAM specification and the same supplied data. SAM is a precomputed tam_ref, not a tinyAM refit. Estimates overlap to numerical precision.",
  "Abundance and recruitment are thousand fish; biomass and yield are tonnes; average weights are kg. SSB is start-of-year mature biomass (native spawning fractions are zero). The 8+ group is unchanged.",
  "For comparison, average F and M are calculated from native SAM states with tinyAM's abundance weighting and selected ages. SAM's own arithmetic Fbar is a different summary. Confidence intervals for these redefined averages and total abundance are unavailable and are left missing.",
  "Native SAM log-scale SEs are displayed as 95% log-Wald intervals (1.96 SD), matching tinyAM's interval probability. State SEs use the saved plsd output including fixed-parameter uncertainty. Derived effective-q intervals are unavailable for SAM; baseline-q intervals are retained in the tables. No missing uncertainty has been borrowed from tinyAM.",
  "SAM has no catch weights for the survey-only year 2026: its catch-weight and yield cells are missing. The tinyAM preparation repeats 2025 catch weights solely to satisfy a complete grid; 2026 yield should not be interpreted as a source comparison.",
  "The reference supplies native predictions and standardized observation residuals in the output tables. The dashboard's residual diagnostics and process diagnostics are calculated for the tinyAM fit; they are not an independent SAM diagnostic calculation. No OSA calculation or refitting is performed.",
  "See analysis/sam_equivalence/README.md for the pinned revision, equivalence tests, scientific advisories and starting-value limitation.")
vis_tam(model_list = list(tinyAM = fit, SAM = ref), output_file = normalizePath(output, mustWork = FALSE),
  open_file = FALSE, render_args = list(quiet = TRUE), background = background)
cat("Dashboard:", normalizePath(output), "\n")
