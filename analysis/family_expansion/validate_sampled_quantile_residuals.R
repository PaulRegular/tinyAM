# Run from the package root:
# Rscript --vanilla analysis/family_expansion/validate_sampled_quantile_residuals.R [histories] [output]
pkgload::load_all(".", helpers = FALSE, quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
n_histories <- if (length(args)) as.integer(args[1L]) else 40L
if (is.na(n_histories) || n_histories < 1) stop("histories must be a positive integer")
out_dir <- if (length(args) > 1L) args[2L] else
  "analysis/family_expansion/results/sampled-quantile-residuals"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

years <- 1981:2012
ages <- 1:4
recruitment <- 100
natural_mortality <- 0.25
fishing_mortality <- 0.24
process_sd <- 0.45
catch_sd <- 0.30
index_sd <- 0.25
y50 <- 45
scenarios <- c("lognormal", "pld_lognormal")
methods <- c("mle-eb", "mle-mvn")

simulate_history <- function(history, family) {
  set.seed(43800L + history)
  mortality <- natural_mortality + fishing_mortality
  abundance <- matrix(recruitment, length(years), length(ages))
  abundance[1L, ] <- recruitment * exp(-mortality * (ages - 1))
  for (i in seq_along(years)[-1L]) {
    source <- abundance[i - 1L, -length(ages)]
    source[length(source)] <- source[length(source)] + abundance[i - 1L, length(ages)]
    abundance[i, -1L] <- source * exp(-mortality + rnorm(length(source), sd = process_sd))
  }
  catch <- as.numeric(abundance * fishing_mortality / mortality * -expm1(-mortality))
  catch <- catch * exp(rnorm(length(catch), sd = catch_sd))
  index <- as.numeric(abundance)
  index_noise <- rnorm(length(index), sd = index_sd)
  if (family == "pld_lognormal") {
    p_positive <- -expm1(-log(2) * index / y50)
    index <- ifelse(runif(length(index)) < p_positive,
      index / p_positive * exp(index_noise - index_sd^2 / 2), 0)
  } else {
    index <- index * exp(index_noise)
  }
  cells <- expand.grid(year = years, age = ages)
  list(catch = transform(cells, obs = catch, unit = "fish"),
    index = transform(cells, obs = index, survey = "RV", samp_time = 0, unit = "fish"),
    weight = transform(cells, obs = 1, unit = "kg"),
    maturity = transform(cells, obs = pmin(1, age / max(ages)), spawn_time = 0))
}

fit_history <- function(data, family) {
  tinyAM::fit_tam(data,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "off", mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~1),
    M_settings = list(process = "off", mu_supplied = ~ I(0.25)),
    catch_settings = list(sd_form = ~0, sd_supplied = ~ I(0.30)),
    index_settings = list(q_form = ~0, sd_form = ~0, sd_supplied = ~ I(0.25),
      family = if (family == "pld_lognormal") tinyAM::pld_lognormal(y50) else tinyAM::lognormal()),
    start_par = list(log_r0 = log(recruitment), rec_beta = c("(Intercept)" = log(recruitment)),
      log_mu_f = c("(Intercept)" = log(fishing_mortality)), log_sd_n = log(process_sd)),
    silent = TRUE)
}

summarize_residuals <- function(r, history, family, method) {
  do.call(rbind, lapply(split(r, r$type), function(x) {
    z <- x$residual[is.finite(x$residual)]
    data.frame(history = history, family = family, method = method,
      type = x$type[1L], n = length(z), zeros = sum(x$obs == 0, na.rm = TRUE),
      mean = mean(z), sd = sd(z), rms = sqrt(mean(z^2)),
      tail = mean(abs(z) > qnorm(0.975)))
  }))
}

one_history <- function(history, family) {
  warnings <- character()
  outcome <- withCallingHandlers(tryCatch({
    data <- simulate_history(history, family)
    fit <- fit_history(data, family)
    if (!isTRUE(fit$is_converged)) {
      checks <- tinyAM::check_tam(fit)$numerical
      stop(paste(checks$detail[checks$status != "pass"], collapse = " | "))
    }
    result <- lapply(methods, function(method) {
      r <- stats::residuals(fit, type = "quantile", method = method,
        residual_seed = 94000L + history)
      stopifnot(sum(is.finite(r$residual)) == sum(fit$dat$is_observed))
      summarize_residuals(r, history, family, method)
    })
    list(scores = do.call(rbind, result), error = "")
  }, error = function(e) list(scores = NULL, error = conditionMessage(e))),
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  outcome$status <- data.frame(history = history, family = family,
    success = !nzchar(outcome$error), error = outcome$error,
    warnings = paste(unique(warnings), collapse = " | "))
  outcome
}

runs <- lapply(scenarios, function(family) lapply(seq_len(n_histories), one_history, family = family))
runs <- do.call(c, runs)
status <- do.call(rbind, lapply(runs, "[[", "status"))
scores <- do.call(rbind, lapply(runs, "[[", "scores"))
write.csv(status, file.path(out_dir, "attempts.csv"), row.names = FALSE)
if (is.null(scores)) stop("No pair completed. See attempts.csv.")
write.csv(scores, file.path(out_dir, "history-scores.csv"), row.names = FALSE)
estimates <- aggregate(cbind(mean, sd, rms, tail) ~ family + method + type, scores, mean)
first <- subset(scores, method == methods[1L])
second <- subset(scores, method == methods[2L])
paired <- merge(first, second, by = c("history", "family", "type"), suffixes = c("_eb", "_mvn"))
paired$sd_change <- paired$sd_mvn - paired$sd_eb
paired$tail_change <- paired$tail_mvn - paired$tail_eb
uncertainty <- do.call(rbind, lapply(split(paired, list(paired$family, paired$type), drop = TRUE),
  function(d) data.frame(family = d$family[1L], type = d$type[1L], pairs = nrow(d),
    sd_change = mean(d$sd_change), sd_change_se = sd(d$sd_change) / sqrt(nrow(d)),
    tail_change = mean(d$tail_change), tail_change_se = sd(d$tail_change) / sqrt(nrow(d)))))
write.csv(estimates, file.path(out_dir, "calibration.csv"), row.names = FALSE)
write.csv(uncertainty, file.path(out_dir, "paired-changes.csv"), row.names = FALSE)
print(with(status, table(family, success)))
print(estimates, row.names = FALSE, digits = 3)
print(uncertainty, row.names = FALSE, digits = 3)
