# Keep the family IDs used in saved study checkpoints.
matched_families <- c(lognormal = "Lognormal", standard_ln = "Logit-link delta lognormal",
  standard_gg = "Logit-link delta generalized gamma", poisson_ln = "Poisson-link delta lognormal",
  poisson_gg = "Poisson-link delta generalized gamma", tweedie = "Tweedie")

matched_family <- function(id, y50 = NULL) {
  switch(id, lognormal = tinyAM::lognormal(),
    standard_ln = tinyAM::delta_lognormal(type = "logit-link"),
    standard_gg = tinyAM::delta_gengamma(type = "logit-link"),
    poisson_ln = tinyAM::delta_lognormal(type = "poisson-link", y50 = y50),
    positive_only = tinyAM::delta_lognormal(type = "poisson-link", y50 = y50),
    poisson_gg = tinyAM::delta_gengamma(type = "poisson-link"), tweedie = tinyAM::tweedie(),
    stop("Unknown family: ", id))
}

# Use the same complete terminal observation for every family. Keep native-age
# biology, including older weights and maturity, for the package's plus group.
matched_pool <- function(data, terminal) {
  for (component in c("catch", "index")) {
    d <- data[[component]]
    if (max(d$age) <= terminal) next
    tail <- d$age >= terminal
    values <- tapply(d$obs[tail], d$year[tail], sum)
    d$obs[d$age == terminal] <- unname(values[as.character(d$year[d$age == terminal])])
    d$obs[d$age > terminal] <- NA_real_
    data[[component]] <- d
  }
  data
}

matched_reference <- function(reference, terminal) {
  if (max(reference$age) <= terminal) return(reference)
  older <- split(reference[reference$age >= terminal, ], reference$year[reference$age >= terminal])
  pooled <- do.call(rbind, lapply(older, function(d) {
    n <- sum(d$N)
    z <- -log(sum(d$N * exp(-d$Z)) / n)
    f <- sum(d$catch_expected) / (n * -expm1(-z) / z)
    data.frame(year = d$year[1], age = terminal, N = n, Z = z, F_realized = f,
      M_realized = z - f, ssb = sum(d$ssb), catch_expected = sum(d$catch_expected),
      index_expected = n)
  }))
  rbind(reference[reference$age < terminal, names(pooled)], pooled)
}

matched_settings <- function(sample, reference, terminal, family, mortality) {
  data <- matched_pool(sample$data, terminal)
  m_center <- mean(reference$M_realized)
  stopifnot(is.finite(m_center), m_center > 0)
  # Validation checks the native grid before prepare_tam() selects modeled
  # ages. Keep native covariate rows and use the pooled rate at the terminal age.
  population <- sample$reference[c("year", "age", "M_realized")]
  modeled <- match(paste(population$year, population$age), paste(reference$year, reference$age))
  keep <- !is.na(modeled)
  population$M_realized[keep] <- reference$M_realized[modeled[keep]]
  data$pop_covariates <- transform(population, M_center = m_center)
  obs <- function(component) {
    if (mortality != "ablation") return(list(family = matched_family(family), sd_form = ~ 1))
    settings <- sample$exact_settings[sample$exact_settings$component == component, ]
    stopifnot(nrow(settings) == 1L, settings$y50 > 0, settings$sd > 0)
    list(family = matched_family(family, settings$y50), sd_form = ~ 0,
      sd_supplied = stats::as.formula(substitute(~ I(s), list(s = settings$sd)), env = baseenv()))
  }
  m <- switch(mortality,
    known = list(process = "off", mu_supplied = ~ M_center),
    estimated = list(process = "ar1", mu_supplied = ~ M_center,
      age_breaks = seq_len(terminal + 1L), first_dev_year = 1L),
    list(process = "off", mu_supplied = ~ M_realized))
  list(data = data, ages = seq_len(terminal), years = 1:30, silent = TRUE,
    N_settings = list(process = "off", init = "exp"),
    R_settings = list(process = "rw", sd = .50),
    F_settings = list(process = "ar1", mu_form = ~ 1, mean_ages = 5:10),
    M_settings = m, start_par = list(log_mu_f = setNames(log(.1), "(Intercept)")),
    catch_settings = obs("catch"), index_settings = c(obs("index"), list(q_form = ~ 0, q_link = "log")))
}

matched_scores <- function(fit, holdout, terminal) {
  holdout <- matched_pool(holdout, terminal)
  rows <- lapply(c("catch", "index"), function(part) {
    d <- fit$obs_pred[[part]]
    d$holdout <- holdout[[part]]$obs[match(paste(d$year, d$age),
      paste(holdout[[part]]$year, holdout[[part]]$age))]
    if (is.null(d$p_zero)) d$p_zero <- 0
    if (is.null(d$pred_positive_median)) d$pred_positive_median <- d$pred
    if (is.null(d$is_structural)) d$is_structural <- d$pred == 0
    rows <- which(is.finite(d$holdout) & !d$is_structural)
    d <- d[rows, ]
    d$component <- part
    pos <- which(d$holdout > 0)
    d$error <- d$covered <- rep(NA_real_, nrow(d))
    d$error[pos] <- log(d$pred_positive_median[pos] / d$holdout[pos])
    if (length(pos)) {
      log_y <- log(d$holdout[pos])
      family <- if (is.null(d$family)) "lognormal" else d$family[pos[1L]]
      cdf <- switch(family, tweedie = {
        log_tail <- tinyAM:::.tweedie_log_cdf(log_y, log(d$pred[pos]), d$sd[pos],
          d$tweedie_power[pos], lower.tail = FALSE)
        -expm1(pmin(0, log_tail - d$log_p_positive[pos]))
      }, delta_gengamma = exp(tinyAM:::.positive_log_cdf(family, log_y,
        d$log_obs_positive_mean[pos], d$sd[pos], d$gengamma_Q[pos])),
        stats::pnorm(log_y, log(d$pred_positive_median[pos]), d$sd[pos]))
      d$covered[pos] <- cdf >= .05 & cdf <= .95
    }
    d[c("component", "year", "age", "holdout", "p_zero", "error", "covered")]
  })
  do.call(rbind, rows)
}

matched_recovery <- function(fit, reference, terminal, mortality) {
  join <- function(d) merge(d, reference, by = c("year", "age"))
  number <- join(fit$pop$N)
  old <- stats::aggregate(cbind(est, N) ~ year, subset(number, age >= 14), sum)
  oldest <- subset(number, age == 20)
  ssb <- tinyAM::get_ssb(fit)
  rec <- tinyAM::get_recruitment(fit)
  true_ssb <- tapply(reference$ssb, reference$year, sum)
  true_rec <- subset(reference, age == 1)
  ssb$truth <- unname(true_ssb[as.character(ssb$year)])
  rec$truth <- true_rec$N[match(rec$year, true_rec$year)]
  low_stock <- ssb$year[order(ssb$truth, ssb$year)[seq_len(10L)]]
  rates <- list(F = join(tinyAM::get_f(fit)), M = join(tinyAM::get_m(fit)), Z = join(fit$pop$Z))
  annual <- lapply(c("F", "M", "Z"), function(name) {
    d <- subset(rates[[name]], age %in% 5:10)
    truth <- switch(name, F = d$F_realized, M = d$M_realized, Z = d$Z)
    d$truth <- truth
    do.call(rbind, lapply(split(d, d$year), function(y) data.frame(year = y$year[1L],
      est = weighted.mean(y$est, y$N), truth = weighted.mean(y$truth, y$N))))
  })
  names(annual) <- c("F", "M", "Z")
  error <- function(label, d, log_scale = TRUE) {
    deviation <- if (log_scale) log(d$est / d$truth) else d$est - d$truth
    covered <- if (all(c("lwr", "upr") %in% names(d))) d$truth >= d$lwr & d$truth <= d$upr else NA
    data.frame(metric = label, absolute = mean(abs(deviation)), bias = mean(deviation),
      coverage = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
      estimate = mean(d$est), truth = mean(d$truth))
  }
  old$truth <- old$N
  oldest$truth <- oldest$N
  if (terminal == 14L) {
    terminal_rows <- subset(number, age == 14L)
    old$lwr <- terminal_rows$lwr[match(old$year, terminal_rows$year)]
    old$upr <- terminal_rows$upr[match(old$year, terminal_rows$year)]
  }
  do.call(rbind, lapply(c("full", "low_stock"), function(window) {
    keep <- function(d) if (window == "low_stock") d[d$year %in% low_stock, ] else d
    mortality_score <- function(name, label) {
      d <- keep(annual[[name]])
      rate <- error(label, d, log_scale = FALSE)
      states <- keep(subset(rates[[name]], age %in% 5:10))
      truth <- switch(name, F = states$F_realized, M = states$M_realized, Z = states$Z)
      if (name != "M" || mortality == "estimated") {
        covered <- truth >= states$lwr & truth <= states$upr
        if (any(!is.na(covered))) rate$coverage <- mean(covered, na.rm = TRUE)
      }
      rate
    }
    estimates <- rbind(error("Age 14+ abundance", keep(old)),
      if (terminal == 20L) error("Age 20+ abundance", keep(oldest)),
      error("SSB", keep(ssb)), error("Recruitment", keep(rec)),
      mortality_score("F", "Fishing mortality"),
      mortality_score("M", "Natural mortality"),
      mortality_score("Z", "Total mortality"))
    estimates$window <- window
    estimates
  }))
}

matched_fit <- function(sample, terminal, family, experiment, mortality) {
  reference <- matched_reference(sample$reference, terminal)
  stopifnot(all(is.finite(reference$M_realized)))
  spec <- matched_settings(sample, reference, terminal, family, mortality)
  fitter <- if (family == "positive_only") delta_occurrence_ablated_fit() else tinyAM::fit_tam
  private <- new.env(parent = environment(fitter))
  # Only distribution summaries are used in the study. Avoid sampling latent
  # states for residual diagnostics that cannot affect fitting or scoring.
  private$tidy_obs_pred <- function(fit, interval, ...) tinyAM:::.tidy_obs_predictions(fit, interval)
  environment(fitter) <- private
  fit <- do.call(fitter, spec)
  check <- tinyAM::check_tam(fit)
  index <- fit$obs_pred$index
  n <- fit$pop$N
  expected <- n$est[match(paste(index$year, index$age), paste(n$year, n$age))]
  stopifnot(all(index$q == 1), isTRUE(all.equal(index$pred, expected, tolerance = 1e-12)))
  if (mortality != "estimated") {
    m <- tinyAM::get_m(fit)
    expected_m <- if (mortality == "known") rep(mean(reference$M_realized), nrow(m)) else
      reference$M_realized[match(paste(m$year, m$age), paste(reference$year, reference$age))]
    stopifnot(isTRUE(all.equal(m$est, expected_m, tolerance = 1e-12)))
  }
  list(converged = isTRUE(check$is_converged), checks = check$numerical,
    gradient = check$max_gradient, objective = fit$opt$objective,
    recovery = if (check$is_converged) matched_recovery(fit, reference, terminal, mortality),
    scores = if (check$is_converged) matched_scores(fit, sample$holdout, terminal),
    parameters = subset(fit$fixed_par, grepl("y50|tweedie_power|gengamma_Q", par)))
}
