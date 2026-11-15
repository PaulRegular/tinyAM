# Controlled experiment for the information in sampled zeros.
# Run from the tinyAM root with simulate_pld.R.
pld_simple_ages <- 1:14
pld_simple_years <- 1:30
pld_simple_weak <- c(5:7, 16:18)
pld_simple_sd <- c(catch = 0.25, index = 0.35)
pld_simple_zero_target <- c(catch = 0.10, index = 0.15)

pld_simple_truth <- function(id, pulse = TRUE) {
  ages <- pld_simple_ages
  years <- pld_simple_years
  typical <- withr::with_seed(82000L + id, stats::rlnorm(1L, log(30000), 0.5))
  recruitment <- typical * ifelse(years %in% pld_simple_weak, 0.4, 1)
  fishing <- outer(c(rep(0.06, 7L), rep(0.22, 5L), rep(0.035, 18L)),
    exp(-0.5 * ((ages - 7.7) / 2.5)^2))
  natural <- ifelse(pulse & years >= 14L, 0.60, 0.25)
  mortality <- fishing + natural
  abundance <- matrix(NA_real_, length(years), length(ages))
  abundance[, 1L] <- recruitment
  abundance[1L, -1L] <- recruitment[1L] * exp(-cumsum(mortality[1L, -length(ages)]))
  for (t in years[-1L]) {
    survivors <- abundance[t - 1L, ] * exp(-mortality[t - 1L, ])
    abundance[t, -1L] <- survivors[-length(ages)]
    abundance[t, length(ages)] <- sum(tail(survivors, 2L))
  }
  weight <- 1e-5 * (120 * (1 - exp(-0.1 * (ages + 0.5))))^3
  maturity <- stats::plogis(ages - 5)
  catch <- abundance * fishing * (-expm1(-mortality)) / mortality
  stopifnot(all(is.finite(abundance) & abundance > 0),
    all(catch > 0 & catch < abundance))
  list(id = id, pulse = pulse, mean = list(catch = catch, index = abundance),
    ssb = drop(abundance %*% (weight * maturity)), recruitment = recruitment,
    weight = weight, maturity = maturity, natural = natural)
}

pld_simple_calibration <- function() {
  ids <- 101:106
  histories <- lapply(ids, pld_simple_truth)
  components <- names(pld_simple_zero_target)
  values <- lapply(components, function(component) {
    means <- lapply(histories, function(x) x$mean[[component]])
    target <- pld_simple_zero_target[[component]]
    old_means <- unlist(lapply(means, function(x) x[, -1L]), use.names = FALSE)
    shares <- unlist(lapply(means, function(x) {
      (x / rowSums(x))[, -1L]
    }), use.names = FALSE)
    probability <- function(y50) mean(exp(-log(2) * old_means / y50))
    y50 <- exp(stats::uniroot(function(x) probability(exp(x)) - target,
      range(log(old_means)) + c(-20, 20))$root)
    finite_probability <- function(n) mean(exp(n * log1p(-shares)))
    continuous <- exp(stats::uniroot(function(x) finite_probability(exp(x)) - target,
      c(0, 20))$root)
    choices <- unique(c(floor(continuous), ceiling(continuous)))
    n <- choices[which.min(abs(vapply(choices, finite_probability, numeric(1)) - target))]
    data.frame(component = component, target = target, y50 = y50,
      exact_zero = probability(y50), sample_n = n, finite_zero = finite_probability(n))
  })
  list(ids = ids, values = do.call(rbind, values))
}

pld_simple_observations <- function(truth, calibration, finite = FALSE) {
  seed <- (if (finite) 86000L else 84000L) + truth$id
  values <- calibration$values
  draw <- function(component) {
    mu <- truth$mean[[component]]
    sigma <- pld_simple_sd[[component]]
    setting <- values[values$component == component, ]
    if (finite) {
      n <- setting$sample_n
      expected_total <- rowSums(mu)
      total <- expected_total * stats::rlnorm(nrow(mu), -sigma^2 / 2, sigma)
      counts <- t(vapply(seq_len(nrow(mu)), function(t) {
        drop(stats::rmultinom(1L, n, mu[t, ]))
      }, integer(ncol(mu))))
      estimate <- counts * total / n
      stopifnot(all(rowSums(counts) == n),
        isTRUE(all.equal(rowSums(estimate), total)))
      return(estimate)
    }
    p <- -expm1(-log(2) * mu / setting$y50)
    # Draw all uniforms and normal errors so pulse/control use common streams.
    seen <- matrix(stats::runif(length(mu)), nrow(mu)) < p
    error <- matrix(stats::rnorm(length(mu), -sigma^2 / 2, sigma), nrow(mu))
    (mu / p) * exp(error) * seen
  }
  obs <- withr::with_seed(seed, lapply(stats::setNames(names(truth$mean),
    names(truth$mean)), draw))
  keys <- expand.grid(year = pld_simple_years, age = pld_simple_ages)
  f_phase <- cut(keys$year, c(0, 7, 12, 30), labels = c("early", "peak", "late"))
  m_phase <- factor(ifelse(keys$year >= 14L, "late", "early"), c("early", "late"))
  rec_phase <- factor(ifelse(keys$year %in% pld_simple_weak, "weak", "usual"),
    c("usual", "weak"))
  data <- tinyAM::as_tam_data(list(
    catch = transform(keys, obs = as.vector(obs$catch), unit = "thousand fish", f_phase = f_phase),
    index = transform(keys, obs = as.vector(obs$index), unit = "thousand fish",
      survey = "RV", samp_time = 0),
    weight = transform(keys, obs = truth$weight[age], unit = "kg"),
    maturity = transform(keys, obs = truth$maturity[age], spawn_time = 0),
    pop_covariates = transform(keys, m_phase = m_phase, rec_phase = rec_phase)))
  tinyAM::check_tam_data(data, quiet = TRUE)
  zero <- do.call(rbind, lapply(names(obs), function(component) {
    data.frame(component = component, age = pld_simple_ages,
      observed = nrow(obs[[component]]), zeros = colSums(obs[[component]] == 0))
  }))
  list(data = data, zeros = zero, seed = seed)
}

pld_simple_settings <- function(input, truth, calibration, arm, stage) {
  model <- function(component) {
    if (arm == "lognormal") return(tinyAM::lognormal())
    scale <- calibration$values$y50[match(component, calibration$values$component)]
    tinyAM::pld_lognormal(y50 = if (stage != "finite") scale else NULL)
  }
  precision <- function(component) {
    if (stage == "finite") return(list(sd_form = ~ 1))
    sigma <- unname(pld_simple_sd[[component]])
    list(sd_form = ~ 0, sd_supplied = stats::as.formula(substitute(~ I(s), list(s = sigma))))
  }
  young <- subset(input$data$index, age == 1L & obs > 0)
  initial <- young$obs[match(1L, young$year)]
  if (is.na(initial)) initial <- stats::median(young$obs)
  mean_rec <- vapply(c(FALSE, TRUE), function(weak) {
    stats::median(young$obs[young$year > 1L & (young$year %in% pld_simple_weak) == weak])
  }, numeric(1))
  stopifnot(is.finite(initial), all(is.finite(mean_rec)))
  start <- list(log_r0 = log(initial), rec_beta = setNames(log(mean_rec),
    c("rec_phaseusual", "rec_phaseweak")),
    log_mu_f = setNames(rep(log(0.08), 3L), paste0("f_phase", c("early", "peak", "late"))))
  if (stage == "mortality") {
    m <- list(process = "off", mu_supplied = NULL, mu_form = ~ 0 + m_phase)
    start$mu_m <- setNames(rep(log(0.3), 2L), c("m_phaseearly", "m_phaselate"))
  } else {
    data <- input$data
    data$pop_covariates$known_m <- truth$natural[match(data$pop_covariates$year, pld_simple_years)]
    input$data <- data
    m <- list(process = "off", mu_supplied = ~ known_m)
  }
  list(data = input$data, silent = TRUE, start_par = start,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "off", mu_form = ~ 0 + rec_phase),
    F_settings = list(process = "off", mu_form = ~ 0 + f_phase + normal(age)),
    M_settings = m,
    catch_settings = c(list(family = model("catch")), precision("catch")),
    index_settings = c(list(family = model("index"), q_form = ~ 0), precision("index")))
}

pld_simple_metrics <- function(fit, truth, stage) {
  older <- subset(fit$pop$N, age == 14L & year %in% 16:23)
  weak <- subset(fit$pop$N, age %in% 2:8 & (year - age + 1L) %in% pld_simple_weak)
  estimated_rec <- tinyAM::get_recruitment(fit)
  estimated_ssb <- tinyAM::get_ssb(fit)
  compare <- function(est, value, label, interval = NULL) {
    delta <- log(est / value)
    data.frame(metric = label, absolute = mean(abs(delta)), bias = mean(delta),
      coverage = if (is.null(interval)) NA_real_ else mean(value >= interval$lwr & value <= interval$upr),
      width = if (is.null(interval)) NA_real_ else mean(log(interval$upr / interval$lwr)))
  }
  metrics <- do.call(rbind, list(
    compare(older$est, truth$mean$index[cbind(older$year, older$age)], "Age 14+ abundance"),
    compare(weak$est, truth$mean$index[cbind(weak$year, weak$age)], "Weak cohort abundance"),
    compare(estimated_rec$est[pld_simple_weak], truth$recruitment[pld_simple_weak],
      "Weak recruitment", estimated_rec[pld_simple_weak, ]),
    compare(estimated_ssb$est, truth$ssb, "SSB", estimated_ssb)))
  change <- NULL
  if (stage == "mortality") {
    m <- which(names(fit$opt$par) == "mu_m")
    stopifnot(length(m) == 2L)
    level <- exp(unname(fit$opt$par[m]))
    change_est <- diff(level)
    slope <- c(-level[1L], level[2L])
    se <- sqrt(drop(crossprod(slope, fit$sdrep$cov.fixed[m, m] %*% slope)))
    interval <- change_est + c(-1, 1) * stats::qnorm(0.975) * se
    target <- truth$natural[14L] - truth$natural[13L]
    change <- data.frame(estimate = change_est, truth = target, lower = interval[1L],
      upper = interval[2L], absolute = abs(change_est - target), bias = change_est - target,
      coverage = interval[1L] <= target && interval[2L] >= target,
      positive = interval[1L] > 0)
  }
  list(metrics = metrics, mortality = change,
    ssb_ratio = stats::median(estimated_ssb$est / truth$ssb))
}

pld_simple_minimize <- function(start, objective, gradient, lower, upper, control, exact_hessian) {
  opt <- stats::nlminb(start, objective, gradient, lower = lower, upper = upper,
    control = control)
  opt$study_newton <- FALSE
  if (opt$convergence != 0L) return(opt)
  g <- as.numeric(gradient(opt$par))
  factor <- tryCatch(chol(exact_hessian(opt$par)), error = function(e) NULL)
  if (is.null(factor) || any(!is.finite(g))) return(opt)
  direction <- backsolve(factor, forwardsolve(t(factor), g))
  for (fraction in 2^-(0:8)) {
    candidate <- opt$par - fraction * direction
    if (any(candidate < lower | candidate > upper | !is.finite(candidate))) next
    value <- objective(candidate)
    if (!is.finite(value) || value >= opt$objective) next
    updated_gradient <- gradient(candidate)
    if (any(!is.finite(updated_gradient)) || max(abs(updated_gradient)) >= max(abs(g))) next
    opt$par <- candidate
    opt$objective <- value
    opt$study_newton <- TRUE
    break
  }
  opt
}

pld_simple_optimizer <- function(fit, optimizer) {
  replaced <- 0L
  replace_optimizer <- function(expr) {
    if (!is.call(expr)) return(expr)
    parts <- as.list(expr)
    if (identical(parts[[1L]], quote(stats::nlminb))) {
      replaced <<- replaced + 1L
      if (optimizer == "tight") {
        parts$control <- substitute(utils::modifyList(old, list(rel.tol = 1e-13)),
          list(old = parts$control))
      } else {
        parts[[1L]] <- as.name("pld_simple_minimize")
        parts$exact_hessian <- quote(obj$he)
      }
      return(as.call(parts))
    }
    for (i in seq_along(parts)[-1L]) {
      if (is.call(parts[[i]])) parts[[i]] <- replace_optimizer(parts[[i]])
    }
    as.call(parts)
  }
  body(fit) <- replace_optimizer(body(fit))
  if (optimizer == "newton") {
    env <- new.env(parent = environment(fit))
    env$pld_simple_minimize <- pld_simple_minimize
    environment(fit) <- env
  }
  stopifnot(replaced == 1L)
  fit
}

pld_simple_fit <- function(input, truth, calibration, arm, stage, optimizer = "default") {
  fitter <- if (arm == "positive") pld_occurrence_ablated_fit() else tinyAM::fit_tam
  if (optimizer != "default") fitter <- pld_simple_optimizer(fitter, optimizer)
  spec <- pld_simple_settings(input, truth, calibration, arm, stage)
  fit <- do.call(fitter, spec)
  checks <- tinyAM::check_tam(fit)
  list(converged = isTRUE(checks$is_converged), numerical = checks$numerical,
    gradient = checks$max_gradient, optimizer = fit$opt$convergence,
    study_newton = isTRUE(fit$opt$study_newton), objective = fit$opt$objective,
    values = if (isTRUE(checks$is_converged)) pld_simple_metrics(fit, truth, stage))
}
