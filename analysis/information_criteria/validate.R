# Run from the repository root. Results are local; every attempt is retained.
# Rscript analysis/information_criteria/validate.R [--quick]

capture_attempt <- function(scenario, seed, run) {
  warnings <- character()
  record <- list(scenario = scenario, seed = seed, numerical_success = FALSE,
    optimizer_code = NA_integer_, optimizer_message = "", pdHess = NA,
    objective = NA_real_, gradient = NA_real_, EDF = NA_real_, penalty = NA_real_,
    optimism = NA_real_, trace_error = NA_real_, error = "")
  tryCatch(withCallingHandlers({
    values <- run(seed)
    record[names(values)] <- values
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) record$error <<- conditionMessage(e))
  record$warnings <- paste(warnings, collapse = " | ")
  as.data.frame(record, stringsAsFactors = FALSE)
}

gaussian_trial <- function(seed) {
  set.seed(seed)
  groups <- 12L
  z <- diag(groups)[rep(seq_len(groups), each = 4L), ]
  truth <- .5 + drop(z %*% rnorm(groups, 0, .45))
  y <- rnorm(length(truth), truth, .3)
  joint <- RTMB::MakeADFun(function(p) {
    -sum(RTMB::dnorm(y, p$beta + z %*% p$b, exp(p$log_sigma), log = TRUE)) -
      sum(RTMB::dnorm(p$b, 0, exp(p$log_tau), log = TRUE))
  }, list(beta = .5, log_sigma = log(.3), log_tau = log(.45), b = rep(0, groups)),
  random = "b", silent = TRUE)
  opt <- nlminb(joint$par, joint$fn, joint$gr)
  joint$fn(opt$par)
  p <- joint$env$parList()
  h <- joint$env$spHess(joint$env$last.par, random = TRUE)
  prior <- Matrix::Diagonal(groups, exp(-2 * p$log_tau))
  trace <- tinyAM:::.aic_trace(h, prior)
  dense <- diag(solve(as.matrix(h), as.matrix(prior)))
  edf <- 2 + groups - sum(trace) # beta and observation SD; tau is prior-only.
  prediction <- p$beta + drop(z %*% p$b)
  train <- -sum(dnorm(y, prediction, exp(p$log_sigma), log = TRUE))
  test <- sum(p$log_sigma + .5 * log(2 * pi) +
    ((truth - prediction)^2 + .3^2) / (2 * exp(2 * p$log_sigma)))
  gradient <- max(abs(joint$gr(opt$par)))
  list(numerical_success = opt$convergence == 0L && gradient < .01,
    optimizer_code = opt$convergence, optimizer_message = opt$message,
    objective = opt$objective, gradient = gradient, EDF = edf,
    penalty = 2 * edf, optimism = 2 * (test - train),
    trace_error = max(abs(trace - dense)))
}

assessment_truth <- function(scenario) {
  grid <- expand.grid(year = 2000:2024, age = 1:6)
  data <- list(catch = transform(grid, obs = 100, unit = "thousand fish"),
    index = transform(grid, obs = 100, survey = "Survey", samp_time = .4,
      q_block = factor(age), unit = "thousand fish"),
    weight = transform(grid, obs = .4 * age, unit = "kg"),
    maturity = transform(grid, obs = plogis(age - 3)))
  settings <- list(N_settings = list(process = if (scenario == "cohort") "iid" else "off"),
    rec_settings = list(process = "iid", mu_form = if (scenario == "BH") ~bh(ssb) else ~1),
    F_settings = list(process = "iid", mu_form = ~1), index_settings = list(q_form = ~1))
  dat <- do.call(prepare_tam, c(list(data = data), settings))
  p <- make_par(dat)
  p$log_r0 <- log(5000)
  p$log_r[] <- log(5000)
  if (!is.null(p$rec_beta)) p$rec_beta[] <- log(5000)
  p$log_mu_f[] <- log(.25)
  p$log_f[] <- log(.25)
  p$log_sd_f <- log(.1)
  p$log_sd_r <- log(.2)
  if (!is.null(p$log_sd_n)) p$log_sd_n <- log(.1)
  p$log_q[] <- log(.6)
  p$log_sd_catch[] <- p$log_sd_index[] <- log(.2)
  if (scenario == "BH") p <- tinyAM:::.initialize_rec_curve(p, dat, initialize_states = TRUE)
  list(dat = dat, par = p, settings = settings)
}

assessment_trial <- function(seed, truth) {
  set.seed(seed)
  d <- truth$dat
  p <- truth$par
  simulation <- nll_fun(p, d, simulate = TRUE)
  states <- intersect(names(p), names(simulation))
  realized <- p
  realized[states] <- simulation[states]
  population <- RTMB::MakeADFun(function(p) nll_fun(p, d), realized, silent = TRUE)$report()
  data <- d$obs
  data$catch$obs <- exp(simulation$log_obs[d$obs_map$type == "catch"])
  data$index$obs <- exp(simulation$log_obs[d$obs_map$type == "index"])
  fit <- do.call(fit_tam, c(list(data = data, start_par = p, silent = TRUE), truth$settings))
  numerical <- list(numerical_success = fit$is_converged, objective = fit$opt$objective,
    optimizer_code = fit$opt$convergence, optimizer_message = fit$opt$message,
    pdHess = isTRUE(fit$sdrep$pdHess), gradient = max(abs(fit$gradient)))
  tryCatch({
    criterion <- cAIC(fit, details = TRUE)
    prediction <- fit$rep$log_pred
    sd <- fit$rep$sd_obs
    train <- -sum(dnorm(simulation$log_obs, prediction, sd, log = TRUE))
    test <- sum(log(sd) + .5 * log(2 * pi) +
      ((population$log_pred - prediction)^2 + population$sd_obs^2) / (2 * sd^2))
    c(numerical, list(EDF = criterion$EDF, penalty = 2 * criterion$EDF,
      optimism = 2 * (test - train)))
  }, error = function(e) c(numerical, list(error = conditionMessage(e))))
}

validation_summary <- function(attempts) {
  summarize <- function(x, subset) {
    ok <- is.finite(x$penalty) & is.finite(x$optimism)
    difference <- x$penalty[ok] - x$optimism[ok]
    data.frame(scenario = x$scenario[1], subset = subset, attempts = nrow(x),
      numerical_success = sum(x$numerical_success), calculated = sum(ok),
      mean_penalty = mean(x$penalty[ok]), mean_optimism = mean(x$optimism[ok]),
      mean_difference = mean(difference), mc_se = sd(difference) / sqrt(length(difference)),
      median_penalty = median(x$penalty[ok]), median_optimism = median(x$optimism[ok]),
      warnings = sum(!is.na(x$warnings) & nzchar(x$warnings)),
      errors = sum(!is.na(x$error) & nzchar(x$error)))
  }
  do.call(rbind, lapply(split(attempts, attempts$scenario), function(x) {
    result <- rbind(summarize(x, "all attempts"),
      summarize(x[x$numerical_success, , drop = FALSE], "numerically converged"))
    result$scenario <- x$scenario[1]
    result
  }))
}

run_validation <- function(quick = FALSE, output = "analysis/information_criteria/results") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(seq_len(if (quick) 5L else 100L), function(i) {
    capture_attempt("Gaussian", 10000L + i, gaussian_trial)
  })
  for (scenario in c("BH", "cohort")) {
    truth <- assessment_truth(scenario)
    for (i in seq_len(if (quick) 2L else 30L)) {
      rows[[length(rows) + 1L]] <- capture_attempt(scenario, 20000L + i,
        function(seed) assessment_trial(seed, truth))
    }
  }
  attempts <- do.call(rbind, rows)
  summary <- validation_summary(attempts)
  write.csv(attempts, file.path(output, "attempts.csv"), row.names = FALSE)
  write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()),
    paste("Source revision:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    "Exact working-source hashes:", capture.output(tools::md5sum(c(
      "R/aic.R", "R/nll.R", "analysis/information_criteria/validate.R")))),
    file.path(output, "session.txt"))
  print(summary, row.names = FALSE)
  invisible(list(attempts = attempts, summary = summary))
}

if (sys.nframe() == 0L) {
  pkgload::load_all(quiet = TRUE)
  run_validation(quick = "--quick" %in% commandArgs(TRUE))
}
