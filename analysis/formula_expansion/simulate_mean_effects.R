args <- commandArgs(trailingOnly = TRUE)
n_simple <- if (length(args)) as.integer(args[1]) else 100L
n_full <- if (length(args) > 1L) as.integer(args[2]) else 30L
stopifnot(n_simple > 0, n_full > 0)
pkgload::load_all(quiet = TRUE)
root <- file.path("analysis", "formula_expansion")
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

process_formula <- function(type, component) {
  term <- call(type, as.name("year"))
  rhs <- if (component == "F") call("+", quote(factor(age)), term) else call("+", 0, term)
  stats::as.formula(call("~", rhs), env = environment())
}

interval_row <- function(estimate, se, truth, name, link = "log") {
  transform <- if (link == "log") exp else plogis
  data.frame(parameter = name, truth = truth, estimate = transform(estimate),
    lower = transform(estimate - qnorm(.975) * se),
    upper = transform(estimate + qnorm(.975) * se))
}

# Replicated log-mortality observations isolate the variance mathematics.
simple_recovery <- function(type, sparse = FALSE, seed) {
  set.seed(seed)
  d <- expand.grid(year = 2000:2039, age = seq_len(if (sparse) 2L else 8L))
  d$obs <- 1
  d$is_proj <- FALSE
  dat <- list(obs = list(catch = d))
  compiled <- tinyAM:::.parse_mean_formula(process_formula(type, "F"), dat, "F")
  dat$F_terms <- compiled$terms
  p <- tinyAM:::.q_term_parameters(dat$F_terms)
  term <- dat$F_terms[[1]]
  p[[term$sd_parameter]][] <- log(.18)
  if (type == "ar1") p[[term$phi_parameter]][] <- qlogis(.65)
  truth <- tinyAM:::.mean_effects(p, dat, "F", TRUE)
  y <- log(.2) + truth$contribution + rnorm(nrow(d), 0, sqrt(.1^2 + .05^2))
  p$baseline <- log(.2)
  p$log_residual_sd <- log(.1)
  obj <- RTMB::MakeADFun(function(p) {
    e <- tinyAM:::.mean_effects(p, dat, "F")
    e$nll - sum(RTMB::dnorm(y, p$baseline + e$contribution,
      sqrt(exp(2 * p$log_residual_sd) + .05^2), log = TRUE))
  }, p, random = term$parameter, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 1000))
  obj$fn(opt$par)
  sdr <- RTMB::sdreport(obj, getReportCovariance = FALSE)
  fixed <- sdr$par.fixed
  se <- sqrt(diag(sdr$cov.fixed))
  tab <- rbind(interval_row(fixed[term$sd_parameter], se[term$sd_parameter], .18, "mean_sd"),
    interval_row(fixed["log_residual_sd"], se["log_residual_sd"], .1, "residual_sd"))
  if (type == "ar1") tab <- rbind(tab, interval_row(fixed[term$phi_parameter],
    se[term$phi_parameter], .65, "phi", "logit"))
  tab$type <- type
  tab$design <- if (sparse) "two ages" else "eight ages"
  tab$optimizer <- opt$convergence
  tab$gradient <- max(abs(obj$gr(opt$par)))
  tab$pdHess <- sdr$pdHess
  list(table = tab, trajectory = data.frame(year = d$year[1:40],
    truth = truth$contribution[1:40], estimate = tinyAM:::.mean_effects(
      obj$env$parList(), dat, "F")$contribution[1:40], type = type,
    design = tab$design[1]))
}

make_observations <- function() {
  grid <- expand.grid(year = 2000:2029, age = 1:8)
  catch <- transform(grid, obs = 1, relative_sd = .08)
  index <- rbind(transform(grid, survey = "early", samp_time = .25),
                 transform(grid, survey = "late", samp_time = .75))
  index$survey <- factor(index$survey)
  index$obs <- 1
  index$relative_sd <- .08
  weight <- transform(grid, obs = .2 * age^1.5)
  maturity <- transform(grid, obs = plogis(age - 4))
  list(catch = catch, index = index, weight = weight, maturity = maturity)
}

# Full assessments estimate F, M, recruitment and q from simulated observations.
full_recovery <- function(component, type, seed, residual = TRUE, mean_sd = .12) {
  set.seed(seed)
  settings <- list(N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "iid", mu_form = ~ factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.25)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd, fill_missing = FALSE),
    index_settings = list(q_form = ~ survey, sd_form = ~ 0,
      sd_supplied = ~ relative_sd, fill_missing = FALSE))
  settings[[paste0(component, "_settings")]]$mu_form <- process_formula(type, component)
  if (component == "M" && residual) {
    settings$M_settings$process <- "iid"
    settings$M_settings$age_breaks <- 1:8
    settings$M_settings$first_dev_year <- 2000
  }
  dat <- do.call(prepare_tam, c(list(data = make_observations()), settings))
  p <- make_par(dat)
  term <- dat[[paste0(component, "_terms")]][[1]]
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.12)
  p$log_sd_f <- log(.08)
  p$log_mu_f[] <- c(log(.15), .09 * seq_len(7))
  p$log_f[] <- matrix(drop(dat$F_modmat %*% p$log_mu_f), 30, 8)
  if (component == "M" && residual) p$log_sd_m <- log(.06)
  p$log_q[] <- c(log(.7), log(.9 / .7))
  p[[term$sd_parameter]][] <- log(mean_sd)
  if (type == "ar1") p[[term$phi_parameter]][] <- qlogis(.65)
  simulated <- nll_fun(p, dat, simulate = TRUE)
  keep <- intersect(names(p), names(simulated))
  p[keep] <- simulated[keep]
  obs <- dat$obs
  values <- split(exp(simulated$log_obs), dat$obs_map$type)
  obs$catch$obs <- values$catch
  obs$index$obs <- values$index
  start <- p
  start[[term$sd_parameter]][] <- log(.1)
  if (type == "ar1") start[[term$phi_parameter]][] <- qlogis(.5)
  start[[term$parameter]][] <- 0
  start$log_sd_f <- log(.1)
  if (component == "M" && residual) start$log_sd_m <- log(.1)
  fit <- suppressWarnings(do.call(fit_tam, c(list(data = obs, silent = TRUE,
    start_par = start, grad_tol = .01), settings)))
  fixed <- fit$fixed_par
  select <- function(parameter, truth, label) {
    tab <- fixed[fixed$par == sub("^log(it)?_", "", parameter), ]
    if (!nrow(tab)) tab <- fixed[fixed$par == parameter, ]
    data.frame(parameter = label, truth = truth, estimate = tab$est,
      lower = tab$lwr, upper = tab$upr)
  }
  tab <- select(term$sd_parameter, mean_sd, "mean_sd")
  if (component == "F" || residual) tab <- rbind(tab,
    select(if (component == "F") "log_sd_f" else "log_sd_m",
           if (component == "F") .08 else .06, "residual_sd"))
  if (type == "ar1") tab <- rbind(tab, select(term$phi_parameter, .65, "phi"))
  tab$component <- component
  tab$type <- type
  tab$design <- if (component == "M" && !residual) "mean only" else "mean + IID residuals"
  if (mean_sd == .3) tab$design <- "strong mean + IID residuals"
  tab$optimizer <- fit$opt$convergence
  tab$gradient <- fit$diagnostics$max_gradient
  tab$pdHess <- isTRUE(fit[["sdrep"]]$pdHess)
  tab$converged <- fit$is_converged
  tab$message <- fit$opt$message
  # Preserve one example, not every fitted assessment.
  actual <- tinyAM:::.mean_effects(p, dat, component)$contribution[1:30]
  estimated <- fit$formula_effects$levels[[term$id]]
  list(table = tab, trajectory = data.frame(year = 2000:2029, truth = actual,
    estimate = estimated$est, lower = estimated$lwr, upper = estimated$upr,
    component = component, type = type, design = tab$design[1]), fit = fit)
}

safe_run <- function(fun, ...) {
  tryCatch(fun(...), error = function(e) list(error = conditionMessage(e)))
}

simple <- list()
full <- list()
failures <- list()
examples <- list()
for (type in c("iid", "rw", "ar1")) for (sparse in c(FALSE, TRUE)) {
  for (i in seq_len(n_simple)) {
    result <- safe_run(simple_recovery, type, sparse, 71000 + i +
      match(type, c("iid", "rw", "ar1")) * 1000 + sparse * 10000)
    key <- paste("simple", type, sparse, i, sep = "_")
    if (!is.null(result$error)) failures[[key]] <- result$error else {
      result$table$replicate <- i
      simple[[key]] <- result$table
      if (i == 1L) examples[[key]] <- result$trajectory
    }
  }
  message("Standalone: ", type, "; sparse = ", sparse)
}
for (component in c("F", "M")) for (type in c("iid", "rw", "ar1")) {
  residuals <- if (component == "M") c(TRUE, FALSE) else TRUE
  for (residual in residuals) for (i in seq_len(n_full)) {
    result <- safe_run(full_recovery, component, type,
      93000 + i + match(type, c("iid", "rw", "ar1")) * 1000 +
        (component == "M") * 10000 + !residual * 20000, residual)
    key <- paste("full", component, type, residual, i, sep = "_")
    if (!is.null(result$error)) failures[[key]] <- result$error else {
      result$table$replicate <- i
      full[[key]] <- result$table
      if (i == 1L) {
        examples[[key]] <- result$trajectory
        if (type == "rw" && residual) {
          saveRDS(result$fit, file.path(root, "results", paste0(key, ".rds")))
        }
      }
    }
    message("Full model: ", key, if (!is.null(result$error)) paste("ERROR", result$error))
  }
}
for (type in c("iid", "rw", "ar1")) for (i in seq_len(n_full)) {
  key <- paste("full_M", type, "strong", i, sep = "_")
  result <- safe_run(full_recovery, "M", type,
    163000 + i + match(type, c("iid", "rw", "ar1")) * 1000, TRUE, .3)
  if (!is.null(result$error)) failures[[key]] <- result$error else {
    result$table$replicate <- i
    full[[key]] <- result$table
    if (i == 1L) examples[[key]] <- result$trajectory
  }
  message("Full model: ", key, if (!is.null(result$error)) paste("ERROR", result$error))
}
saveRDS(list(simple = simple, full = full, failures = failures, examples = examples,
  n_simple = n_simple, n_full = n_full, session = sessionInfo()),
  file.path(root, "results", "mean_effects_recovery.rds"))
