test_that("whitening matches every process likelihood and stationary boundary", {
  x <- matrix(sin(1:20), 5, 4, dimnames = list(2000:2004, 2:5))
  sd <- c(.2, .4, .3, .6)
  for (process in c("iid", "rw", "cor_rw", "ar1")) {
    for (phi in list(c(0, .6), c(.4, 0), c(.4, .6))) {
      density <- function(x) {
        if (process == "cor_rw") dprocess_cor_rw(x, sd, .5) else
          .dprocess_scaled(x, sd, process, phi)
      }
      z <- .process_whiten(x, sd, process, phi, rho = .5)
      expect_equal(sum(z^2) / 2, as.numeric(density(x * 0) - density(x)), tolerance = 1e-10)
      expect_equal(nrow(z), nrow(x) - (process %in% c("rw", "cor_rw")))
    }
  }
  expect_equal(.process_whiten(x[, 1, drop = FALSE], .2, "ar1", c(0, .6))[1, 1], x[1, 1] * .8 / .2)
})

.diagnostic_fixture <- function(process = "ar1", rec = "rw", start = 1986,
                                curve = FALSE, mapped = FALSE, age_sd = FALSE) {
  d <- prepare_tam(cod_obs, years = 1983:1991, ages = 2:5,
    N_settings = list(process = process, sd_form = if (age_sd) ~I(age >= 4) else ~1),
    rec_settings = list(process = rec, mu_form = if (curve) ~bh(ssb, lag = 2) else if (rec == "rw") NULL else ~1),
    F_settings = list(process = if (process == "rw") "cor_rw" else process,
      mu_form = if (process == "rw") NULL else ~1, sd_form = if (age_sd) ~I(age >= 4) else ~1),
    M_settings = list(process = process, first_dev_year = start, age_breaks = c(2, 4, 5),
      sd_form = if (age_sd) ~I(age >= 4) else ~1),
    proj_settings = list(n_proj = 2))
  p <- make_par(d)
  p$log_f[] <- sin(seq_along(p$log_f)) / 4
  p$log_n[] <- p$log_n + cos(seq_along(p$log_n)) / 10
  p$log_m[] <- log(.2) + cos(seq_along(p$log_m)) / 4
  if (age_sd) for (component in c("n", "f", "m")) p[[paste0("sd_beta_", component)]] <- c(log(.3), .7)
  map <- if (mapped) list(log_sd_n = factor(NA)) else list()
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, map = map, silent = TRUE)
  mode <- obj$par
  Q <- Matrix::Diagonal(length(mode), x = 100)
  Q <- Matrix::forceSymmetric(methods::as(Q, "CsparseMatrix"))
  dimnames(Q) <- list(names(mode), names(mode))
  structure(list(dat = d, obj = obj, parameter_values = mode, parameter_map = map,
    rep = obj$report(mode), sdrep = list(jointPrecision = Q)), class = "tam_fit")
}

test_that("joint perturbations have analytic covariance and restore RNG", {
  Q <- Matrix::forceSymmetric(Matrix::Matrix(matrix(c(4, 1, 0, 1, 3, .8, 0, .8, 2), 3)))
  set.seed(792)
  before <- .Random.seed
  J <- matrix(c(1, 0, 2, -1, .5, 1), 2)
  draws <- vapply(1:4000, function(seed) as.numeric(J %*% .process_perturbation(Q, seed)), numeric(2))
  expect_identical(.Random.seed, before)
  expect_equal(stats::cov(t(draws)), as.matrix(J %*% solve(Q) %*% t(J)), tolerance = .045)
  expect_identical(.process_perturbation(Q, 123456), .process_perturbation(Q, 123456))
  rm(".Random.seed", envir = .GlobalEnv)
  .process_perturbation(Q, 2)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("AD joint propagation agrees with directional differences of shared predictions", {
  fit <- .diagnostic_fixture(curve = TRUE, rec = "ar1", age_sd = TRUE)
  Q <- fit$sdrep$jointPrecision * 1e10
  fit$sdrep$jointPrecision <- Q
  sampled <- .tidy_process_residuals(fit, 721L)
  expect_equal(sampled$status$status, "ok")
  delta <- .process_perturbation(Q, 721L)
  dat <- fit$dat; dat$.likelihood_component <- "process_innovations"
  g <- function(values) nll_fun(fit$obj$env$parList(par = values), dat)
  directional <- as.numeric((g(fit$parameter_values + delta) - g(fit$parameter_values - delta)) / 2)
  fitted <- unlist(lapply(sampled[c("recruitment", "N", "F", "M")], `[[`, "fitted"), use.names = FALSE)
  residual <- unlist(lapply(sampled[c("recruitment", "N", "F", "M")], `[[`, "est"), use.names = FALSE)
  expect_equal(residual - fitted, directional, tolerance = 1e-8)
  withr::with_preserve_seed({
    RNGkind("L'Ecuyer-CMRG"); set.seed(72)
    state <- .Random.seed; kind <- RNGkind()
    expect_equal(sampled, .tidy_process_residuals(fit, 721L))
    expect_identical(.Random.seed, state)
    expect_identical(RNGkind(), kind)
  })
})

test_that("shared predictions handle recruitment, plus groups, M start and maps", {
  for (process in c("iid", "rw", "ar1")) for (rec in c("iid", "rw", "ar1")) {
    fit <- .diagnostic_fixture(process, rec, curve = rec != "rw", mapped = TRUE)
    mode <- fit$obj$env$last.par
    uncertainty <- serialize(fit$sdrep, NULL)
    x <- .tidy_process_residuals(fit, 19L)
    expect_equal(x$status$status, "ok")
    expect_identical(fit$obj$env$last.par, mode)
    expect_identical(serialize(fit$sdrep, NULL), uncertainty)
    expect_true(all(x$M$year >= 1986 + (process == "rw")))
    expect_equal(length(unique(x$M$age_block)), ncol(fit$obj$env$parList(par = fit$parameter_values)$log_m))
    expect_false(any(x$M$is_proj))
    expect_false(any(x$F$is_proj))
    expect_true(all(x$N$age %in% fit$dat$ages[-1]))
    expect_equal(min(x$N$year), 1984 + (process == "rw"))
    p <- fit$obj$env$parList(par = fit$parameter_values)
    pop <- .population_states(p, fit$dat, fit$rep$Z, F = fit$rep$F, M = fit$rep$M)
    i <- fit$dat$rec$eligible
    s <- .rec_scale(p, fit$dat)
    u <- pop$log_recruitment - pop$log_mu_R
    expected <- if (rec == "rw") (u[i] - u[i - 1L]) / s$sd else if (rec == "iid") u[i] / s$sd else
      c(u[i[1]] * sqrt(1 - s$phi^2), u[i[-1]] - s$phi * u[i[-length(i)]]) / s$sd
    expect_equal(x$recruitment$fitted, unname(expected[!fit$dat$is_proj[i]]))
    expect_equal(x, .tidy_process_residuals(fit, 19L))
  }
})

test_that("invalid uncertainty is explicit and never replaced by fitted innovations", {
  fit <- .diagnostic_fixture()
  fit$sdrep <- NULL
  expect_warning(x <- .tidy_process_residuals(fit, 2), "unavailable")
  expect_named(x, "status")
  expect_equal(x$status$status, "unavailable")
  expect_true(length(.tidy_process_deviations(fit)) > 0L)
  fit <- .diagnostic_fixture()
  names(fit$parameter_values)[1] <- "wrong"
  expect_warning(.tidy_process_residuals(fit, 2), "ordering")
  fit <- .diagnostic_fixture()
  fit$sdrep$jointPrecision <- -fit$sdrep$jointPrecision
  expect_warning(x <- .tidy_process_residuals(fit, 2), "unavailable")
  expect_equal(x$status$status, "unavailable")
  expect_equal(.tidy_process_residuals(fit, 2, FALSE)$status$status, "disabled")
  expect_equal(.tidy_process_residuals(structure(list(), class = "tam_ref"), 2)$status$status, "reference")
  for (seed in list(NA, -1, 1.2, "a", c(1, 2))) expect_error(.validate_process_seed(seed), "integer")
})

test_that("tidy sampling is model-order independent and opt-out preserves outputs", {
  skip_if(is.null(default_fit))
  one <- tidy_tam(model_list = list(A = default_fit), process_seed = 42)
  two <- tidy_tam(model_list = list(B = default_fit, A = default_fit), process_seed = 42)
  for (component in c("recruitment", "N", "F")) {
    x <- two$process_residuals[[component]]
    expect_equal(one$process_residuals[[component]]$est, x$est[x$model == "A"])
  }
  off <- tidy_tam(default_fit, process_residuals = FALSE)
  on <- tidy_tam(default_fit)
  expect_equal(off[names(off) != "process_residuals"], on[names(on) != "process_residuals"])
  expect_equal(off$process_residuals$status$status, "disabled")
})
