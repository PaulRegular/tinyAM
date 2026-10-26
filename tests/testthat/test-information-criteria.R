# A conditional mode at fixed outer values isolates criterion mathematics from
# optimizer performance. Integration tests below also use the fitted fixtures.
ic_test_fit <- function(..., map = list()) {
  args <- list(data = cod_obs, years = 1983:1992, ages = 2:5,
    N_settings = list(process = "off"),
    F_settings = list(process = "iid", mu_form = ~1))
  supplied <- list(...)
  args[names(supplied)] <- supplied
  d <- do.call(prepare_tam, args)
  p <- make_par(d)
  random <- c("log_f", "log_r", tinyAM:::.formula_random_parameters(d))
  if (d$N_settings$init == "latent") random <- c(random, "log_n0")
  if (d$N_settings$process != "off") random <- c(random, "log_n")
  if (d$M_settings$process != "off") random <- c(random, "log_m")
  o <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, random = random, map = map, silent = TRUE)
  value <- as.numeric(o$fn(o$par))
  structure(list(dat = d, obj = o, parameter_values = o$env$last.par,
    parameter_map = map, opt = list(par = o$par, objective = value),
    is_converged = FALSE), class = "tam_fit")
}

test_that("marginal likelihood and AIC use active outer parameters", {
  f <- ic_test_fit(map = list(log_mu_f = factor(NA)))
  ll <- suppressWarnings(logLik(f))
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -f$opt$objective)
  expect_equal(attr(ll, "df"), length(f$opt$par))
  expect_equal(attr(ll, "nobs"), sum(f$dat$is_observed))
  expect_equal(suppressWarnings(AIC(f)), 2 * f$opt$objective + 2 * length(f$opt$par))
  expect_equal(suppressWarnings(AIC(f, k = 3)), 2 * f$opt$objective + 3 * length(f$opt$par))
  multiple <- suppressWarnings(AIC(f, f))
  expect_named(multiple, c("df", "AIC"))
  expect_equal(multiple$AIC, rep(suppressWarnings(AIC(f)), 2))
  expect_error(AIC(f, list()), "tam_fit")
  expect_error(AIC(f, k = NA_real_), "finite")
})

test_that("Laplace marginal likelihood agrees with an exact Gaussian integral", {
  y <- c(-.3, .2, .6, 1.4)
  z <- diag(2)[c(1, 1, 2, 2), ]
  obj <- RTMB::MakeADFun(function(p) {
    -sum(RTMB::dnorm(y, p$beta + z %*% p$b, .4, log = TRUE)) -
      sum(RTMB::dnorm(p$b, 0, .6, log = TRUE))
  }, list(beta = 0, b = c(0, 0)), random = "b", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  covariance <- .6^2 * tcrossprod(z) + diag(.4^2, 4)
  residual <- y - opt$par[1]
  expected <- -.5 * (length(y) * log(2 * pi) +
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
    drop(crossprod(residual, solve(covariance, residual))))
  fit <- ic_test_fit()
  fit$obj <- obj
  fit$opt <- opt
  expect_equal(as.numeric(suppressWarnings(logLik(fit))), expected, tolerance = 1e-8)
  expect_equal(suppressWarnings(AIC(fit)), -2 * expected + 2, tolerance = 1e-8)
})

test_that("comparison checks observations and units, allowing row reordering", {
  f <- ic_test_fit()
  changed <- f
  changed$dat$observed[1] <- changed$dat$observed[1] + .01
  expect_error(AIC(f, changed), "same included")
  expect_error(cAIC(f, changed), "same included")
  changed <- f
  changed$dat$obs_map$unit <- "different units"
  expect_error(cAIC(f, changed), "units")
  changed <- f
  rows <- rev(seq_len(nrow(f$dat$obs_map)))
  changed$dat$obs_map <- f$dat$obs_map[rows, ]
  changed$dat$is_observed <- f$dat$is_observed[rows]
  changed$dat$observed <- f$dat$log_obs[rows][changed$dat$is_observed]
  expect_identical(tinyAM:::.aic_observations(f), tinyAM:::.aic_observations(changed))
})

test_that("sparse trace agrees with independently calculated Gaussian EDF", {
  s <- c(.2, .4, .6, .8, 1)
  covariance <- outer(s, s) * .4^abs(outer(1:5, 1:5, "-"))
  prior <- solve(covariance)
  z <- rbind(diag(5), diag(5), diag(5))
  joint <- prior + crossprod(z) / .7^2
  expected <- diag(solve(joint, prior))
  actual <- tinyAM:::.aic_trace(Matrix::Matrix(joint, sparse = TRUE), Matrix::Matrix(prior, sparse = TRUE))
  expect_equal(actual, expected, tolerance = 1e-10)
  # Prior entries beyond the joint sparsity pattern must still be included.
  joint <- diag(3, 5)
  joint[row(joint) == col(joint) + 1L | col(joint) == row(joint) + 1L] <- -.4
  prior <- diag(1, 5)
  prior[1, 5] <- prior[5, 1] <- .3
  actual <- tinyAM:::.aic_trace(Matrix::Matrix(joint, sparse = TRUE), Matrix::Matrix(prior, sparse = TRUE))
  expect_equal(actual, diag(solve(joint, prior)), tolerance = 1e-10)
  expect_equal(sum(1 - actual), 5 - sum(diag(solve(joint, prior))))
  prior[1, ] <- prior[, 1] <- 0
  expect_true(all(is.finite(tinyAM:::.aic_trace(Matrix::Matrix(joint, sparse = TRUE), Matrix::Matrix(prior, sparse = TRUE)))))
  expect_error(suppressWarnings(tinyAM:::.aic_trace(Matrix::Matrix(matrix(1, 5, 5), sparse = TRUE), Matrix::Diagonal(5))), "positive definite")
  expect_error(suppressWarnings(tinyAM:::.aic_trace(Matrix::Diagonal(5, -1), Matrix::Diagonal(5))), "positive definite")
  expect_warning(tinyAM:::.aic_trace(Matrix::Diagonal(2, c(1, 1e-14)), Matrix::Diagonal(2)), "numerically sensitive")
})

test_that("likelihood components preserve values and derivatives", {
  f <- ic_test_fit()
  process <- tinyAM:::.fitted_objective(f, "process")
  observation <- tinyAM:::.fitted_objective(f, "observation")
  p <- f$parameter_values
  expect_equal(f$obj$env$f(p), process$env$f(p) + observation$env$f(p))
  expect_equal(f$obj$env$f(p, order = 1),
    process$env$f(p, order = 1) + observation$env$f(p, order = 1), tolerance = 1e-8)
  expect_true(all(is.finite(observation$env$f(p, order = 1))))
  # Fixed F-mean and recruitment-SD coefficients only enter the process density.
  expected <- sum(names(f$opt$par) %in% c("log_r0", "log_q", "log_sd_catch", "log_sd_index"))
  expect_equal(tinyAM:::.aic_conditional_count(f, process, f$obj$env$lrandom()), expected)
})

test_that("cAIC reports reconciled EDF without changing the fitted object", {
  f <- ic_test_fit(index_settings = list(q_form = ~ factor(age) + rw(year, sd = .1)))
  before <- list(f$parameter_values, f$obj$env$last.par, f$obj$env$last.par.best, f$opt, f$dat)
  table <- suppressWarnings(cAIC(f, details = TRUE))
  expect_true(all(is.finite(unlist(table))))
  expect_equal(table$delta, 0)
  expect_gt(table$edf_other, 0)
  expect_equal(table$edf_random, unname(rowSums(table[paste0("edf_", c("rec", "N", "F", "M", "other"))])))
  expect_equal(table$EDF, table$edf_fixed + table$edf_random)
  expect_equal(table$cAIC, -2 * table$logLik + 2 * table$EDF)
  expect_equal(suppressWarnings(cAIC(f)), table$cAIC)
  expect_identical(before, list(f$parameter_values, f$obj$env$last.par, f$obj$env$last.par.best, f$opt, f$dat))
  # A later objective evaluation must not replace the stored fitted mode.
  f$obj$fn(f$opt$par + .02)
  later <- f$obj$env$last.par
  expect_equal(suppressWarnings(cAIC(f)), table$cAIC)
  expect_identical(f$obj$env$last.par, later)
  multiple <- suppressWarnings(cAIC(f, f))
  expect_named(multiple, c("EDF", "cAIC", "delta"))
  expect_equal(multiple$delta, c(0, 0))
  expect_error(cAIC(f, details = NA), "TRUE or FALSE")
  alternative <- ic_test_fit(index_settings = list(q_form = ~1))
  ranked <- suppressWarnings(cAIC(f, alternative))
  expect_true(all(diff(ranked$cAIC) >= 0))
  expect_equal(ranked$delta, ranked$cAIC - min(ranked$cAIC))
})

test_that("criterion warnings distinguish intrinsic and anchored walks", {
  f <- ic_test_fit()
  f$is_converged <- TRUE
  expect_no_warning(logLik(f)) # Anchored recruitment RW is proper.
  f$dat$F_settings$process <- "rw"
  expect_warning(logLik(f), "improper")
  f$dat$F_settings$process <- "iid"
  f$dat$N_settings$init <- "latent"
  expect_warning(logLik(f), "improper")
  f$dat$N_settings$init <- "exp"
  f$dat$N_settings$process <- "iid"
  expect_warning(tinyAM:::.aic_warnings(f, conditional = TRUE), "nonlinear")
  f$dat$N_settings$process <- "off"
  f$dat$rec$curve <- list(type = "bh")
  expect_warning(tinyAM:::.aic_warnings(f, conditional = TRUE), "nonlinear")
  f$dat$rec$curve <- NULL
  f$bounds <- list(lower = f$opt$par, upper = rep(Inf, length(f$opt$par)))
  expect_warning(logLik(f), "bound")
})

test_that("fitted N and M processes are included in conditional EDF", {
  expect_true(is.finite(suppressWarnings(AIC(default_fit))))
  n <- suppressWarnings(cAIC(default_fit, details = TRUE))
  m <- suppressWarnings(cAIC(M_dev, details = TRUE))
  expect_gt(n$edf_N, 0)
  expect_gt(m$edf_M, 0)
  expect_equal(n$n_fixed, length(default_fit$opt$par))
  expect_equal(m$n_random, sum(M_dev$obj$env$lrandom()))
})

test_that("latent initial abundance and fixed M outside process ages are counted", {
  f <- ic_test_fit(N_settings = list(process = "off", init = "latent"))
  table <- suppressWarnings(cAIC(f, details = TRUE))
  expect_gt(table$edf_other, 0)
  expect_equal(table$n_random, sum(f$obj$env$lrandom()))
  partial <- ic_test_fit(M_settings = list(process = "iid", mu_form = ~0 + factor(age),
    age_breaks = 3:5, first_dev_year = 1983))
  process <- tinyAM:::.fitted_objective(partial, "process")
  count <- tinyAM:::.aic_conditional_count(partial, process, partial$obj$env$lrandom())
  plain <- ic_test_fit()
  ordinary <- tinyAM:::.aic_conditional_count(plain, tinyAM:::.fitted_objective(plain, "process"), plain$obj$env$lrandom())
  expect_equal(count, ordinary + 1L) # M mean controls age 2, outside fitted blocks.
})

test_that("forecast-only states do not add artificial AIC or EDF", {
  historical <- ic_test_fit()
  projected <- ic_test_fit(proj_settings = list(n_proj = 3))
  expect_equal(suppressWarnings(AIC(projected)), suppressWarnings(AIC(historical)), tolerance = 1e-7)
  h <- suppressWarnings(cAIC(historical, details = TRUE))
  p <- suppressWarnings(cAIC(projected, details = TRUE))
  expect_gt(p$n_random, h$n_random)
  expect_equal(p$EDF, h$EDF, tolerance = 1e-7)
  expect_equal(p$cAIC, h$cAIC, tolerance = 1e-7)
})

test_that("fixed and shared maps change the appropriate counts", {
  f <- ic_test_fit()
  q <- length(f$obj$env$parList()$log_q)
  shared <- ic_test_fit(map = list(log_q = factor(rep(1, q))))
  a <- suppressWarnings(cAIC(f, details = TRUE))
  b <- suppressWarnings(cAIC(shared, details = TRUE))
  expect_equal(a$n_fixed - b$n_fixed, q - 1)
  expect_equal(a$edf_fixed - b$edf_fixed, q - 1)
  fixed <- ic_test_fit(map = list(log_mu_f = factor(NA)))
  c <- suppressWarnings(cAIC(fixed, details = TRUE))
  expect_equal(a$n_fixed - c$n_fixed, 1)
  expect_equal(a$edf_fixed, c$edf_fixed)
})
