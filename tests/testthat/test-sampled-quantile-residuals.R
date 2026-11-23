.sampled_quantile_fixture <- function(stochastic = TRUE, mixed = TRUE) {
  data <- cod_obs
  for (name in names(data)) {
    if ("age" %in% names(data[[name]])) data[[name]] <- subset(data[[name]], age <= 5)
  }
  data$catch$obs[data$catch$age == 2] <- 0
  data$catch$obs[data$catch$year == 1992 & data$catch$age == 4] <- NA_real_
  data$index$obs[data$index$year == 1992 & data$index$age == 3] <- 0
  data$index$obs[data$index$year == 1993 & data$index$age == 4] <- NA_real_
  dat <- prepare_tam(data, years = 1991:1994, ages = 2:5,
    N_settings = list(process = "off"),
    rec_settings = list(process = if (stochastic) "iid" else "off",
      sd = if (stochastic) 0.3 else NULL, mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~1, age_breaks = 3:6),
    catch_settings = list(cor_struct = "ar1"),
    index_settings = list(q_form = ~1, family = if (mixed) delta_lognormal(type = "poisson-link", y50 = 100) else lognormal()),
    proj_settings = list(n_proj = 1, n_mean = 1, F_mult = 1))
  parameters <- make_par(dat)
  map <- list(log_mu_f = factor(NA))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), parameters, map = map,
    random = if (stochastic) "log_r" else NULL, silent = TRUE)
  obj$fn(obj$par)
  mode <- obj$env$last.par
  random <- obj$env$lrandom()
  Q <- Matrix::Diagonal(length(mode), x = 20)
  Q <- Matrix::forceSymmetric(methods::as(Q, "CsparseMatrix"))
  dimnames(Q) <- list(names(mode), names(mode))
  if (sum(random) > 1L) {
    ids <- which(random)[1:2]
    Q[ids[1L], ids[2L]] <- Q[ids[2L], ids[1L]] <- 5
  }
  fit <- structure(list(dat = dat, obj = obj, parameter_values = mode, parameter_map = map,
    rep = obj$report(mode), sdrep = list(jointPrecision = Q, pdHess = TRUE,
      par.fixed = mode[!random], par.random = mode[random])), class = "tam_fit")
  fit$obs_pred <- tidy_obs_pred(fit, quantile_method = "mle-eb")
  fit
}

test_that("draw-based quantile residuals hold fixed coefficients at their estimates", {
  fit <- .sampled_quantile_fixture()
  ref <- residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 712L)
  random <- fit$obj$env$lrandom()
  q <- fit$sdrep$jointPrecision
  q[!random, !random] <- q[!random, !random] * 7
  fixed_id <- which(!random)[1L]
  random_id <- which(random)[1L]
  q[fixed_id, random_id] <- q[random_id, fixed_id] <- 0.5
  fit$sdrep$jointPrecision <- q
  expect_identical(residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 712L), ref)

  # At infinite conditional precision the state draw approaches the fitted mode.
  fit$sdrep$jointPrecision <- q * 1e20
  near <- residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 712L)
  eb <- residuals(fit, type = "quantile", method = "mle-eb", residual_seed = 712L)
  expect_equal(near$residual, eb$residual, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(ref$residual, eb$residual)))
})

test_that("a conditional draw retains observation identity, gaps and the supplied fit", {
  fit <- .sampled_quantile_fixture()
  local_mocked_bindings(fit_tam = function(...) stop("A diagnostic must not refit"))
  state <- serialize(list(mode = fit$parameter_values, uncertainty = fit$sdrep,
    observations = fit$obs_pred, last = fit$obj$env$last.par,
    best = fit$obj$env$last.par.best), NULL)
  before <- fit$obj$env$last.par
  withr::with_preserve_seed({
    RNGkind("L'Ecuyer-CMRG")
    set.seed(41)
    rng <- .Random.seed
    kind <- RNGkind()
    r <- residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 593L)
    expect_identical(r, residuals(fit, type = "quantile", method = "mle-mvn",
      residual_seed = 593L))
    expect_identical(.Random.seed, rng)
    expect_identical(RNGkind(), kind)
  })
  another <- residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 594L)
  expect_false(isTRUE(all.equal(r$residual, another$residual)))
  expect_identical(fit$obj$env$last.par, before)
  expect_identical(serialize(list(mode = fit$parameter_values, uncertainty = fit$sdrep,
    observations = fit$obs_pred, last = fit$obj$env$last.par,
    best = fit$obj$env$last.par.best), NULL), state)

  map <- fit$dat$obs_map
  expect_equal(r[c("type", "year", "age")], map[c("type", "year", "age")], ignore_attr = TRUE)
  expect_identical(is.finite(r$residual), as.logical(map$is_observed))
  expect_true(any(map$is_structural))
  expect_true(any(fit$dat$obs_correlation$gap > 1))
  zero <- which(map$type == "index" & fit$dat$is_observed & fit$dat$obs_values == 0)
  expect_gt(length(zero), 0)
  expect_true(all(is.finite(r$residual[zero])))

  retained <- add_residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = 593L)
  metadata <- attr(retained$obs_pred, "quantile", exact = TRUE)
  expect_identical(metadata$method, "mle-mvn")
  expect_equal(metadata$seed, 593L)
  expect_equal(r$residual, residuals(retained, type = "quantile")$residual)
  tidied <- tidy_obs_pred(fit, quantile_method = "mle-mvn", residual_seed = 593L)
  expect_equal(tidied$catch$quantile_res, retained$obs_pred$catch$quantile_res)
  expect_equal(tidied$index$quantile_res, retained$obs_pred$index$quantile_res)
  expect_equal(r$residual, unlist(lapply(retained$obs_pred, function(d) d$quantile_res),
    use.names = FALSE))
  for (name in names(retained$obs_pred)) {
    expect_equal(retained$obs_pred[[name]]$pred, fit$obs_pred[[name]]$pred)
    expect_equal(retained$obs_pred[[name]]$std_res, fit$obs_pred[[name]]$std_res)
    expect_true(all(unname(retained$obs_pred[[name]]$quantile_method[
      is.finite(retained$obs_pred[[name]]$quantile_res)]) == "mle-mvn"))
  }
})

test_that("the explicit sampled method rejects unavailable or misaligned state uncertainty", {
  fit <- .sampled_quantile_fixture()
  bad <- fit
  bad$sdrep$jointPrecision <- NULL
  expect_error(residuals(bad, type = "quantile", method = "mle-mvn"), "precision|uncertainty")
  bad$sdrep$jointPrecision <- -fit$sdrep$jointPrecision
  expect_error(residuals(bad, type = "quantile", method = "mle-mvn"), "precision|Cholesky|positive")
  bad <- fit
  rownames(bad$sdrep$jointPrecision)[1L] <- "mismatched"
  expect_error(residuals(bad, type = "quantile", method = "mle-mvn"), "precision|ordering")
  expect_error(residuals(fit, type = "quantile", method = "mle-mvn", residual_seed = -1), "seed")
  expect_error(residuals(fit, type = "quantile", method = "undocumented"), "method")

  deterministic <- .sampled_quantile_fixture(stochastic = FALSE, mixed = FALSE)
  expect_false(any(deterministic$obj$env$lrandom()))
  expect_error(residuals(deterministic, type = "quantile", method = "mle-mvn"), "random")
  ordinary <- residuals(deterministic, type = "standardized")
  expect_equal(residuals(deterministic, type = "quantile", method = "mle-eb")$residual,
    ordinary$residual)
})

test_that("ordinary tidying selects a supported method and marks uncertainty fallback", {
  fit <- .sampled_quantile_fixture()
  fit$obs_pred <- NULL
  auto <- tidy_obs_pred(fit, residual_seed = 593L)
  metadata <- attr(auto, "quantile", exact = TRUE)
  expect_identical(metadata$requested, "auto")
  expect_identical(metadata$method, "mle-mvn")
  expect_identical(metadata$status, "ok")
  expect_equal(auto$index$quantile_res, tidy_obs_pred(fit, quantile_method = "mle-mvn",
    residual_seed = 593L)$index$quantile_res)

  fit$sdrep$jointPrecision <- NULL
  fallback <- tidy_obs_pred(fit, residual_seed = 593L)
  metadata <- attr(fallback, "quantile", exact = TRUE)
  expect_identical(metadata$requested, "auto")
  expect_identical(metadata$method, "mle-eb")
  expect_identical(metadata$status, "fallback")
  expect_match(metadata$detail, "precision|uncertainty")
  expect_equal(fallback$index$quantile_res, tidy_obs_pred(fit, quantile_method = "mle-eb",
    residual_seed = 593L)$index$quantile_res)

  fit <- .sampled_quantile_fixture(stochastic = FALSE, mixed = FALSE)
  fit$obs_pred <- NULL
  deterministic <- tidy_obs_pred(fit)
  metadata <- attr(deterministic, "quantile", exact = TRUE)
  expect_identical(metadata$requested, "auto")
  expect_identical(metadata$method, "mle-eb")
  expect_identical(metadata$status, "ok")
  expect_equal(deterministic$index$quantile_res, deterministic$index$std_res)
})
