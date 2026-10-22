hindcast_fixture <- function() {
  d <- expand.grid(year = 2000:2005, age = 2:4)
  d$obs <- d$age
  w <- d; w$obs <- 2
  index <- d[d$age <= 3, ]; index$survey <- "Survey"
  index$obs <- exp(index$year - 2000 + index$age)
  index$pred <- index$obs; index$is_proj <- FALSE
  classed <- list(dat = list(years = 2000:2005, ages = 2:3, plus_ages = 3:4,
    is_proj = rep(FALSE, 6)), grad_tol = .01,
    refit_args = list(data = list(catch = d, weight = w)),
    obs_pred = list(index = index, catch = index), is_converged = TRUE)
  class(classed) <- c("tam_fit", "list")
  classed
}

mock_hindcast_update <- function(object, years, proj_settings, ...) {
  r <- object
  future <- if (is.null(proj_settings)) numeric() else max(years) + seq_len(proj_settings$n_proj)
  r$dat$years <- c(years, future)
  r$dat$is_proj <- c(rep(FALSE, length(years)), rep(TRUE, length(future)))
  r$used_projection <- proj_settings
  r$obs_pred <- lapply(object$obs_pred, function(d) {
    d <- d[d$year %in% c(years, future), ]
    d$is_proj <- d$year %in% future
    d$pred <- d$obs * 1.1
    d$obs[d$is_proj] <- NA_real_
    d
  })
  r
}
mock_hindcast_tidy <- function(model_list, ...) {
  stack <- function(type) do.call(rbind, lapply(names(model_list), function(n) {
    d <- model_list[[n]]$obs_pred[[type]]; d$fold <- as.numeric(n); d
  }))
  ssb <- do.call(rbind, lapply(names(model_list), function(n) {
    d <- model_list[[n]]$dat
    data.frame(year = d$years, is_proj = d$is_proj, est = 100, fold = as.numeric(n))
  }))
  list(obs_pred = list(index = stack("index"), catch = stack("catch")), pop = list(ssb = ssb))
}

test_that("hindcast controls preserve raw plus ages and recycle by horizon", {
  f <- hindcast_fixture(); before <- f
  expect_equal(.hindcast_controls(f, 2004:2005, .retro_projection(NULL))$target, c(18, 18))
  expect_identical(f, before)
  f$refit_args$data$catch$obs[1] <- NA
  expect_match(.hindcast_controls(f, 2000:2001, .retro_projection(NULL))$error[1], "incomplete")
  s <- .retro_projection(list(n_proj = 2, catch = c(10, 20)))
  expect_equal(.hindcast_controls(f, 2000:2001, s)$target, c(10, 20))
  expect_equal(.retro_projection(list(n_proj = 2, F_mult = .8))$F_mult, c(.8, .8))
  expect_error(.retro_projection(list(catch = 10, F_mult = 1)), "not both")
  expect_error(.retro_projection(list(n_proj = 3, catch = c(1, 2))), "n_proj")
})

test_that("single and multiple model hindcasts retain forecasts and diagnostics", {
  local_mocked_bindings(update = mock_hindcast_update, tidy_tam = mock_hindcast_tidy,
    check_tam = function(...) list(is_converged = TRUE))
  future::plan(future::sequential)
  f <- hindcast_fixture(); before <- f
  h <- fit_hindcast(f, folds = 2, progress = FALSE)
  expect_s3_class(h, "tam_retro")
  expect_identical(f, before)
  overall <- subset(h$scores, type == "index" & level == "overall" & is.na(horizon))
  expect_equal(overall$RMSE, log(1.1))
  expect_equal(overall$MASE, log(1.1))
  expect_equal(h$fold_diagnostics$n_scored_index, c(2, 2, 0))
  expect_equal(h$projection_controls$target, c(18, 18))
  expect_null(h$fits$f[["2005"]]$used_projection)
  expect_equal(h$fits$f[["2003"]]$used_projection, list(n_proj = 1, n_mean = 1, catch = 18))
  expect_output(print(h), "log_RMSE")
  multi <- fit_hindcast(f, same = f, folds = 2, progress = FALSE,
    proj_settings = list(n_proj = 3, n_mean = 3, F_mult = c(.8, 1, 1.2)))
  expect_equal(names(multi$fits), c("f", "same"))
  expect_true(multi$comparable)
  expect_equal(multi$fold_diagnostics$available_n_proj, rep(c(2, 1, 0), 2))
  expect_equal(multi$fits$f[["2003"]]$used_projection,
    list(n_proj = 2, n_mean = 3, F_mult = c(.8, 1)))
  expect_true(all(c(1, 2) %in% multi$forecast_pairs$index$horizon))
  expect_true(all(c("model", "fold", "horizon") %in% names(multi$pop$ssb)))
  expect_true(all(multi$forecast_pairs$index$baseline_year <= multi$forecast_pairs$index$fold))
  retro <- fit_retro(f, folds = 1, progress = FALSE)
  expect_true(all(vapply(retro$fits$f, function(x) !any(x$dat$is_proj), logical(1))))
  expect_equal(nrow(retro$scores), 0L)
  expect_error(fit_retro(f, proj_settings = list(n_mean = 3)), "only when")
  expect_error(fit_hindcast(f, F_mult = 1), "tam_fit")
})

test_that("missing targets and all-failed runs remain inspectable", {
  local_mocked_bindings(update = mock_hindcast_update, tidy_tam = mock_hindcast_tidy,
    check_tam = function(...) list(is_converged = TRUE))
  future::plan(future::sequential)
  f <- hindcast_fixture()
  f$refit_args$data$catch$obs[which(f$refit_args$data$catch$year == 2004)[1]] <- NA_real_
  expect_warning(h <- fit_hindcast(f, folds = 2, progress = FALSE), "failed")
  expect_match(h$fold_diagnostics$error[1], "incomplete")
  expect_null(h$fits$f[["2003"]])
  expect_true("2003" %in% names(h$fits$f))
  local_mocked_bindings(check_tam = function(...) list(is_converged = FALSE))
  expect_warning(h <- fit_hindcast(f, folds = 1, progress = FALSE), "failed")
  expect_equal(nrow(h$scores), 0L)
  expect_equal(nrow(h$fold_diagnostics), 2L)
  expect_output(print(h), "NA")
})

test_that("validation cases are compared without ranking", {
  local_mocked_bindings(update = mock_hindcast_update, tidy_tam = mock_hindcast_tidy,
    check_tam = function(...) list(is_converged = TRUE))
  future::plan(future::sequential)
  f <- hindcast_fixture(); other <- f
  other$obs_pred$index$obs[1] <- other$obs_pred$index$obs[1] * 2
  # Change an actual validation observation, not an unused historical observation.
  other$obs_pred$index$obs[other$obs_pred$index$year == 2005] <- 10
  expect_warning(h <- fit_hindcast(f, other = other, folds = 1, progress = FALSE), "not directly comparable")
  expect_false(h$comparable)
  expect_false("rank" %in% names(h$scores))
  a <- formals(fit_retro); b <- formals(fit_hindcast)
  a$hindcast <- b$hindcast <- NULL
  expect_identical(a, b)
})
