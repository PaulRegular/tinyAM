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
  future <- if (is.null(proj_settings)) numeric() else max(years) + 1
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

test_that("hindcast targets preserve raw plus ages and allow explicit overrides", {
  f <- hindcast_fixture()
  before <- f
  d <- .hindcast_controls(f, 2004:2005, "observed", NULL)
  expect_equal(d$target, c(18, 18))
  expect_identical(f, before)
  f$refit_args$data$catch$obs[1] <- NA
  d <- .hindcast_controls(f, 2000:2001, "observed", NULL)
  expect_match(d$error[1], "incomplete")
  expect_true(is.na(d$target[1]))
  expect_equal(.hindcast_controls(f, 2000:2001, c(`2000` = 12), NULL)$target, c(12, 18))
  expect_equal(.hindcast_controls(f, 2000:2001, 10, NULL)$target, c(10, 10))
  expect_equal(.hindcast_controls(f, 2000:2001, NULL, .8)$target, c(.8, .8))
  expect_error(.hindcast_controls(f, 2000:2001, c(1, 2), NULL), "Name")
  expect_error(.hindcast_controls(f, 2000:2001, "observed", 1), "not both")
})

test_that("hindcasts score indices, preserve inputs and retain the historical reference", {
  local_mocked_bindings(update = mock_hindcast_update, tidy_tam = mock_hindcast_tidy,
    check_tam = function(...) list(is_converged = TRUE))
  future::plan(future::sequential)
  f <- hindcast_fixture(); before <- f
  h <- fit_hindcast(f, folds = 2, progress = FALSE)
  expect_identical(f, before)
  expect_equal(h$hindcast_rmse, log(1.1))
  expect_equal(h$fold_diagnostics$n_scored_index, c(2, 2, 0))
  expect_equal(h$fold_diagnostics$target[1:2], c(18, 18))
  expect_null(h$fits[["2005"]]$used_projection)
  expect_false(any(h$fits[["2005"]]$dat$is_proj))
  expect_equal(h$fits[["2003"]]$used_projection, list(n_proj = 1, n_mean = 1, catch = 18))
  quo <- fit_hindcast(f, folds = 1, F_mult = .8, progress = FALSE)
  expect_equal(quo$fits[["2004"]]$used_projection$F_mult, .8)
  expect_error(fit_hindcast(f, folds = 0, catch = 10, F_mult = .8), "not both")
  retro <- fit_retro(f, folds = 1, progress = FALSE)
  expect_true(all(vapply(retro$fits, function(x) !any(x$dat$is_proj), logical(1))))
  expect_null(retro$hindcast_rmse)
})

test_that("missing targets and failed folds retain informative diagnostics", {
  local_mocked_bindings(update = mock_hindcast_update, tidy_tam = mock_hindcast_tidy,
    check_tam = function(...) list(is_converged = TRUE))
  future::plan(future::sequential)
  f <- hindcast_fixture()
  i <- which(f$refit_args$data$catch$year == 2004)[1]
  f$refit_args$data$catch$obs[i] <- NA_real_
  expect_warning(h <- fit_hindcast(f, folds = 2, progress = FALSE), "2003")
  expect_match(h$fold_diagnostics$error[1], "incomplete")
  expect_false(h$fold_diagnostics$is_converged[1])
  expect_equal(names(h$fits), c("2004", "2005"))
  supplied <- fit_hindcast(f, folds = 2, catch = c(`2004` = 18), progress = FALSE)
  expect_true(all(supplied$fold_diagnostics$is_converged))
  expect_equal(supplied$fold_diagnostics$target_source[1], "supplied catch")
  # Historical-only retros strip even a reference fit's inherited projections.
  f$dat$years <- 2000:2007; f$dat$is_proj <- c(rep(FALSE, 6), TRUE, TRUE)
  retro <- fit_retro(f, folds = 1, progress = FALSE)
  expect_false(any(retro$fits[["2005"]]$dat$is_proj))
})

test_that("hindcast scoring retains unavailable forecasts without an empty-table error", {
  data <- data.frame(year = 2001L, age = 2L, survey = "survey", pred = 10,
    obs = 11, is_proj = FALSE)
  paired <- tinyAM:::.hindcast_index_data(data, data)
  expect_equal(nrow(paired), 1L)
  expect_false(any(paired$is_proj))
  expect_true(is.nan(compute_hindcast_rmse(paired)))
  empty <- tinyAM:::.hindcast_index_data(data[FALSE, ], data[FALSE, ])
  expect_equal(nrow(empty), 0L)
  expect_true(is.nan(compute_hindcast_rmse(empty)))
})
