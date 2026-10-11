test_that("prepare_tam and fit_tam expose the same data and model arguments", {
  settings <- c("data", "years", "ages", "N_settings", "F_settings",
                "M_settings", "catch_settings", "index_settings", "proj_settings")
  expect_identical(names(formals(prepare_tam)), settings)
  expect_true(all(settings %in% names(formals(fit_tam))))
  expect_false(any(c("obs", "...") %in% names(formals(fit_tam))))
  expect_identical(as.list(formals(fit_tam))[settings], as.list(formals(prepare_tam)))
  expect_false("make_dat" %in% getNamespaceExports("tinyAM"))

  dat <- prepare_tam(data = cod_obs, years = 2000:2005, ages = 2:6)
  expect_true("obs" %in% names(dat))
  expect_false("data" %in% names(dat))
  expect_true(all(vapply(dat$obs, function(x) "obs" %in% names(x), logical(1))))
  expect_error(prepare_tam(obs = cod_obs), "unused argument")
  expect_error(fit_tam(obs = cod_obs), "unused argument")
})

test_that("check_tam_data is exactly the existing validator", {
  expect_identical(check_tam_data, check_obs)
  expect_invisible(check_tam_data(cod_obs))
  invalid <- cod_obs
  invalid$index$samp_time[1] <- 2
  expect_error(check_tam_data(invalid), "samp_time")
})

test_that("update stores evaluated data under the public argument name", {
  expect_identical(default_fit$refit_args$data, cod_obs)
  expect_false("obs" %in% names(default_fit$refit_args))
  new_data <- cod_obs
  new_data$index <- new_data$index[new_data$index$age <= 10, ]
  call <- update(default_fit, data = new_data, evaluate = FALSE)
  expect_identical(call$data, quote(new_data))
  expect_false("obs" %in% names(as.list(call)))
  local_mocked_bindings(fit_tam = function(...) list(refit_args = list(...)))
  fit <- update(default_fit, data = new_data, years = 2000:2005,
                ages = 2:6, silent = TRUE)
  expect_identical(fit$refit_args$data, new_data)
  expect_equal(fit$refit_args$years, 2000:2005)
  expect_equal(fit$refit_args$ages, 2:6)
})
