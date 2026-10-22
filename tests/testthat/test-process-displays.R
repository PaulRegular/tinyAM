test_that("residual extraction preserves rows, masks and fitted objects", {
  d <- data.frame(year = 2000:2003, age = 2, obs = c(2, 0, NA, 3),
    pred = 1, std_res = c(.5, NA, NA, 1), is_proj = c(FALSE, FALSE, FALSE, TRUE))
  fit <- structure(list(obs_pred = list(catch = d, index = transform(d, survey = "S"))), class = "tam_fit")
  original <- fit
  x <- residuals(fit)
  expect_equal(nrow(x), 8L)
  expect_equal(x$residual[c(1, 5)], c(.5, .5))
  expect_true(all(is.na(x$residual[c(2:4, 6:8)])))
  expect_setequal(x$type, c("catch", "index"))
  expect_equal(residuals(fit, type = "response")$residual[1], 1)
  expect_equal(residuals(fit, type = "log")$residual[1], log(2))
  expect_error(residuals(fit, type = "osa"), "not been computed")
  fit$obs_pred$catch$osa_res <- rep(.3, 4)
  expect_equal(residuals(fit, type = "osa")$residual[1], .3)
  fit$obs_pred$catch$osa_res <- NULL
  expect_identical(fit, original)
})

test_that("deviations use absolute states, unique M blocks and RW increments", {
  for (process in c("iid", "ar1", "rw")) {
    d <- prepare_tam(cod_obs, years = 1983:1990, ages = 2:5,
      N_settings = list(process = process),
      F_settings = list(process = process, mu_form = if (process == "rw") NULL else ~1),
      M_settings = list(process = process, age_breaks = c(2, 4, 5)))
    p <- make_par(d)
    p$log_f[] <- sin(seq_along(p$log_f)) / 4
    p$log_m[] <- log(.2) + cos(seq_along(p$log_m)) / 4
    o <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    r <- o$report(o$par)
    fit <- structure(list(dat = d, rep = r, obj = o, parameter_values = o$par), class = "tam_fit")
    mode <- o$env$last.par
    values <- .tidy_process_deviations(fit)
    expected <- p$log_f
    if (process == "rw") expected <- expected[-1, ] - expected[-nrow(expected), ]
    expect_equal(values$F$est, as.vector(expected), ignore_attr = TRUE)
    expect_equal(nrow(values$M), (nrow(p$log_m) - (process == "rw")) * ncol(p$log_m))
    expect_true("age_block" %in% names(values$M))
    expect_false("age" %in% names(values$M))
    expect_true(all(values$N$age %in% d$ages[-1]))
    if (process == "rw") expect_equal(min(values$N$year), 1985)
    expect_identical(o$env$last.par, mode)
  }
})

test_that("process bubbles use signed values and explicit M block axes", {
  d <- data.frame(year = 2000:2003, age_block = c("2-4", "2-4", "5-8", "5-8"),
    est = c(-1, 1, -.5, .5), kind = "Log deviation from mean")
  original <- d
  p <- plotly::plotly_build(plot_bubbles(d, res = ~est, y = ~age_block, ylab = "M state block"))
  expect_equal(p$x$layout$yaxis$title, "M state block")
  expect_equal(p$x$layout$title, "Process deviations")
  expect_identical(d, original)
  expect_length(p$x$data, 2L)
  expect_true(any(grepl("-1", unlist(lapply(p$x$data, `[[`, "text")))))
})

test_that("mean F and initial N do not partially match absent process states", {
  d <- prepare_tam(cod_obs, years = 1983:1990, ages = 2:5,
    N_settings = list(process = "off", init = "free"),
    F_settings = list(process = "ar1", mu_form = ~logistic(age, by = F_y_block)))
  p <- make_par(d)
  o <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  fit <- structure(list(dat = d, rep = o$report(o$par), obj = o, parameter_values = o$par), class = "tam_fit")
  x <- .tidy_process_deviations(fit)
  expect_named(x, "F")
  expect_equal(nrow(x$F), length(d$years) * length(d$ages))
})
