typed_weights <- function(data, multiplier = 2) {
  stock <- data$weight
  stock$type <- "stock"
  catch <- stock
  catch$type <- "catch"
  catch$obs <- stock$obs * multiplier
  data$weight <- rbind(catch, stock)
  data
}

test_that("typed weights require complete valid weights and a common unit", {
  d <- typed_weights(cod_obs)
  expect_true(check_tam_data(d))
  d$weight <- d$weight[sample(nrow(d$weight)), ]
  expect_true(check_tam_data(d))
  for (type in c("invalid", NA_character_)) {
    bad <- d; bad$weight$type[1] <- type
    expect_error(check_tam_data(bad), "weight\\$type")
  }
  bad <- d; bad$weight <- subset(bad$weight, type == "catch")
  expect_error(check_tam_data(bad), "Stock weights are required")
  bad <- d; bad$weight <- bad$weight[-1, ]
  expect_error(check_tam_data(bad), "coverage is incomplete")
  bad <- d; bad$weight <- rbind(bad$weight, bad$weight[1, ])
  expect_error(check_tam_data(bad), "duplicate")
  bad <- d; bad$weight$obs[1] <- NA_real_
  expect_error(check_tam_data(bad), "NA")
  bad <- d; bad$weight$obs[1] <- -1
  expect_error(check_tam_data(bad), "non-negative")
  bad <- d; bad$weight$unit[bad$weight$type == "catch"] <- "g"
  expect_error(check_tam_data(bad), "unit")
})

test_that("catch weights leave historical numbers likelihood and stock biology unchanged", {
  untyped <- make_test_dat()
  stock <- cod_obs; stock$weight$type <- "stock"
  stock <- make_test_dat(data = stock)
  both <- make_test_dat(data = typed_weights(cod_obs, 1))
  changed <- make_test_dat(data = typed_weights(cod_obs, 2))
  expect_equal(untyped$W, stock$W)
  expect_equal(untyped$W, both$W_catch)
  expect_equal(untyped$W, changed$W)
  expect_equal(changed$W_catch, 2 * changed$W)
  par <- make_par(untyped)
  objective <- function(d) RTMB::MakeADFun(function(p) nll_fun(p, d), par, silent = TRUE)
  a <- objective(untyped); b <- objective(changed)
  expect_equal(a$fn(a$par), b$fn(b$par), tolerance = 1e-10)
  expect_equal(a$gr(a$par), b$gr(b$par), tolerance = 1e-10)
  ra <- a$report(); rb <- b$report()
  expect_equal(rb$N, ra$N)
  expect_equal(rb$ssb, ra$ssb)
  expect_equal(rb$biomass, ra$biomass)
  expect_equal(rb$total_yield, 2 * ra$total_yield)
  expect_equal(rb$total_yield_pred, 2 * ra$total_yield_pred)
  expect_equal(rb$W_catch, 2 * rb$W)
  expect_identical(tinyAM:::.population_unit(changed, "W_catch"), "kg")
})

test_that("plus-group catch weights use the same abundance shares as stock weights", {
  d <- make_test_dat(data = typed_weights(cod_obs), ages = 2:10)
  par <- make_par(d)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), par, silent = TRUE)
  r <- obj$report()
  shares <- r$N_plus / rowSums(r$N_plus)
  expect_equal(as.numeric(r$W[, ncol(r$W)]), unname(rowSums(shares * d$W_plus_input)))
  expect_equal(as.numeric(r$W_catch[, ncol(r$W)]), unname(rowSums(shares * d$W_catch_plus_input)))
  expect_equal(r$W_catch, 2 * r$W)
  shuffled <- typed_weights(cod_obs)
  shuffled$weight <- shuffled$weight[sample(nrow(shuffled$weight)), ]
  s <- make_test_dat(data = shuffled, ages = 2:10)
  expect_equal(s$W_catch, d$W_catch)
  expect_equal(s$W_catch_plus_input, d$W_catch_plus_input)
})

test_that("projection averaging and held-out catch targets use separate weights", {
  d <- typed_weights(cod_obs)
  d$weight$obs[d$weight$type == "catch"] <- d$weight$obs[d$weight$type == "catch"] *
    (d$weight$year[d$weight$type == "catch"] - min(d$weight$year) + 1)
  years <- 2010:2020
  prep <- make_test_dat(data = d, years = years,
    proj_settings = list(n_proj = 2, n_mean = 3, F_mult = 1))
  for (type in c("stock", "catch")) {
    w <- tinyAM:::.select_weights(d$weight, type)
    expected <- aggregate(obs ~ age, w[w$year %in% 2018:2020, ], mean)
    input <- if (type == "stock") prep$W else prep$W_catch
    expect_equal(as.numeric(input["2021", ]), expected$obs[match(AGES, expected$age)])
  }
  f <- default_fit; f$refit_args$data <- d
  ages <- sort(unique(c(f$dat$ages, f$dat$plus_ages)))
  targets <- tinyAM:::.hindcast_controls(f, 2021:2022, list(catch = "observed"))
  for (i in 1:2) {
    counts <- subset(d$catch, year == targets$year[i] & age %in% ages)
    w <- subset(d$weight, type == "catch" & year == targets$year[i] & age %in% ages)
    expect_equal(targets$target[i], sum(counts$obs * w$obs[match(counts$age, w$age)]))
  }
})

test_that("catch targets, simulation and refits retain catch weights", {
  d <- make_test_dat(data = typed_weights(cod_obs), years = 2000:2005, ages = 2:6,
    N_settings = list(process = "off"),
    proj_settings = list(n_proj = 2, n_mean = 2, catch = c(100, 80)))
  p <- make_par(d); p$log_r0 <- log(1e5); p$log_r[] <- log(1e5); p$log_f[] <- log(.15)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  expect_equal(unname(r$total_yield_pred[d$is_proj]), c(100, 80), tolerance = 1e-7)
  expect_true(all(is.finite(obj$gr())))
  set.seed(81)
  simulated <- nll_fun(p, d, simulate = TRUE)
  p$log_r[] <- simulated$log_r; p$log_f[] <- simulated$log_f
  simulated_report <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$report()
  expect_equal(unname(simulated_report$total_yield_pred[d$is_proj]), c(100, 80), tolerance = 1e-7)
  expect_equal(simulated_report$W_catch, 2 * simulated_report$W)
  f <- update(default_fit, data = typed_weights(cod_obs), start_par = default_fit, silent = TRUE)
  expect_equal(f$opt$objective, default_fit$opt$objective, tolerance = 1e-6)
  expect_equal(f$rep$N, default_fit$rep$N, tolerance = 1e-5)
  expect_true("type" %in% names(f$refit_args$data$weight))
  expect_true("W_catch" %in% names(tidy_pop(f)))
  sims <- sim_tam(f, n = 1, par_uncertainty = "none", redraw_random = FALSE)
  expect_true("W_catch" %in% names(sims))
  h <- fit_hindcast(f, folds = 1, start_from_fit = TRUE, progress = FALSE)
  successful <- Filter(function(x) inherits(x, "tam_fit"), h$fits[[1]])
  expect_gt(length(successful), 0)
  for (fit in successful) expect_equal(fit$dat$W_catch, 2 * fit$dat$W)
})
