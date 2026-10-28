power_data <- function() {
  d <- cod_obs
  d$index$power_group <- ifelse(d$index$age <= 3, paste0("age", d$index$age), NA_character_)
  d
}

test_that("power parses categorical groups without dropping inactive observations", {
  d <- make_test_dat(data = power_data(), years = 2000:2005,
    index_settings = list(q_form = ~q_block + power(N, by = power_group)))
  p <- make_par(d)
  expect_setequal(names(p$log_q_power), c("age2", "age3"))
  expect_equal(length(d$q_power$index), nrow(d$obs$index))
  expect_equal(d$q_power$index == 0L, is.na(d$obs$index$power_group))
  expect_false(any(grepl("sd.*power|eta.*power", names(p))))
  all <- make_test_dat(index_settings = list(q_form = ~q_block + power(N)))
  expect_identical(all$q_power$levels, "all")
  z <- power_data(); z$index$power_group <- match(z$index$power_group, c("age2", "age3"))
  numeric <- make_test_dat(data = z,
    index_settings = list(q_form = ~q_block + power(N, by = power_group)))
  expect_equal(length(make_par(numeric)$log_q_power), 2L)
  z$index$power_group[] <- NA
  inactive <- make_test_dat(data = z,
    index_settings = list(q_form = ~q_block + power(N, by = power_group)))
  expect_false("log_q_power" %in% names(make_par(inactive)))
})

test_that("power uses latent abundance at survey time and positive exponents", {
  d <- make_test_dat(data = power_data(), years = 2000:2005,
    index_settings = list(q_form = ~q_block + power(N, by = power_group)))
  p <- make_par(d)
  p$log_q_power[] <- log(c(.7, 1.4))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  r <- obj$report()
  ii <- d$obs_map$type == "index"
  cells <- cbind(as.character(d$obs$index$year), as.character(d$obs$index$age))
  log_a <- log(r$N[cells]) - r$Z[cells] * d$obs$index$samp_time
  expected <- exp(c(0, p$log_q_power))[d$q_power$index + 1L]
  expect_equal(as.numeric(r$q_power), as.numeric(expected))
  expect_equal(as.numeric(r$log_pred[ii]), as.numeric(r$log_q_base + expected * log_a))
  expect_equal(as.numeric(r$log_q_obs), as.numeric(r$log_q_base + (expected - 1) * log_a))
  expect_true(all(is.finite(obj$gr())))
  expect_true(all(is.finite(obj$he())))
  simulated <- nll_fun(p, d, simulate = TRUE)
  next_p <- p; names_to_copy <- intersect(names(p), names(simulated))
  next_p[names_to_copy] <- simulated[names_to_copy]
  next_r <- RTMB::MakeADFun(function(p) nll_fun(p, d), next_p, silent = TRUE)$report()
  expect_equal(as.numeric(next_r$q_power), as.numeric(expected))
  expect_false(isTRUE(all.equal(r$log_q_obs, next_r$log_q_obs)))
})

test_that("p equal to one preserves the original likelihood and other derivatives", {
  plain <- make_test_dat(years = 2000:2005)
  powered <- make_test_dat(years = 2000:2005,
    index_settings = list(q_form = ~q_block + power(N)))
  a <- RTMB::MakeADFun(function(p) nll_fun(p, plain), make_par(plain), silent = TRUE)
  b <- RTMB::MakeADFun(function(p) nll_fun(p, powered), make_par(powered), silent = TRUE)
  expect_equal(a$fn(), b$fn(), tolerance = 1e-10)
  expect_equal(a$gr(), b$gr()[names(b$par) != "log_q_power"], ignore_attr = TRUE, tolerance = 1e-9)
  expect_equal(a$report()$log_pred, b$report()$log_pred, tolerance = 1e-10)
})

test_that("unsupported power specifications fail explicitly", {
  expect_error(power(N), "formula marker")
  for (form in list(~power(obs), ~power(log(N)), ~power(N):age,
    ~I(power(N)), ~power(N) + power(N, by = survey))) {
    expect_error(make_test_dat(index_settings = list(q_form = form)), "power|Structured")
  }
  expect_error(make_test_dat(index_settings = list(q_form = ~power(N), q_link = "logit")), "q_link")
  expect_error(make_test_dat(catch_settings = list(sd_form = ~power(N))), "only supported")
  expect_error(make_test_dat(rec_settings = list(mu_form = ~power(N))), "only supported")
  expect_error(make_test_dat(F_settings = list(mu_form = ~power(N))), "only supported")
  z <- power_data(); z$index$power_group[z$index$obs > 0] <- NA
  expect_error(make_test_dat(data = z,
    index_settings = list(q_form = ~q_block + power(N, by = power_group))), "positive observations")
})

test_that("power fits report effective q with uncertainty and align starts", {
  f <- update(default_fit, index_settings = list(q_form = ~q_block + power(N)),
    rec_settings = list(init = "latent"), start_par = default_fit, silent = TRUE)
  z <- f$obs_pred$index
  expect_true(all(c("q", "q_base", "q_power", "q_base_lwr", "q_base_se_scale") %in% names(z)))
  expect_equal(z$pred, z$q * f$rep$N[cbind(as.character(z$year), as.character(z$age))] *
    exp(-f$rep$Z[cbind(as.character(z$year), as.character(z$age))] * z$samp_time))
  expect_true(all(z$q_power > 0))
  expect_false(any(c("q_base", "q_power") %in% names(f$pop)))
  expect_true(any(is.finite(z$q_lwr)))
  expect_true(any(is.finite(z$q_base_lwr)))
  expect_true(all(is.na(z$q_base_unit)))
  expect_identical(unique(z$q_base_se_scale), "log")
  expect_true("q_power" %in% f$fixed_par$par)
  expect_false("log_q_power" %in% f$obj$env$.random)
  expect_true("log_r0" %in% f$obj$env$.random)
  expect_equal(tinyAM:::.start_parameters(f, f$dat)$log_q_power,
    as.list(f$sdrep, "Estimate")$log_q_power)
  switched <- f$dat; switched$q_power$levels <- "different"
  expect_length(tinyAM:::.start_parameters(f, switched)$log_q_power, 0L)
  projected <- update(f, proj_settings = list(n_proj = 2, F_mult = 1), start_par = f, silent = TRUE)
  pz <- subset(projected$obs_pred$index, is_proj)
  expect_true(nrow(pz) > 0)
  expect_true(all(is.finite(pz$q) & pz$q_power > 0))
  r <- residuals(f, type = "standardized")
  expect_equal(sum(is.finite(r$residual)), sum(f$dat$is_observed))
  retro <- fit_hindcast(f, folds = 1, start_from_fit = TRUE, progress = FALSE,
    proj_settings = list(F_mult = 1))
  expect_equal(nrow(retro$fold_diagnostics), 2L)
  expect_identical(f$dat$index_settings$q_form, ~q_block + power(N))
  before <- f$obj$env$last.par
  testthat::local_mocked_bindings(oneStepPredict = function(obj, ...) {
    expect_equal(as.numeric(obj$fn(f$opt$par)), f$opt$objective, tolerance = 1e-6)
    list(residual = rep(0, sum(f$dat$is_observed)))
  }, .package = "RTMB")
  osa <- residuals(f, type = "osa", trace = FALSE)
  expect_equal(sum(is.finite(osa$residual)), sum(f$dat$is_observed))
  expect_identical(f$obj$env$last.par, before)
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available() &&
      requireNamespace("flexdashboard", quietly = TRUE)) {
    file <- tempfile(fileext = ".html")
    vis_tam(f, output_file = file, open_file = FALSE, render_args = list(quiet = TRUE))
    expect_true(file.exists(file))
    con <- file(file, "rb"); html <- paste(readLines(con, warn = FALSE), collapse = "\n"); close(con)
    expect_match(html, "<h3>Initial recruitment</h3>", fixed = TRUE)
    expect_match(html, "q_power", fixed = TRUE)
  }
})
