pkgload::load_all(quiet = TRUE)
source("analysis/formula_expansion/recruitment_validation.R")

testthat::test_that("simulated observations retain the prepared survey/year/age row mapping", {
  d <- do.call(prepare_tam, c(list(data = rec_inputs()), rec_settings(~ rw(year))))
  generated <- list(log_obs = log(seq_len(nrow(d$obs_map))))
  obs <- rec_observations(generated, d)
  again <- do.call(prepare_tam, c(list(data = obs), rec_settings(~ rw(year))))
  testthat::expect_equal(again$log_obs, generated$log_obs)
  testthat::expect_equal(again$obs_map[c("type", "year", "age", "survey")],
                         d$obs_map[c("type", "year", "age", "survey")])
})

testthat::test_that("stock-recruit derivatives with respect to parent SSB are exact", {
  S <- c(1, 30, 100)
  p <- list(log_sr_alpha = log(4), log_sr_beta = log(.02))
  for (type in c("bh", "ricker")) {
    obj <- RTMB::MakeADFun(function(s) sum(exp(tinyAM:::.rec_log_curve(log(s$ssb), p, type))),
                          list(ssb = S), silent = TRUE)
    expected <- if (type == "bh") 4 / (1 + .02 * S)^2 else 4 * exp(-.02 * S) * (1 - .02 * S)
    testthat::expect_equal(unname(obj$gr(obj$par)), matrix(expected, 1L), tolerance = 1e-10)
  }
})

testthat::test_that("the isolated study uses the exact recruitment likelihood", {
  for (case in rec_cases) {
    result <- rec_isolated_recovery(case, "well informed", 123L)
    testthat::expect_true(is.finite(result$diagnostics$objective))
    testthat::expect_true(all(is.finite(result$curve$estimate)))
    testthat::expect_lt(result$curve_rmse, .2)
    testthat::expect_equal(result$parameters$truth[1L], log(.35))
  }
})

testthat::test_that("simulated stock-recruit curves use the current biological surface", {
  for (case in rec_cases[1:4]) {
    d <- do.call(prepare_tam, c(list(data = rec_inputs()), rec_settings(rec_formula(case))))
    p <- make_par(d)
    set.seed(952L)
    simulated <- nll_fun(p, d, simulate = TRUE)
    shared <- intersect(names(p), names(simulated))
    p[shared] <- simulated[shared]
    report <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$report()
    expected <- tinyAM:::.rec_log_curve(log(report$ssb[d$rec$eligible - 1L]), p, d$rec$curve$type)
    testthat::expect_equal(unname(report$rec_log_mean), unname(expected))
    testthat::expect_equal(unname(report$rec_residual), unname(log(report$recruitment[d$rec$eligible]) - expected))
  }
})
