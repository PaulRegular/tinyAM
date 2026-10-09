pkgload::load_all(quiet = TRUE)
source(file.path("analysis", "formula_expansion", "q_validation.R"))

testthat::test_that("recovery designs use native q curves for both links and all terms", {
  for (link in c("log", "logit")) for (spec in q_specs) {
    dat <- do.call(prepare_tam, c(list(data = q_full_observations()), q_full_settings(spec, link)))
    p <- q_truth_parameters(make_par(dat), dat, link)
    term <- dat$q_terms[[1]]
    if (term$type != "logistic") p[[term$parameter]][] <- .15
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, dat), p, silent = TRUE)
    testthat::expect_equal(as.numeric(obj$report()$log_q_obs),
      as.numeric(q_log_curve(p, dat)), tolerance = 1e-10)
    testthat::expect_true(is.finite(obj$fn(obj$par)))
    testthat::expect_true(all(is.finite(obj$gr(obj$par))))
    if (link == "logit") testthat::expect_true(all(exp(q_log_curve(p, dat)) < 1))
    if (term$type != "logistic") {
      testthat::expect_identical(names(p[[term$parameter]]),
        unlist(lapply(term$groups, `[[`, "states"), use.names = FALSE))
    }
  }
})

testthat::test_that("sparse log-link IID designs are rejected before fitting", {
  testthat::expect_error(do.call(prepare_tam,
    c(list(data = q_full_observations(TRUE)), q_full_settings("iid", "log"))), "observation error")
})

testthat::test_that("uncertainty is transformed consistently with optimized parameters", {
  truth <- list(log_sd_q_iid_year = log(.18), logit_phi_q_ar1_year = qlogis(.65),
    q_a50_logistic_age = c(3.5, 4.5), log_q_slope_logistic_age = log(c(1.2, 1)))
  errors <- lapply(truth, function(x) rep(.1, length(x)))
  tab <- q_parameter_intervals(truth, errors, truth)
  testthat::expect_equal(tab$estimate, tab$truth)
  testthat::expect_true(all(tab$lower < tab$truth & tab$upper > tab$truth))
  testthat::expect_equal(tab$lower[1], .18 * exp(-qnorm(.975) * .1))
  testthat::expect_equal(tab$upper[2], plogis(qlogis(.65) + qnorm(.975) * .1))
})
