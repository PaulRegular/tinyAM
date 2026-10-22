test_that("fit starts use stored modes without requiring uncertainty", {
  dat <- prepare_tam(cod_obs, years = 1983:1990, ages = 2:5)
  p <- make_par(dat)
  source <- structure(list(parameter_values = 123,
    obj = list(env = list(parList = function(par) { expect_equal(par, 123); p })),
    sdrep = NULL, dat = dat), class = "tam_fit")
  before <- source
  expect_equal(.start_parameters(source, dat), p)
  expect_identical(source, before)
  expect_equal(.merge_start_par(p, .start_parameters(source, dat)), p)
  changed <- dat; changed$index_settings$q_link <- "logit"
  expect_message(starts <- .start_parameters(source, changed), "q_link changed")
  expect_false("log_q" %in% names(starts))
  expect_equal(starts$log_f, p$log_f)
  expect_error(.start_parameters(1, dat), "named parameter list")
})
