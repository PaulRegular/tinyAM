test_that("unit metadata validates without changing observations", {
  d <- cod_obs
  d$catch$unit <- "thousand fish"
  d$weight$unit <- "kg/fish"
  d$index$unit <- "fish per tow"
  original <- d
  expect_true(check_tam_data(d))
  expect_identical(d, original)
  expect_equal(tinyAM:::.data_units(d)$mass, "tonnes")
  d$catch$unit[1] <- "fish"
  expect_error(check_tam_data(d), "inconsistent units")
  d <- original
  d$weight$unit[1] <- NA_character_
  expect_error(check_tam_data(d), "missing values")
  d <- original
  d$index$survey[1] <- "Other"
  d$index$unit[1] <- "fish"
  expect_true(check_tam_data(d))
})

test_that("units preserve designs, parameters, likelihoods and derivatives", {
  plain <- cod_obs
  for (nm in names(plain)) plain[[nm]]$unit <- NULL
  labeled <- plain
  labeled$catch$unit <- "thousand fish"
  labeled$weight$unit <- "kg/fish"
  labeled$index$unit <- "thousand fish"
  a <- prepare_tam(plain, years = 2000:2005, ages = 2:6)
  b <- prepare_tam(labeled, years = 2000:2005, ages = 2:6)
  p <- make_par(a)
  expect_identical(make_par(b), p)
  for (nm in c("F_modmat", "M_modmat", "q_modmat", "sd_catch_modmat", "sd_index_modmat")) expect_identical(a[[nm]], b[[nm]])
  oa <- RTMB::MakeADFun(function(p) nll_fun(p, a), p, silent = TRUE)
  ob <- RTMB::MakeADFun(function(p) nll_fun(p, b), p, silent = TRUE)
  expect_equal(oa$fn(), ob$fn())
  expect_equal(oa$gr(), ob$gr())
  expect_equal(oa$report(), ob$report())
  pop <- tidy_rep(list(dat = b, rep = ob$report()))
  expect_true(all(pop$N$unit == "thousand fish"))
  expect_true(all(pop$ssb$unit == "tonnes"))
  expect_true(all(pop$F$unit == "per year"))
  patterns <- tinyAM:::.tidy_patterns(pop)
  expect_equal(subset(patterns$population, metric == "Abundance")$est, pop$N$est * 1000)
  expect_equal(subset(patterns$population, metric == "SSB")$est, pop$ssb_mat$est * 1000)
  expect_identical(labeled$catch$obs, plain$catch$obs)
})

test_that("default plot labels include units and explicit labels win", {
  d <- data.frame(year = 2000:2002, est = 1:3, unit = "tonnes")
  expect_equal(plotly::plotly_build(plot_trend(d))$x$layout$yaxis$title, "Estimate (tonnes)")
  expect_equal(plotly::plotly_build(plot_trend(d, ylab = "My label"))$x$layout$yaxis$title, "My label")
  d$unit[1] <- "kg"
  expect_error(plot_trend(d), "incompatible units")
})
