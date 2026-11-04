test_that("input validation reports success while preserving its logical return", {
  expect_message(expect_invisible(check_tam_data(cod_obs)), "Input data passed validation")
  expect_message(check_tam_data(cod_obs), "fleets: 1; surveys: 1")
  expect_silent(result <- check_tam_data(cod_obs, quiet = TRUE))
  expect_identical(result, TRUE)
  for (quiet in list(NA, 1, "yes", logical())) {
    expect_error(check_tam_data(cod_obs, quiet = quiet), "quiet.*TRUE or FALSE")
  }
  bad <- cod_obs
  bad$catch$obs[1] <- -1
  expect_error(check_tam_data(bad, quiet = TRUE), "non-negative")
})

test_that("validation summaries count fleet and survey identities", {
  d <- cod_obs
  d$catch <- rbind(transform(d$catch, fleet = "Trawl"), transform(d$catch, fleet = "Gillnet"))
  d$index <- rbind(transform(d$index, survey = "RV"), transform(d$index, survey = "Sentinel"))
  expect_message(check_tam_data(d), "fleets: 2; surveys: 2")
  expect_silent(as_tam_data(d))
  expect_silent(make_tam_data(d$catch, d$index, d$weight, d$maturity))
})

test_that("dashboard rendering checks Pandoc before writing output", {
  local_mocked_bindings(pandoc_available = function(...) FALSE, .package = "rmarkdown")
  local_mocked_bindings(render = function(...) stop("Must not render"), .package = "rmarkdown")
  file <- tempfile(fileext = ".html")
  expect_error(vis_tam(default_fit, output_file = file, open_file = FALSE),
    "Pandoc is required")
  expect_false(file.exists(file))
})
