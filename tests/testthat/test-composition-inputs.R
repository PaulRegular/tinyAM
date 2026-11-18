composition_input_fixture <- function() {
  grid <- expand.grid(age = 1:4, year = 2000:2003)
  index <- grid[grid$age > 1L, ]
  list(
    catch = transform(grid, obs = c(0, 0.2, 0.35, 0.45)[age]),
    index = transform(index, obs = c(0.21, 0.33, 0.464)[age - 1L],
      survey = "RV", samp_time = 0.8),
    weight = transform(grid, obs = age / 10),
    maturity = transform(grid, obs = c(0, 0.05, 0.25, 0.7)[age],
      unit = "proportion at age")
  )
}

test_that("explicit composition units are rejected for catch and index", {
  for (unit in c(" AGE_COMPOSITION ", "composition-at-age", "Proportion at age",
    "Proportions at age", "percent_at_age", "Percentage at age")) {
    catch <- composition_input_fixture()
    catch$catch$unit <- unit
    expect_error(check_tam_data(catch, quiet = TRUE), "catch")

    index <- composition_input_fixture()
    index$index$unit <- unit
    expect_error(check_tam_data(index, quiet = TRUE), "RV")
  }
  expect_error(as_tam_data(index), "composition")
  expect_error(fit_tam(index, years = 2000:2003, ages = 1:4), "composition")
})

test_that("repeated unit sums identify the relevant fleet and survey", {
  dat <- composition_input_fixture()
  dat$catch <- rbind(
    transform(dat$catch, fleet = "Gillnet"),
    transform(dat$catch, fleet = "Trawl", obs = obs / 2)
  )
  dat$index <- rbind(
    transform(dat$index, unit = "unknown"),
    transform(dat$index, survey = "Sentinel", unit = "fish per tow")
  )
  dat$index <- subset(dat$index, survey != "RV" | year != 2003 | age != 4)
  before <- dat

  warned <- character()
  result <- withCallingHandlers(suppressMessages(check_tam_data(dat)),
    warning = function(warning) {
      warned <<- c(warned, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  expect_length(warned, 1L)
  expect_match(warned, "age composition")
  expect_match(warned, "Gillnet")
  expect_match(warned, "RV")
  expect_false(grepl("Trawl|Sentinel", warned))
  expect_identical(result, TRUE)
  expect_identical(dat, before)
  expect_silent(check_tam_data(dat, quiet = TRUE))
})

test_that("unlabelled catch is checked as one series", {
  dat <- composition_input_fixture()
  dat$index <- dat$index[0, ]
  expect_warning(suppressMessages(check_tam_data(dat)), "catch")
})

test_that("known units and incomplete or weak evidence do not warn", {
  dat <- composition_input_fixture()
  dat$catch$unit <- "thousand fish"
  dat$index$unit <- "fish per tow"
  expect_no_warning(check_tam_data(dat))

  dat <- composition_input_fixture()
  dat$catch$obs[dat$catch$year == 2003] <- dat$catch$obs[dat$catch$year == 2003] / 10
  dat$index <- subset(dat$index, year != 2000 | age != 4)
  dat$index$obs[dat$index$year == 2001 & dat$index$age == 4] <- NA_real_
  expect_no_warning(check_tam_data(dat))

  dat$catch$obs <- as.numeric(dat$catch$age == 4)
  dat$index <- subset(dat$index, age == 3)
  dat$index$obs <- 1
  expect_no_warning(check_tam_data(dat))

  dat$index <- dat$index[0, ]
  expect_no_warning(check_tam_data(dat))
})
