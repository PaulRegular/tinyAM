
## compute_mohns_rho ---

test_that("compute_mohns_rho: basic case matches hand calculation", {
  # terminal (fold = 2015) + two peels (2013, 2014)
  df <- data.frame(
    year = c(2013, 2014, 2015,
             2013, 2014, 2015,
             2013, 2014, 2015),
    fold = c(2013, 2013, 2013,   # peel 2013
             2014, 2014, 2014,   # peel 2014
             2015, 2015, 2015),  # terminal
    est  = c( 90,  NA,  NA,
              NA,  95,  NA,
              100, 100, 100),
    is_proj = FALSE
  )

  # Hand calc:
  # 2013 peel: (90 - 100)/100 = -0.10
  # 2014 peel: (95 - 100)/100 = -0.05
  # Mean = -0.075
  rho <- compute_mohns_rho(df)
  expect_equal(rho, -0.075, tolerance = 1e-12)
})

test_that("compute_mohns_rho: projection rows are ignored", {
  # Same as basic case but add distracting projection rows that should be dropped
  df <- data.frame(
    year = c(2013, 2014, 2015,
             2013, 2014, 2015,
             2013, 2014, 2015,
             # extra rows flagged as projections (should not affect rho)
             2016, 2016, 2016),
    fold = c(2013, 2013, 2013,
             2014, 2014, 2014,
             2015, 2015, 2015,
             2016, 2016, 2016),
    est  = c( 90,  NA,  NA,
              NA,  95,  NA,
              100, 100, 100,
              123,  50,  10),
    is_proj = c(rep(FALSE, 9), rep(TRUE, 3))
  )

  rho <- compute_mohns_rho(df)
  expect_equal(rho, -0.075, tolerance = 1e-12)
})

test_that("compute_mohns_rho: zero terminal estimates can yield Inf/NaN (documented)", {
  df <- data.frame(
    year = c(2015, 2015, 2015),
    fold = c(2014, 2015, 2016),  # terminal is max(fold) == 2016
    est  = c(5, 10, 0),          # terminal est_t at 2015 = 0
    is_proj = FALSE
  )
  rho <- compute_mohns_rho(df)
  expect_true(is.infinite(rho) || is.nan(rho))
})

test_that("compute_mohns_rho: returns NA when no valid peel/terminal pairs", {
  # Only terminal run (no peels) -> mean(numeric(0), na.rm=TRUE) = NA
  df <- data.frame(
    year    = 2010:2012,
    fold    = rep(2012, 3),
    est     = c(100, 110, 120),
    is_proj = FALSE
  )
  rho <- compute_mohns_rho(df)
  expect_true(is.na(rho))
})

test_that("compute_mohns_rho: missing values in est are ignored via na.rm=TRUE", {
  # Same as basic, but make one peel value NA at the used year — result should
  # be based on the remaining valid peel(s).
  df <- data.frame(
    year = c(2013, 2014, 2015,
             2013, 2014, 2015,
             2013, 2014, 2015),
    fold = c(2013, 2013, 2013,
             2014, 2014, 2014,
             2015, 2015, 2015),
    est  = c( NA,  NA,  NA,     # peel 2013 now NA at 2013 (drops from mean)
              NA,  95,  NA,     # peel 2014 valid at 2014
              100, 100, 100),    # terminal
    is_proj = FALSE
  )
  # Only the 2014 peel contributes: (95-100)/100 = -0.05
  rho <- compute_mohns_rho(df)
  expect_equal(rho, -0.05, tolerance = 1e-12)
})


## compute_hindcast_rmse ---

test_that("compute_hindcast_rmse: basic case matches hand calculation (natural scale)", {
  # Two ages, one fold; obs at year==fold; projections flagged by is_proj==TRUE
  d <- data.frame(
    year    = c(2019, 2019, 2019, 2019),
    age     = c(2,    3,    2,    3),
    obs     = c(10,   20,   NA,   NA),
    pred    = c(NA,   NA,   11,   18),
    is_proj = c(FALSE,FALSE,TRUE, TRUE),
    fold    = c(2019, 2019, 2019, 2019)
  )

  # Hand calc: errors = (10-11)^2 = 1 ; (20-18)^2 = 4 ; mean = 2.5 ; RMSE = sqrt(2.5)
  rmse <- compute_hindcast_rmse(d, log = FALSE)
  expect_equal(rmse, sqrt(2.5), tolerance = 1e-12)
})

test_that("compute_hindcast_rmse: log-scale drops exact zeros (documented behavior)", {
  # One pair has a zero (dropped before log), the other is positive
  d <- data.frame(
    year    = c(2019, 2019, 2019, 2019),
    age     = c(2,    3,    2,    3),
    obs     = c(0,    10,   NA,   NA),   # zero observed at age 2 -> dropped when log=TRUE
    pred    = c(NA,   NA,   5,    20),
    is_proj = c(FALSE,FALSE,TRUE, TRUE),
    fold    = c(2019, 2019, 2019, 2019)
  )

  # Remaining pair (age 3): RMSE on log scale = |log(10) - log(20)| = log(2)
  rmse <- compute_hindcast_rmse(d, log = TRUE)
  expect_equal(rmse, log(2), tolerance = 1e-12)
})

test_that("compute_hindcast_rmse: returns NA when no valid matches exist", {
  # No projected rows -> merge is empty -> mean(numeric(0), na.rm=TRUE) -> NA
  d <- data.frame(
    year    = c(2019, 2019),
    age     = c(2,    3),
    obs     = c(10,   20),
    pred    = c(NA,   NA),
    is_proj = c(FALSE,FALSE),
    fold    = c(2019, 2019)
  )
  rmse <- compute_hindcast_rmse(d, log = FALSE)
  expect_true(is.na(rmse))
})

test_that("compute_hindcast_rmse: aggregates across multiple ages correctly", {
  # Three ages, same fold; natural scale
  d <- data.frame(
    year    = c(2019, 2019, 2019,   2019, 2019, 2019),
    age     = c(2,    3,    4,      2,    3,    4),
    obs     = c(10,   20,   30,     NA,   NA,   NA),
    pred    = c(NA,   NA,   NA,     12,   18,   27),
    is_proj = c(FALSE,FALSE,FALSE,  TRUE, TRUE, TRUE),
    fold    = c(2019, 2019, 2019,   2019, 2019, 2019)
  )
  # SEs: (10-12)^2=4, (20-18)^2=4, (30-27)^2=9 -> mean=17/3 -> RMSE=sqrt(17/3)
  rmse <- compute_hindcast_rmse(d, log = FALSE)
  expect_equal(rmse, sqrt(17/3), tolerance = 1e-12)
})


