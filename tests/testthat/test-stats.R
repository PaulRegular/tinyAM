
test_that("compute_mohns_rho: basic case matches hand calculation", {
  # terminal (fold = 2015) + two peels (2013, 2014)
  df <- data.frame(
    year = c(2013, 2014, 2015,
             2013, 2014, 2015,
             2013, 2014, 2015),
    fold = c(2013, 2013, 2013,  # peel 2013
                   2014, 2014, 2014,  # peel 2014
                   2015, 2015, 2015), # terminal
    est = c( 90,  NA,  NA,   # peel 2013: only value used is at year=2013
             NA,  95,  NA,   # peel 2014: only value used is at year=2014
             100, 100, 100)  # terminal: values at all years
  )

  # For 2013 peel: (90 - 100)/100 = -0.10
  # For 2014 peel: (95 - 100)/100 = -0.05
  # Mean = (-0.10 + -0.05)/2 = -0.075
  rho <- compute_mohns_rho(df)
  expect_equal(rho, -0.075, tolerance = 1e-12)
})

test_that("compute_mohns_rho: stepped-forward convention (+1 year shift) handled", {
  # Here fold labels are one year behind the estimate year used.
  # Function will add +1 to fold if max(year) > max(fold).
  df <- data.frame(
    year = c(2013, 2014, 2015,
             2013, 2014, 2015,
             2013, 2014, 2015),
    fold = c(2012, 2012, 2012,  # peel labeled as 2012 -> shifted to 2013
                   2013, 2013, 2013,  # peel 2013 -> shifted to 2014
                   2014, 2014, 2014), # terminal labeled as 2014 -> shifted to 2015
    est = c( 90,  NA,  NA,   # peel that will match year 2013 after shift
             NA,  95,  NA,   # peel that will match year 2014 after shift
             100, 100, 100)  # terminal series (after shift -> 2015 terminal)
  )
  rho <- compute_mohns_rho(df)
  # Same numbers as the basic case, just shifted; expect same result
  expect_equal(rho, -0.075, tolerance = 1e-12)
})

test_that("compute_mohns_rho: zero terminal estimates produce Inf/NaN (documented behavior)", {
  df2 <- data.frame(
    year = c(2015, 2015, 2015),
    fold = c(2014, 2015, 2016),  # terminal is fold == 2016, est_t at 2015 = 0
    est = c(5, 10, 0)                  # peel 2014: est_r at year 2015 = 5; terminal 2016: est_t at 2015 = 0
  )
  rho <- compute_mohns_rho(df2)
  expect_true(is.infinite(rho) || is.nan(rho))
})

test_that("compute_mohns_rho: returns NA when no valid peel/terminal pairs", {
  # Only terminal run present; no peels -> mean(numeric(0)) is NA
  df <- data.frame(
    year = 2010:2012,
    fold = rep(2012, 3),
    est = c(100, 110, 120)
  )
  rho <- compute_mohns_rho(df)
  expect_true(is.na(rho))
})

test_that("compute_mohns_rho: fold coercion works (factor/character)", {
  df_char <- data.frame(
    year = c(2013, 2014, 2013, 2014),
    fold = c("2013", "2013", "2014", "2014"),
    est = c(90, NA, 100, 100)  # terminal fold = "2014"
  )
  expect_silent(rho_char <- compute_mohns_rho(df_char))
  expect_true(is.finite(rho_char))

  df_fac <- df_char
  df_fac$fold <- factor(df_fac$fold)
  expect_silent(rho_fac <- compute_mohns_rho(df_fac))
  expect_equal(rho_char, rho_fac)
})

