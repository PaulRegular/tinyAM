# tests/testthat/test-check_obs.R

valid_obs <- function() {
  list(
    catch = data.frame(year = 2000:2002, age = rep(2:3, length.out = 3), obs = 1:3),
    index = data.frame(year = 2000:2002, age = rep(2:3, length.out = 3), obs = 1:3, samp_time = c(0.2, 0.8, NA)),
    weight = data.frame(year = 2000:2002, age = rep(2:3, length.out = 3), obs = c(0.5, 0.6, 0.7)),
    maturity = data.frame(year = 2000:2002, age = rep(2:3, length.out = 3), obs = c(0.2, 0.3, 0.4))
  )
}

test_that("check_obs passes on valid input and returns TRUE invisibly", {
  v <- valid_obs()
  expect_invisible(expect_true(check_obs(v)))
})

test_that("missing required tables aborts", {
  v <- valid_obs(); v$index <- NULL
  expect_error(check_obs(v), "Missing required tables", fixed = TRUE)
})

test_that("each table must have year/age/obs with numeric types", {
  v <- valid_obs(); v$catch$obs <- as.character(v$catch$obs)
  expect_error(check_obs(v), "`catch\\$obs` must be numeric", perl = TRUE)

  v <- valid_obs(); v$weight$year <- factor(v$weight$year)
  expect_error(check_obs(v), "`weight\\$year` must be numeric", perl = TRUE)
})

test_that("index samp_time must be numeric in [0,1] if present", {
  v <- valid_obs(); v$index$samp_time <- c(-0.1, 1.2, NA)
  expect_error(check_obs(v), "`index\\$samp_time` must be numeric in \\[0, 1\\]", perl = TRUE)

  v <- valid_obs(); v$index$samp_time <- c(0, 1, NA)  # boundary ok
  expect_invisible(expect_true(check_obs(v)))
})

test_that("duplicate (year, age) warns, not aborts", {
  v <- valid_obs()
  v$catch <- rbind(v$catch, v$catch[1, , drop = FALSE])
  expect_warning(check_obs(v), "duplicate \\(year, age\\) rows", perl = TRUE)
})


