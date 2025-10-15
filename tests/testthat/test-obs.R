
# ---- check_obs ----

test_that("check_obs passes on valid input and returns TRUE invisibly", {
  expect_invisible(expect_true(check_obs(cod_obs)))
})

test_that("check_obs aborts when required tables are missing", {
  obs <- cod_obs
  obs$index <- NULL
  expect_error(check_obs(obs), "Missing required tables", fixed = TRUE)
})

test_that("check_obs enforces numeric types for year/age/obs", {
  obs <- cod_obs
  obs$catch$obs <- as.character(obs$catch$obs)
  expect_error(check_obs(obs), "`catch\\$obs` must be numeric", perl = TRUE)

  obs <- cod_obs
  obs$weight$year <- factor(obs$weight$year)
  expect_error(check_obs(obs), "`weight\\$year` must be numeric", perl = TRUE)

  obs <- cod_obs
  obs$maturity$age <- as.factor(obs$maturity$age)
  expect_error(check_obs(obs), "`maturity\\$age` must be numeric", perl = TRUE)
})

test_that("check_obs enforces integer-ish years (no fractional years)", {
  obs <- cod_obs
  obs$catch$year <- as.numeric(obs$catch$year) + 0.5
  expect_error(check_obs(obs), "`catch\\$year` must be numeric", ignore.case = TRUE)
})

test_that("check_obs requires index$survey present", {
  obs <- cod_obs
  obs$index$survey <- NULL
  expect_error(check_obs(obs), "index\\$survey", perl = TRUE)
})

test_that("check_obs requires index$samp_time present, numeric in [0,1], and no NAs", {
  # Out-of-range values
  obs <- cod_obs
  obs$index$samp_time[] <- c(-0.1, 1.2)[(seq_len(nrow(obs$index)) %% 2) + 1]
  expect_error(check_obs(obs), "`index\\$samp_time` must be numeric in \\[0, 1\\]", perl = TRUE)

  # Missing column
  obs <- cod_obs
  obs$index$samp_time <- NULL
  expect_error(check_obs(obs), "index\\$samp_time", perl = TRUE)

  # NA not permitted
  obs <- cod_obs
  obs$index$samp_time[1] <- NA_real_
  expect_error(check_obs(obs), "index\\$samp_time", perl = TRUE)

  # Valid boundary values 0 and 1 pass
  obs <- cod_obs
  obs$index$samp_time[] <- 0
  expect_invisible(expect_true(check_obs(obs)))
  obs$index$samp_time[] <- 1
  expect_invisible(expect_true(check_obs(obs)))
})

test_that("check_obs checks that weight and maturity have no NA in obs; catch may have NA", {
  obs <- cod_obs
  obs$weight$obs[1] <- NA_real_
  expect_error(check_obs(obs), "weight\\$obs.*no NA", perl = TRUE)

  obs <- cod_obs
  obs$maturity$obs[1] <- NA_real_
  expect_error(check_obs(obs), "maturity\\$obs.*no NA", perl = TRUE)

  # catch may include NA in obs
  obs <- cod_obs
  obs$catch$obs[1] <- NA_real_
  expect_invisible(expect_true(check_obs(obs)))
})

test_that("check_obs checks that catch/weight/maturity covers full modeled (year, age) grid; index can be partial", {
  # Drop a few rows from weight -> should fail coverage check
  obs <- cod_obs
  drop_n <- min(3L, nrow(obs$weight))
  obs$weight <- obs$weight[-seq_len(drop_n), , drop = FALSE]
  expect_error(check_obs(obs), "`weight` coverage is incomplete", fixed = TRUE)

  obs <- cod_obs
  drop_n <- min(3L, nrow(obs$index))
  obs$index <- obs$index[-seq_len(drop_n), , drop = FALSE]
  expect_invisible(expect_true(check_obs(obs)))
})

