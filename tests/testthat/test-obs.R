
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

test_that("check_obs aborts when required tables are missing", {
  v <- valid_obs(); v$index <- NULL
  expect_error(check_obs(v), "Missing required tables", fixed = TRUE)
})

test_that("check_obs checks that each table must have year/age/obs with numeric types", {
  v <- valid_obs(); v$catch$obs <- as.character(v$catch$obs)
  expect_error(check_obs(v), "`catch\\$obs` must be numeric", perl = TRUE)

  v <- valid_obs(); v$weight$year <- factor(v$weight$year)
  expect_error(check_obs(v), "`weight\\$year` must be numeric", perl = TRUE)
})

test_that("check_obs checks that index samp_time is present and numeric in [0,1]", {
  v <- valid_obs(); v$index$samp_time <- c(-0.1, 1.2, NA)
  expect_error(check_obs(v), "`index\\$samp_time` must be numeric in \\[0, 1\\]", perl = TRUE)

  v <- valid_obs(); v$index$samp_time <- c(0, 1, NA)
  expect_error(check_obs(v))

  v <- valid_obs(); v$index$samp_time <- NULL
  expect_error(check_obs(v))
})

test_that("check_obs warns if duplicate (year, age) rows are present", {
  v <- valid_obs()
  v$catch <- rbind(v$catch, v$catch[1, , drop = FALSE])
  expect_warning(check_obs(v), "duplicate \\(year, age\\) rows", perl = TRUE)
})




test_that("add_proj_rows appends projection rows and adds is_proj", {
  obs <- northern_cod_data

  out <- add_proj_rows(obs, n_proj = 2, n_mean = 3, quiet = TRUE)

  for (nm in c("catch","index","weight","maturity")) {
    expect_true("is_proj" %in% names(out[[nm]]))
    maxy <- max(obs[[nm]]$year, na.rm = TRUE)
    expect_true(all((maxy + 1):(maxy + 2) %in% out[[nm]]$year))
    expect_true(any(out[[nm]]$is_proj))
  }
})

test_that("add_proj_rows sets index$obs to NA", {
  out <- add_proj_rows(northern_cod_data, n_proj = 2, n_mean = 3, quiet = TRUE)

  idx_proj <- subset(out$index, is_proj)
  expect_true(all(is.na(idx_proj$obs)))
})

test_that("add_proj_rows copies aux columns from terminal year", {
  out <- add_proj_rows(northern_cod_data, n_proj = 1, n_mean = 3, quiet = TRUE)

  c_in  <- northern_cod_data$catch
  c_out <- out$catch
  maxy  <- max(c_in$year, na.rm = TRUE)

  term <- unique(subset(c_in, year == maxy,
                        select = c(age, sd_obs_block, F_y_block, F_a_block)))
  proj <- unique(subset(c_out, year == maxy + 1 & is_proj,
                        select = c(age, sd_obs_block, F_y_block, F_a_block)))

  m <- merge(term, proj, by = "age", suffixes = c("_term","_proj"), sort = FALSE)
  expect_true(all(m$sd_obs_block_term == m$sd_obs_block_proj))
  expect_true(all(as.character(m$F_y_block_term) == as.character(m$F_y_block_proj)))
  expect_true(all(as.character(m$F_a_block_term) == as.character(m$F_a_block_proj)))
})

test_that("add_proj_rows uses average obs for is_proj rows", {
  obs <- northern_cod_data
  out <- add_proj_rows(obs, n_proj = 1, n_mean = 3, quiet = TRUE)

  max_year <- max(obs$catch$year, na.rm = TRUE)
  mean_years <- (max_year - 3 + 1):max_year

  mean_catch <- aggregate(obs ~ age,
                          data = obs$catch,
                          FUN = mean,
                          subset = year %in% mean_years)
  added_rows <- subset(out$catch, year == max_year + 1 & is_proj)

  expect_equal(added_rows$obs, mean_catch$obs, tolerance = 1e-12)
})

test_that("add_proj_rows aborts if is_proj already present", {
  obs <- northern_cod_data
  obs$catch$is_proj <- FALSE
  expect_error(add_proj_rows(obs, quiet = TRUE), regexp = "is_proj.*already exists", ignore.case = TRUE)
})

test_that("add_proj_rows emits an info message when quiet = FALSE", {
  expect_message(add_proj_rows(northern_cod_data, quiet = FALSE), regexp = "Adding .* projection year", perl = TRUE)
})

