
obs <- northern_cod_data

# ---- check_obs ----

test_that("check_obs passes on valid input and returns TRUE invisibly", {
  expect_invisible(expect_true(check_obs(obs)))
})

test_that("check_obs aborts when required tables are missing", {
  v <- obs
  v$index <- NULL
  expect_error(check_obs(v), "Missing required tables", fixed = TRUE)
})

test_that("check_obs enforces numeric types for year/age/obs", {
  v <- obs
  v$catch$obs <- as.character(v$catch$obs)
  expect_error(check_obs(v), "`catch\\$obs` must be numeric", perl = TRUE)

  v <- obs
  v$weight$year <- factor(v$weight$year)
  expect_error(check_obs(v), "`weight\\$year` must be numeric", perl = TRUE)

  v <- obs
  v$maturity$age <- as.factor(v$maturity$age)
  expect_error(check_obs(v), "`maturity\\$age` must be numeric", perl = TRUE)
})

test_that("check_obs requires index$samp_time present, numeric in [0,1], and no NAs", {
  # Out-of-range values
  v <- obs
  v$index$samp_time[] <- c(-0.1, 1.2)[(seq_len(nrow(v$index)) %% 2) + 1]
  expect_error(check_obs(v), "`index\\$samp_time` must be numeric in \\[0, 1\\]", perl = TRUE)

  # Missing column
  v <- obs
  v$index$samp_time <- NULL
  expect_error(check_obs(v), "index\\$samp_time", perl = TRUE)

  # NA not permitted
  v <- obs
  v$index$samp_time[1] <- NA_real_
  expect_error(check_obs(v), "index\\$samp_time", perl = TRUE)

  # Valid boundary values 0 and 1 pass
  v <- obs
  v$index$samp_time[] <- 0
  expect_invisible(expect_true(check_obs(v)))
  v$index$samp_time[] <- 1
  expect_invisible(expect_true(check_obs(v)))
})

test_that("check_obs warns on duplicate (year, age) rows", {
  v <- obs
  v$catch <- rbind(v$catch, v$catch[1, , drop = FALSE])
  expect_warning(check_obs(v), "duplicate \\(year, age\\) rows", perl = TRUE)
})

# ---- add_proj_rows ----

test_that("add_proj_rows appends projection rows and adds is_proj", {
  out <- add_proj_rows(obs, n_proj = 2, n_mean = 3, quiet = TRUE)

  for (nm in c("catch","index","weight","maturity")) {
    expect_true("is_proj" %in% names(out[[nm]]))
    maxy <- max(obs[[nm]]$year, na.rm = TRUE)
    expect_true(all((maxy + 1):(maxy + 2) %in% out[[nm]]$year))
    expect_true(any(out[[nm]]$is_proj))
  }
})

test_that("add_proj_rows sets index$obs to NA for projected rows", {
  out <- add_proj_rows(obs, n_proj = 2, n_mean = 3, quiet = TRUE)
  idx_proj <- subset(out$index, is_proj)
  expect_true(all(is.na(idx_proj$obs)))
})

test_that("add_proj_rows copies aux columns from terminal year (catch blocks)", {
  out  <- add_proj_rows(obs, n_proj = 1, n_mean = 3, quiet = TRUE)
  c_in <- obs$catch
  c_out <- out$catch
  maxy <- max(c_in$year, na.rm = TRUE)

  term <- unique(subset(c_in, year == maxy,
                        select = c(age, sd_obs_block, F_y_block, F_a_block)))
  proj <- unique(subset(c_out, year == maxy + 1 & is_proj,
                        select = c(age, sd_obs_block, F_y_block, F_a_block)))

  m <- merge(term, proj, by = "age", suffixes = c("_term","_proj"), sort = FALSE)
  expect_true(all(m$sd_obs_block_term == m$sd_obs_block_proj))
  expect_true(all(as.character(m$F_y_block_term) == as.character(m$F_y_block_proj)))
  expect_true(all(as.character(m$F_a_block_term) == as.character(m$F_a_block_proj)))
})

test_that("add_proj_rows uses average obs of last n_mean years for projections (catch)", {
  out <- add_proj_rows(obs, n_proj = 1, n_mean = 3, quiet = TRUE)

  max_year   <- max(obs$catch$year, na.rm = TRUE)
  mean_years <- (max_year - 3 + 1):max_year

  mean_catch <- aggregate(obs ~ age,
                          data = subset(obs$catch, year %in% mean_years),
                          FUN = function(z) mean(z, na.rm = TRUE))
  added_rows <- subset(out$catch, year == max_year + 1 & is_proj)

  # Compare by age after merging to ensure alignment
  cmp <- merge(mean_catch, added_rows[, c("age","obs")],
               by = "age", suffixes = c("_mean","_proj"), sort = FALSE)
  expect_equal(cmp$obs_proj, cmp$obs_mean, tolerance = 1e-12)
})

test_that("add_proj_rows aborts if is_proj already present", {
  v <- obs
  v$catch$is_proj <- FALSE
  expect_error(add_proj_rows(v, quiet = TRUE), regexp = "is_proj.*already exists", ignore.case = TRUE)
})

test_that("add_proj_rows emits an info message when quiet = FALSE", {
  expect_message(add_proj_rows(obs, quiet = FALSE), regexp = "Adding .* projection year", perl = TRUE)
})
