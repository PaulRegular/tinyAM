terminal_pooling_data <- function() {
  grid <- expand.grid(year = 2000:2004, age = 2:6)
  grid$obs <- 5
  native_tail <- cbind(
    `4` = c(2, 0, 0, 7, NA), `5` = c(3, 0, 0, NA, NA), `6` = c(4, 0, NA, 0, NA)
  )
  tail <- grid$age >= 4L
  grid$obs[tail] <- native_tail[cbind(grid$year[tail] - 1999, grid$age[tail] - 3)]
  pooled <- grid
  pooled$obs[pooled$age == 4L] <- c(11, 0, 15, 12, NA)
  pooled$obs[pooled$age > 4L] <- NA_real_
  survey <- function(d, name) {
    transform(d, survey = name, samp_time = .25, label = paste0(name, "-", age))
  }
  list(
    catch = rbind(transform(grid, fleet = "Native"), transform(pooled, fleet = "Pooled"),
      transform(grid, fleet = "Legacy")),
    index = rbind(survey(subset(grid, !(year == 2002 & age == 5L)), "Native"),
      survey(subset(pooled, age <= 4L), "Pooled"), survey(grid, "Legacy")),
    weight = transform(grid, obs = age / 2),
    maturity = transform(grid, obs = as.numeric(age >= 4))
  )
}

test_that("native PLD terminal sums require every sampled age within each series", {
  raw <- terminal_pooling_data()
  original <- raw
  families <- list(Native = pld_gamma(10), Pooled = pld_gengamma(10, Q = 0))
  dat <- prepare_tam(raw, ages = 2:4,
    catch_settings = list(family = families, sd_form = ~ 0 + fleet),
    index_settings = list(family = families, sd_form = ~ 0 + survey, q_form = ~ 0 + survey))
  expect_identical(raw, original)
  expect_equal(dat$plus_ages, 4:6)
  for (component in c("catch", "index")) {
    d <- subset(dat$obs[[component]], age == 4L)
    series <- if (component == "catch") d$fleet else d$survey
    expect_equal(d$obs[series == "Native"], c(9, 0, NA, NA, NA))
    expect_equal(d$obs[series == "Pooled"], c(11, 0, 15, 12, NA))
    expect_equal(d$obs[series == "Legacy"], c(9, 0, 0, 7, NA))
  }
  d <- subset(dat$obs$index, age == 4L & survey == "Native")
  expect_equal(d$samp_time, rep(.25, 5L))
  expect_equal(d$label, rep("Native-4", 5L))
  map <- dat$obs_map
  zeros <- dat$pld$zero
  expect_true(any(map$age[zeros] == 4L & map$year[zeros] == 2001))
  expect_false(any(map$age[zeros] == 4L & map$year[zeros] == 2002))

  raw$index <- rbind(raw$index,
    subset(raw$index, year == 2002 & age == 4L & survey == "Native"))
  expect_error(prepare_tam(raw, ages = 2:4,
    index_settings = list(family = families, sd_form = ~ 0 + survey, q_form = ~ 0 + survey)),
    "Native PLD.*one row per year, age and series")
})

test_that("a training fold retains native PLD support when only later tail ages were sampled", {
  d <- expand.grid(year = 2000:2004, age = 2:6)
  d$obs <- rep(c(20, 15, 0, 2, 1), each = 5)
  d$obs[d$age > 4L & d$year <= 2002] <- NA_real_
  raw <- list(catch = d,
    index = transform(d, survey = "RV", samp_time = .25),
    weight = transform(d, obs = age),
    maturity = transform(d, obs = as.numeric(age >= 4)))
  args <- list(data = raw, ages = 2:4,
    index_settings = list(family = pld_lognormal(10), q_form = ~1),
    catch_settings = list(family = pld_lognormal(10)))
  full <- do.call(prepare_tam, args)
  fold <- do.call(prepare_tam, c(args, list(years = 2000:2002,
    proj_settings = list(n_proj = 1, n_mean = 1, F_mult = 1))))
  for (type in c("catch", "index")) {
    expect_equal(subset(full$obs[[type]], age == 4L)$obs, c(NA, NA, NA, 3, 3))
    expect_true(all(is.na(subset(fold$obs[[type]], age == 4L)$obs)))
  }
  expect_length(fold$pld$zero, 0L)

  reference <- list(dat = full, refit_args = list(data = raw))
  controls <- .hindcast_controls(reference, 2003:2004, .retro_projection(NULL))
  expected <- sum(subset(raw$catch, year == 2003)$obs * subset(raw$weight, year == 2003)$obs)
  expect_equal(controls$target, rep(expected, 2L))
  expect_true(all(is.na(controls$error)))
})

test_that("Northern cod native pooling matches the earlier 14+ inputs", {
  raw <- cod_obs
  prepooled <- raw
  for (type in c("catch", "index")) {
    d <- prepooled[[type]]
    sums <- vapply(split(subset(d, age %in% 14:20), subset(d, age %in% 14:20)$year),
      function(p) if (setequal(p$age, 14:20)) sum(p$obs) else NA_real_, numeric(1))
    d$obs[d$age == 14L] <- unname(sums[as.character(d$year[d$age == 14L])])
    d$obs[d$age > 14L] <- NA_real_
    prepooled[[type]] <- if (type == "index") subset(d, age <= 14L) else d
  }
  prep <- function(x) prepare_tam(x, years = 1983:2025, ages = 2:14,
    index_settings = list(q_form = ~ mono(q_block), family = pld_lognormal()),
    catch_settings = list(family = pld_lognormal()))
  native <- prep(raw)
  legacy <- prep(prepooled)
  expect_identical(native$obs, legacy$obs)
  expect_identical(native$obs_map, legacy$obs_map)
  expect_identical(native$log_obs, legacy$log_obs)
  expect_equal(sum(native$obs$catch$obs == 0, na.rm = TRUE), 57)
  expect_equal(sum(native$obs$index$obs == 0, na.rm = TRUE), 72)
  expect_equal(sum(!is.na(native$obs$catch$obs)), 559)
  expect_equal(sum(!is.na(native$obs$index$obs)), 520)
})
