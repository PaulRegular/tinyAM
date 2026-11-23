test_that("MASE uses identical cases and handles zero baseline error", {
  d <- data.frame(obs = c(4, 8, 9), pred = c(3, 6, 7), baseline = c(2, 4, NA))
  expect_equal(compute_hindcast_mase(d, FALSE), 3 / 6)
  expect_equal(compute_hindcast_mase(d), sum(abs(log(c(4, 8)) - log(c(3, 6)))) / (2 * log(2)))
  d$baseline <- d$obs
  expect_true(is.na(compute_hindcast_mase(d)))
  expect_match(.score_forecasts(d)$mase_reason, "Zero persistence")
  expect_error(compute_hindcast_mase(data.frame()), "baseline")
})

test_that("persistence never looks beyond the forecast origin", {
  actual <- data.frame(year = c(2000, 2002, 2003, 2004), age = 2, survey = "S",
    obs = c(10, 20, 0, 100), is_proj = FALSE)
  pred <- data.frame(year = c(2004, 2004), age = 2, survey = "S",
    pred = c(30, 40), is_proj = TRUE, fold = c(2002, 2003), horizon = c(2, 1))
  d <- .forecast_pairs(pred, actual)
  expect_equal(d$baseline, c(20, 20))
  expect_equal(d$baseline_year, c(2002, 2002))
  expect_equal(d$baseline_gap, c(2, 2))
  expect_equal(nrow(d), 2L)
  expect_equal(compute_hindcast_rmse(d), sqrt(mean((log(d$obs) - log(d$pred))^2)))
})

test_that("survey aggregates sum index predictions and require complete ages", {
  actual <- expand.grid(year = 2000:2002, age = 1:2)
  actual$survey <- "S"; actual$obs <- actual$year - 1990 + actual$age; actual$is_proj <- FALSE
  p <- subset(actual, year == 2002); p$pred <- p$obs * 2
  p$is_proj <- TRUE; p$fold <- 2001; p$horizon <- 1
  result <- .index_score_tables(.forecast_pairs(p, actual), actual)
  expect_equal(result$aggregate$pred, sum(p$pred))
  expect_equal(result$aggregate$obs, sum(p$obs))
  expect_equal(result$aggregate$baseline, sum(actual$obs[actual$year == 2001]))
  p <- p[1, ]
  result <- .index_score_tables(.forecast_pairs(p, actual), actual)
  expect_equal(nrow(result$aggregate), 0L)
  expect_true(any(result$scores$level == "age"))
})

test_that("composition scores retain zeros and exclude incomplete vectors", {
  actual <- expand.grid(year = 2000:2002, age = 1:3)
  actual$obs <- rep(c(0, 1, 2), each = 3); actual$is_proj <- FALSE
  p <- subset(actual, year == 2002); p$pred <- c(1, 1, 2)
  p$is_proj <- TRUE; p$fold <- 2001; p$horizon <- 1
  result <- .catch_compositions(.forecast_pairs(p, actual), actual, 1:3)
  expect_equal(result$pairs$obs, c(0, 1/3, 2/3))
  expect_equal(sum(result$pairs$pred), 1)
  expect_equal(result$mean_age$obs, 8/3)
  expect_equal(result$mean_age$pred, 9/4)
  expect_true(all(c("composition", "proportion", "mean_age") %in% result$scores$metric))
  actual$obs[actual$year == 2002 & actual$age == 1] <- NA
  expect_equal(nrow(.catch_compositions(.forecast_pairs(p, actual), actual, 1:3)$pairs), 0L)
})

test_that("overall MASE weights ages, surveys and horizons separately", {
  d <- data.frame(survey = c("A", "A", "B", "A", "B"), age = c(1, 2, 1, 1, 1),
    fold = 2000, horizon = c(1, 1, 1, 2, 2), obs = 1,
    pred = exp(-c(1, 3, 10, 4, 2)), baseline = exp(-1))
  d$year <- d$fold + d$horizon
  actual <- d; actual$is_proj <- FALSE
  x <- .index_score_tables(d, actual)$scores
  overall <- subset(x, level == "overall" & is.na(horizon))
  expect_equal(overall$MASE, mean(c(mean(c(2, 10)), mean(c(4, 2)))))
  expect_equal(overall$RMSE, sqrt(mean(c(1, 3, 10, 4, 2)^2)))
})

test_that("delta log scores and the zero Brier score use their own validation cases", {
  actual <- data.frame(year = rep(2000:2001, each = 6), age = rep(1:6, 2),
    survey = rep(c(rep("S", 5), "T"), 2), obs = c(rep(1, 6), 6, 0, NA, 0, 4, 4), is_proj = FALSE,
    family = c(rep("delta_lognormal", 5), "lognormal"))
  p <- subset(actual, year == 2001)
  p$is_proj <- TRUE
  p$fold <- 2000
  p$horizon <- 1
  p$pred <- c(2, 2, 2, 0, 2, 4)
  p$log_pred_positive <- log(rep(3, 6))
  p$p_zero <- c(.25, .4, .1, 1, .7, NA)
  p$is_structural <- p$age == 4
  matched <- .forecast_pairs(p, actual)
  scores <- .index_score_tables(matched, actual)$scores
  overall <- subset(scores, level == "overall" & is.na(horizon))
  expect_equal(overall$n, 3L)
  expect_equal(overall$RMSE, sqrt(mean(log(c(6/3, 4/3, 1))^2)))
  expect_equal(overall$MASE, mean(c(mean(log(c(6/3, 4/3)) / log(c(6, 4))), 0)))
  expect_equal(compute_hindcast_rmse(matched), overall$RMSE)
  expect_equal(compute_hindcast_rmse(.bind_rows(list(actual, p))), overall$RMSE)
  expect_equal(.index_score_tables(matched, actual)$aggregate$survey, "T")

  occurrence <- .occurrence_forecasts(matched)
  expect_equal(sort(occurrence$age), c(1, 2, 5))
  brier <- subset(.occurrence_score_table(occurrence, "survey"), level == "overall" & is.na(horizon))
  expect_equal(brier$n, 3L)
  expect_equal(brier$n_zero, 1L)
  expect_equal(brier$Brier, mean(c(.25^2, (1 - .4)^2, .7^2)))
})
