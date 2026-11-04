source("analysis/comp_assessments/R/pool_source_states.R")
testthat::test_that("source pooling sums counts and fleet F before weighting mortality", {
  n <- expand.grid(year = 2000L, age = 1:2, sex = c("female", "male"), stringsAsFactors = FALSE)
  n$type <- "population"; n$measure <- "numbers_at_age"
  n$value <- c(10, 20, 30, 40); n$unit <- "thousand fish"
  n$fleet <- NA_character_; n$region <- "coast"; n$season <- "annual"
  n$se <- n$lwr <- n$upr <- 1; n$notes <- "Source"
  f <- n; f$type <- "mortality"; f$measure <- "fishing_mortality_at_age"
  f$value <- c(0.1, 0.2, 0.3, 0.4); f$unit <- "per year"; f$fleet <- "fleet1"
  f2 <- f; f2$fleet <- "fleet2"; f2$value <- f$value * 2
  m <- f; m$measure <- "natural_mortality_at_age"; m$fleet <- NA_character_; m$value <- 0.2
  x <- rbind(n, f, f2, m); before <- x
  pooled <- pool_source_states(x, plus_age = 1L)
  testthat::expect_identical(x, before)
  testthat::expect_equal(pooled$value[pooled$measure == "numbers_at_age"], 100)
  testthat::expect_equal(pooled$value[pooled$measure == "fishing_mortality_at_age"], 0.9)
  testthat::expect_equal(pooled$value[pooled$measure == "natural_mortality_at_age"], 0.2)
  testthat::expect_true(all(is.na(pooled$se)))
  testthat::expect_true(all(is.na(pooled$lwr)))
  testthat::expect_true(all(pooled$age == 1))
  mixed <- rbind(x, n[1, ]); mixed$sex[nrow(mixed)] <- "combined"
  testthat::expect_error(pool_source_states(mixed), "Mixed combined")
  testthat::expect_error(pool_source_states(rbind(x, n[1, ])), "unique year-age-sex")
  conflict <- rbind(x, f[1, ]); conflict$fleet[nrow(conflict)] <- NA_character_
  testthat::expect_error(pool_source_states(conflict), "total F or component")
  missing_fleet <- pool_source_states(x[-which(x$fleet %in% "fleet2")[1], ])
  testthat::expect_true(any(is.na(missing_fleet$value[missing_fleet$measure == "fishing_mortality_at_age"])))
  testthat::expect_error(pool_source_states(x[-1, ]), "Both source sexes")
})
