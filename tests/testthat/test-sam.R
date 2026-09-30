sam_fixture <- function() read_sam_files(test_path("fixtures", "sam"))

test_that("ICES matrix formats and configuration parse without SAM", {
  x <- sam_fixture()
  expect_equal(dim(x$data$catch[[1]]), c(2, 3))
  expect_equal(unname(x$data$catch[[1]][1, ]), c(10, 20, 30))
  expect_equal(unname(x$data$sw[3, ]), c(1, 2, 4))
  expect_true(all(x$data$cw == 2))
  expect_equal(unname(x$data$nm[, 2]), c(0.2, 0.3, 0.4))
  expect_equal(x$conf$keyLogFsta[1, ], 0:2)
  expect_identical(as.character(x$conf$obsCorStruct), rep("ID", 3))
  expect_identical(as.character(x$conf$obsLikelihoodFlag), rep("LN", 3))
  expect_equal(dim(x$conf$keyParScaledYA), c(0, 0))
  expect_length(x$conf$keyScaledYears, 0)
  expect_true(all(nchar(x$files$md5) == 32))
  f <- tempfile()
  writeLines(c("test", "test", "2000 2001", "1 2", "1", "1 2", "3"), f)
  expect_error(tinyAM:::.read_sam_ices(f), "dimensions|widths")
  writeLines(c("$keyVarObs", "0 1", "2"), f)
  expect_error(tinyAM:::.read_sam_conf(f), "row widths")
  writeLines(c("$minAge", "1", "$minAge", "2"), f)
  expect_error(tinyAM:::.read_sam_conf(f), "Duplicate")
  writeLines(c("$obsCorStruct", "ARBITRARY"), f)
  expect_error(tinyAM:::.read_sam_conf(f), "Unknown")
  writeLines(c("$keyScaledYears", "2000 2001", "2002"), f)
  expect_identical(tinyAM:::.read_sam_conf(f)$keyScaledYears, as.numeric(2000:2002))
  unlink(f)
})

test_that("conversion preserves observation values, effort, timing and biology", {
  x <- sam_fixture()
  obs <- sam_to_tam_obs(x)
  expect_true(check_obs(obs))
  expect_equal(nrow(obs$catch), 9)
  expect_true(all(is.na(obs$catch$obs[obs$catch$year == 2002])))
  expect_equal(obs$catch$obs[obs$catch$year == 2001 & obs$catch$age == 1], 0)
  a <- obs$index[obs$index$survey == "Survey A", ]
  expect_equal(a$obs[a$year == 2000], c(5, 0))
  expect_true(is.na(a$obs[a$year == 2001 & a$age == 1]))
  expect_true(all(a$samp_time == 0.125))
  expect_equal(a$effort[a$year == 2000], c(2, 2))
  b <- obs$index[obs$index$survey == "Survey B", ]
  expect_equal(b$obs, c(14, 15))
  expect_true(all(b$samp_time == 0.5))
  expect_equal(obs$weight$obs[obs$weight$year == 2000], c(1, 2, 4))
  expect_true(all(obs$weight$catch_weight == 2))
  expect_equal(obs$weight$M_assumption[obs$weight$age == 1], c(0.2, 0.3, 0.4))
  expect_equal(obs$maturity$obs[obs$maturity$year == 2000], c(0, 0.5, 1))
  expect_true(all(obs$weight$propF == 0 & obs$weight$propM == 0))
  expect_true(all(obs$catch$is_plus_group == (obs$catch$age == 3)))
  expect_false(any(obs$index$is_plus_group))
  expect_equal(as.character(obs$index$q_block[obs$index$age == 1]), rep("0", 5))
  expect_equal(as.character(obs$catch$sd_block[obs$catch$age > 1]), rep("1", 6))
})

test_that("formula block designs exactly preserve q, SD and M surfaces", {
  obs <- sam_to_tam_obs(sam_fixture())
  dat <- make_dat(obs, ages = 1:3, years = 2000:2002,
    F_settings = list(process = "rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~ 0 + q_block, sd_form = ~ 0 + sd_block, fill_missing = FALSE))
  q <- c(0.2, 0.6)
  expect_equal(unname(exp(drop(dat$q_modmat %*% log(q)))), q[dat$obs$index$q_key + 1L])
  c_sd <- c(0.1, 0.3)
  i_sd <- c(0.2, 0.4, 0.5)
  expect_equal(unname(exp(drop(dat$sd_catch_modmat %*% log(c_sd)))), c_sd[dat$obs$catch$sd_key + 1L])
  expect_equal(unname(exp(drop(dat$sd_index_modmat %*% log(i_sd)))), i_sd[match(dat$obs$index$sd_key, 2:4)])
  expect_equal(exp(dat$log_mu_supplied_m), dat$obs$weight$M_assumption)
  expect_true(any(dat$is_missing)) # explicit zeros reach tinyAM's existing handling
})

test_that("assumption audit separates formula sharing from state and covariance", {
  x <- sam_fixture()
  status <- function(x, setting) {
    tab <- sam_tam_assumptions(x)
    tab$tam_status[match(setting, tab$sam_setting)]
  }
  for (nm in c("keyLogFsta", "keyLogFpar", "keyVarObs", "keyVarF", "keyVarLogN", "corFlag", "obsCorStruct", "stockRecruitmentModelCode", "propF/propM")) {
    expect_identical(status(x, nm), "supported", info = nm)
  }
  expect_identical(status(x, "fbarRange"), "partially_supported")
  expect_identical(status(x, "unknownOption"), "not_checked")
  x$data$nm[1, 1] <- 0
  expect_identical(status(x, "nm.dat"), "unsupported")
  x$conf$keyLogFsta[1, 3] <- 1
  expect_identical(status(x, "keyLogFsta"), "unsupported")
  x$conf$corFlag <- 2
  expect_identical(status(x, "corFlag"), "unsupported")
  x$conf$keyVarF[1, 3] <- 1
  expect_identical(status(x, "keyVarF"), "supported") # duplicated state's later variance key is ignored by SAM
  x$conf$keyLogFsta[1, 3] <- 2
  expect_identical(status(x, "keyVarF"), "unsupported")
  x$conf$keyVarLogN[3] <- 2
  expect_identical(status(x, "keyVarLogN"), "unsupported")
  x$conf$keyVarObs[2, 1] <- 0
  expect_identical(status(x, "keyVarObs"), "partially_supported")
  x$conf$obsCorStruct[2] <- "AR"
  expect_identical(status(x, "obsCorStruct"), "unsupported")
  expect_identical(status(x, "keyCorObs"), "unsupported")
  x$conf$obsLikelihoodFlag[2] <- "ALN"
  expect_identical(status(x, "obsLikelihoodFlag"), "unsupported")
  x$conf$stockRecruitmentModelCode <- 1
  expect_identical(status(x, "stockRecruitmentModelCode"), "unsupported")
  x$data$pm[,] <- 0.5
  expect_identical(status(x, "propF/propM"), "unsupported")
  x$conf <- list()
  expect_identical(status(x, "keyLogFpar"), "not_checked")
})

test_that("fixed q cells also have an exact formula mapping", {
  x <- sam_fixture()
  x$conf$keyLogFpar[2, 2] <- -1
  obs <- sam_to_tam_obs(x)
  design <- stats::model.matrix(~ 0 + q_key_0, obs$index)
  expect_equal(unname(exp(drop(design %*% log(0.4)))), ifelse(obs$index$q_key == -1, 1, 0.4))
  audit <- sam_tam_assumptions(x)
  row <- audit[audit$sam_setting == "keyLogFpar", ]
  expect_identical(row$tam_status, "supported")
  expect_identical(row$tam_mapping, "index_settings$q_form = ~ 0 + q_key_0")
  x$conf$keyLogFpar[,] <- -1
  audit <- sam_tam_assumptions(x)
  expect_identical(audit$tam_mapping[audit$sam_setting == "keyLogFpar"], "index_settings$q_form = ~ 0")
})

test_that("unsupported fleets and executable attributes cannot be silently converted", {
  x <- sam_fixture()
  x$fleets$fleet_type[2] <- 3
  expect_error(sam_to_tam_obs(x), "Unsupported SAM fleet types")
  expect_identical(sam_tam_assumptions(x)$tam_status[4], "unsupported")
  x$fleets$fleet_type[2] <- 7
  x$conf$obsCorStruct[2] <- NA
  expect_identical(sam_tam_assumptions(x)$tam_status[4], "unsupported")
  x <- sam_fixture()
  x$fleets$fleet_type[2] <- 0
  expect_error(sam_to_tam_obs(x), "multiple fleets")
  x <- sam_fixture()
  x$data$sw <- x$data$sw[-3, , drop = FALSE]
  expect_error(sam_to_tam_obs(x), "complete modeled biological grid")
  x <- sam_fixture()
  x$conf$keyLogFpar[2, 1] <- 0.5
  expect_error(sam_to_tam_obs(x), "keys must be integers")
  f <- tempfile()
  writeLines(c("test", "test", "2000 2001", "1 2", "3", "1", "#' @weight", "#' stop('execute')"), f)
  expect_error(tinyAM:::.read_sam_ices(f), "custom R-expression")
  unlink(f)
})

test_that("saved reference extraction maps fitted states and stored predictions", {
  x <- sam_fixture()
  fit <- list(conf = x$conf, opt = list(convergence = 0),
    data = list(years = 2000:2002, fleetTypes = c(0, 2, 2),
      aux = matrix(c(2000, 1, 1, 2000, 2, 1), 2, 3, byrow = TRUE),
      logobs = log(c(10, 5))),
    pl = list(logN = log(matrix(1:9, 3)), logF = log(matrix(rep(c(0.1, 0.2, 0.3), 3), 3)), logFpar = log(c(0.2, 0.6))),
    rep = list(predObs = log(c(11, 6))),
    sdrep = list(value = stats::setNames(log(1:9), rep(c("logssb", "logR", "logfbar"), each = 3))))
  ref <- sam_reference(fit, "synthetic fitted fixture")
  expect_equal(ref$tables$N$obs, as.vector(exp(t(fit$pl$logN))))
  expect_equal(ref$tables$F$obs, rep(c(0.1, 0.2, 0.3), each = 3))
  expect_equal(ref$tables$q$est, c(0.2, 0.2, 0.6))
  expect_equal(ref$tables$SSB$est, 1:3)
  expect_equal(ref$tables$recruitment$est, 4:6)
  expect_equal(ref$tables$Fbar$est, 7:9)
  expect_equal(ref$tables$catch$pred, 11)
  expect_equal(ref$tables$index$pred, 6)
  expect_true(all(ref$availability$available))
  fit$rep <- NULL
  fit$sdrep <- NULL
  ref <- sam_reference(fit)
  expect_false(ref$availability$available[ref$availability$quantity == "SSB"])
  expect_true(all(is.na(ref$tables$catch$pred)))
  expect_false(ref$availability$available[ref$availability$quantity == "observation_predictions"])
  expect_error(sam_reference(list(logN = 1)), "initial parameters")
})
