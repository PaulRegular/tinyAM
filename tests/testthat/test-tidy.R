
test_that("tidy_array errors on non-matrix/array", {
  expect_error(tidy_array(1:3), "`x` must be a matrix or array")
})

test_that("tidy_array requires dimnames by default", {
  m <- matrix(1:4, 2, 2)
  dimnames(m) <- NULL
  expect_error(
    tidy_array(m),
    "All dimensions must have names"
  )
})

test_that("tidy_array works with require_dimnames = FALSE and sets fallback names", {
  m <- matrix(1:4, 2, 2)
  dimnames(m) <- NULL
  out <- tidy_array(m, require_dimnames = FALSE)
  expect_s3_class(out, "data.frame")
  expect_setequal(names(out), c("Var1", "Var2", "x"))
  expect_equal(nrow(out), length(m))
})

test_that("tidy_array respects value_name and converts numeric-like dimnames", {
  m <- matrix(
    1:6, nrow = 2,
    dimnames = list(age = c("2", "3"), year = c("2001", "2002", "2003"))
  )
  out <- tidy_array(m, value_name = "val")
  expect_setequal(names(out), c("age", "year", "val"))
  expect_true(is.numeric(out$age))
  expect_true(is.numeric(out$year))
  expect_equal(nrow(out), length(m))
})

test_that("tidy_mat is an alias of tidy_array", {
  m <- matrix(
    1:6, nrow = 2,
    dimnames = list(age = c("2", "3"), year = c("2001", "2002", "2003"))
  )
  a <- tidy_array(m)
  b <- tidy_mat(m)
  expect_identical(a, b)
})

test_that("tidy_array handles 3D arrays, converting only numeric-like dims", {
  a <- array(
    1:8,
    dim = c(2, 2, 2),
    dimnames = list(age = c("2", "3"), year = c("2001", "2002"), area = c("N", "S"))
  )
  out <- tidy_array(a)
  expect_equal(nrow(out), length(a))
  expect_setequal(names(out), c("age", "year", "area", "x"))
  expect_true(is.numeric(out$age))
  expect_true(is.numeric(out$year))
  expect_false(is.numeric(out$area))
})

test_that("trans_est applies transform and scale to est/lwr/upr", {
  d <- data.frame(est = log(100), lwr = log(80), upr = log(120))
  out <- trans_est(d, transform = exp, scale = 2)
  expect_equal(out$est, 100 / 2)
  expect_equal(out$lwr,  80 / 2)
  expect_equal(out$upr, 120 / 2)
})

test_that("trans_est respects transform = NULL", {
  d <- data.frame(est = 1, lwr = 0.5, upr = 2)
  out <- trans_est(d, transform = NULL, scale = 10)
  expect_equal(out$est, 0.1)
  expect_equal(out$lwr, 0.05)
  expect_equal(out$upr, 0.2)
})

fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14, silent = TRUE)
vals <- as.list(fit$sdrep, "Estimate", report = TRUE)
sds <- as.list(fit$sdrep, "Std. Error", report = TRUE)

test_that("tidy_obs_pred builds residual diagnostics for catch and index", {

  out <- tidy_obs_pred(fit)

  # Structure
  expect_true(is.list(out))
  expect_true(all(c("catch", "index") %in% names(out)))

  # Pred lengths
  expect_equal(nrow(out$catch), nrow(fit$dat$obs$catch))
  expect_equal(nrow(out$index), nrow(fit$dat$obs$index))

  # q recovered by exp
  expect_equal(out$index$q, unname(exp(fit$rep$log_q_obs)))

  # Residuals: NA where obs == 0; otherwise (log(obs) - log(pred)) / sd
  expected_catch_res <- with(out$catch, ifelse(obs == 0, NA_real_, (log(obs) - log(pred)) / sd))
  expected_index_res <- with(out$index, ifelse(obs == 0, NA_real_, (log(obs) - log(pred)) / sd))

  expect_equal(out$catch$std_res, expected_catch_res)
  expect_equal(out$index$std_res, expected_index_res)
})


test_that("tidy_rep_mats tidies key report matrices via tidy_mat", {
  mats <- tidy_rep_mats(fit)

  expect_true(is.list(mats))
  expect_setequal(names(mats), c("N","M","mu_M","F","mu_F","Z","ssb_mat"))

  # Each element is a data frame with dims + "est"
  expect_true(all(vapply(mats, inherits, logical(1), what = "data.frame")))
  expect_true(all(vapply(mats, function(df) all(c("year","age","est") %in% names(df)), logical(1))))
  # Row counts equal length of matrix
  expect_true(all(vapply(mats, nrow, integer(1)) == length(fit$rep$N)))
})


test_that("tidy_sdrep extracts, transforms, and renames series", {
  trends <- tidy_sdrep(fit, interval = 0.95)
  expect_setequal(names(trends), c("ssb"))
  expect_true(all(grepl("^log_", names(vals)))) # tidy_sdrep assumes all ADREPORTEd values are in log space
  expect_true(all(sapply(vals, length) == length(fit$dat$years))) # tidy_sdreport assumes all ADREPORTed values have a length = n_years

  # Check one series numerically
  z <- qnorm(0.975)
  d <- data.frame(est = exp(vals$log_ssb),
                  lwr = exp(vals$log_ssb - z * sds$log_ssb),
                  upr = exp(vals$log_ssb + z * sds$log_ssb))
  expect_equal(trends$ssb$year, fit$dat$years)
  expect_equal(trends$ssb$est, d$est)
  expect_equal(trends$ssb$lwr, d$lwr)
  expect_equal(trends$ssb$upr, d$upr)
})

test_that("tidy_pop concatenates tidy_sdrep and tidy_rep_mats", {
  out <- tidy_pop(fit)
  # Should include names from tidy_sdrep (log_ prefix removed) and report matrices
  expect_true(all(c("ssb", "N", "M", "mu_M", "F", "mu_F", "Z", "ssb_mat") %in% names(out)))
  expect_s3_class(out$N, "data.frame")
})


test_that("tidy_sam supplied single model via ... adds no label column", {
  out  <- tidy_tam(fit)  # default label = 'model'
  expect_false("model" %in% names(out$obs_pred$catch))
  expect_false("model" %in% names(out$obs_pred$index))
  expect_false("model" %in% names(out$pop$N))
  expect_false("model" %in% names(out$pop$ssb))
})


fit_2024 <- fit
fit_2023 <- update(fit, years = 1983:2023)
fit_2022 <- update(fit, years = 1983:2022)

test_that("tidy_sam supplied multiple models via ... use object names as labels", {
  out <- tidy_tam(fit_2023, fit_2024)  # label = 'model'
  expect_true("model" %in% names(out$obs_pred$catch))
  expect_equal(unique(out$obs_pred$catch$model), c("fit_2023", "fit_2024"))
  expect_true("model" %in% names(out$pop$N))
  expect_equal(unique(out$pop$N$model), c("fit_2023", "fit_2024"))
})

test_that("tidy_sam supplied a model_list must be named and names are used as labels (even length 1)", {
  # Unnamed -> error
  expect_error(tidy_tam(model_list = list(fit)))

  # Named (length 1) -> label present and equals the name
  out1 <- tidy_tam(model_list = list(`2024` = fit_2024))
  expect_true("model" %in% names(out1$pop$N))
  expect_equal(unique(out1$pop$N$model), "2024")

  # Named (length >1) -> labels in order of names
  out2 <- tidy_tam(model_list = list(`2023` = fit_2023, `2024` = fit_2024))
  expect_equal(unique(out2$obs_pred$catch$model), c("2023", "2024"))
})

test_that("tidy_sam supplied custom label name is respected", {
  out <- tidy_tam(fit_2024, fit_2023, label = "retro_year")
  expect_true("retro_year" %in% names(out$obs_pred$catch))
})




