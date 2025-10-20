
fit <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14, silent = TRUE)
vals <- as.list(fit$sdrep, "Estimate", report = TRUE)
sds <- as.list(fit$sdrep, "Std. Error", report = TRUE)


## tidy_array ----

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


## tidy_mat ----

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


## tidy_est ----

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


## tidy_obs_pred ----

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

test_that("tidy_obs_pred adds osa residuals", {
  n_observed <- length(fit$dat$observed)
  out <- tidy_obs_pred(fit, add_osa_res = TRUE, trace = FALSE, parallel = TRUE, reverse = FALSE)
  expect_true("osa_res" %in% names(out$catch))
  expect_true("osa_res" %in% names(out$index))
})


## tidy_rep ----

test_that("tidy_rep tidies matrices and vectors, and adds is_proj", {
  tr <- tidy_rep(fit)

  expect_type(tr, "list")
  expect_true(length(tr) > 0)

  # Should include known report matrices and at least one vector (e.g., ssb)
  expect_true(all(c("N","M","mu_M","F","mu_F","Z","ssb_mat") %in% names(tr)))
  expect_true("ssb" %in% names(tr))

  # Matrices → data frames with (year, age, est, is_proj)
  expect_s3_class(tr$N, "data.frame")
  expect_true(all(c("year","age","est","is_proj") %in% names(tr$N)))
  expect_equal(nrow(tr$N), length(fit$rep$N))  # row count equals matrix length

  # Vectors (length = #years) → data frames with (year, est, is_proj)
  expect_s3_class(tr$ssb, "data.frame")
  expect_true(all(c("year","est","is_proj") %in% names(tr$ssb)))
  expect_false("age" %in% names(tr$ssb))
  expect_equal(nrow(tr$ssb), length(fit$dat$years))

  # is_proj should match the year-level is_proj from the fit
  expect_identical(tr$ssb$is_proj, fit$dat$is_proj)
  expect_true(all(tr$N$is_proj == (tr$N$year %in% fit$dat$years[fit$dat$is_proj])))
})



## tidy_sdrep ----

test_that("tidy_sdrep extracts, transforms, and renames series", {
  trends <- tidy_sdrep(fit, interval = 0.95)
  expect_true(all(c("ssb","recruitment","abundance") %in% names(trends)))
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


## tidy_pop ----

test_that("tidy_pop concatenates tidy_sdrep and tidy_rep without duplicates", {
  out <- tidy_pop(fit)

  # Should include key names from rep and sdrep outputs
  expect_true(all(c("ssb", "N", "M", "mu_M", "F", "mu_F", "Z", "ssb_mat") %in% names(out)))
  expect_s3_class(out$N, "data.frame")
  expect_true(all(c("year","age","est") %in% names(out$N)))

  # No duplicate names
  expect_equal(length(names(out)), length(unique(names(out))))
})


## tidy_par ----

test_that("tidy_par returns nested fixed/random with expected columns", {

  tp <- tidy_par(fit)

  # Top-level structure
  expect_type(tp, "list")
  expect_true(all(c("fixed","random") %in% names(tp)))

  # Fixed table has core columns
  expect_s3_class(tp$fixed, "data.frame")
  expect_true(all(c("par","est","se","lwr","upr") %in% names(tp$fixed)))

  # Random names match RTMB random blocks
  expect_setequal(names(tp$random), fit$obj$env$.random)
  # Each random entry is a data.frame with core columns
  for (nm in names(tp$random)) {
    expect_s3_class(tp$random[[nm]], "data.frame")
    expect_true(all(c("par","est","se","lwr","upr") %in% names(tp$random[[nm]])))
  }
})

test_that("tidy_par applies back-transforms for log_ and logit_ prefixes", {

  tp <- tidy_par(fit)

  # sd_* rows should be positive (exp back-transform from log_sd_*)
  sd_rows <- subset(tp$fixed, grepl("^sd_", par))
  if (nrow(sd_rows)) {
    expect_true(all(is.finite(sd_rows$est)))
    expect_true(all(sd_rows$est > 0))
  }

  # phi_f (from logit_phi_f) should be within (0,1)
  phi_rows <- subset(tp$fixed, par == "phi_f")
  if (nrow(phi_rows)) {
    expect_true(all(phi_rows$est > 0 & phi_rows$est < 1))
    expect_true(all(phi_rows$lwr >= 0 & phi_rows$upr <= 1))
  }
})

test_that("tidy_par preserves coefficient names for named vectors and dims for matrices", {
  tp <- tidy_par(fit)

  # Named vectors like q → expect coef column present and not all NA
  q_rows <- subset(tp$fixed, par == "q")
  if (nrow(q_rows)) {
    expect_true("coef" %in% names(q_rows))
    expect_true(any(!is.na(q_rows$coef)))
  }

  # Matrix random effect (e.g., log_f) → expect year & age columns present
  if ("log_f" %in% names(tp$random)) {
    df <- tp$random$log_f
    expect_true(all(c("year","age") %in% names(df)))
    # row count equals number of cells in the reported F matrix
    expect_equal(nrow(df), length(fit$rep$F))
  }
})


## stack_nested ----

test_that("stack_nested adds id column first and type.converts labels (label_type = 'auto')", {
  x <- list(
    `2023` = list(ssb = data.frame(year = 1:2, est = 1:2),
                  N   = data.frame(year = 1:2, age = 2:3, est = 5:6)),
    `2024` = list(ssb = data.frame(year = 1:2, est = 11:12),
                  N   = data.frame(year = 1:2, age = 2:3, est = 15:16))
  )
  out <- stack_nested(x, label = "retro_year", label_type = "auto")

  # id column exists, is first, and numeric (coerced)
  expect_true("retro_year" %in% names(out$ssb))
  expect_identical(names(out$ssb)[1], "retro_year")
  expect_type(out$ssb$retro_year, "integer")  # type.convert made it numeric
  expect_equal(sort(unique(out$ssb$retro_year)), c(2023, 2024))
})

test_that("stack_nested label_type = 'character' keeps labels as character", {
  x <- list(`A` = list(ssb = data.frame(est = 1)), `B` = list(ssb = data.frame(est = 2)))
  out <- stack_nested(x, label = "sim", label_type = "character")

  expect_true("sim" %in% names(out$ssb))
  expect_type(out$ssb$sim, "character")
  expect_equal(sort(unique(out$ssb$sim)), c("A", "B"))
})

test_that("stack_nested label_type = 'factor' coerces to factor", {
  x <- list(`one` = list(ssb = data.frame(est = 1)), `two` = list(ssb = data.frame(est = 2)))
  out <- stack_nested(x, label = "grp", label_type = "factor")

  expect_true(is.factor(out$ssb$grp))
  expect_equal(levels(out$ssb$grp), c("one", "two"))
})

test_that("stack_nested label_type = 'numeric' coerces to numeric (with NA for non-numeric)", {
  x <- list(`2020` = list(ssb = data.frame(est = 1)), `abc` = list(ssb = data.frame(est = 2)))
  out <- stack_nested(x, label = "run", label_type = "numeric")

  expect_true(is.numeric(out$ssb$run))
  # first label converts to number, second becomes NA
  expect_true(all(is.na(out$ssb$run) | out$ssb$run == 2020))
})

test_that("stack_nested stacks all subtables across outer list", {
  x <- list(
    A = list(ssb = data.frame(year = 1:2, est = 1:2),
             N   = data.frame(year = 1:2, age = 2:3, est = 5:6)),
    B = list(ssb = data.frame(year = 1:2, est = 11:12)) # 'N' missing here
  )
  out <- stack_nested(x, label = "model")
  expect_true("ssb" %in% names(out))
  expect_true("N" %in% names(out))  # retained for that model even if missing in others
})


## tidy_tam ----

# Convenience variants for tests
fit_2024 <- fit
fit_2023 <- update(fit, years = min(fit$dat$years):(max(fit$dat$years) - 1))
fit_2022 <- update(fit, years = min(fit$dat$years):(max(fit$dat$years) - 2))

test_that("tidy_tam: single model via ... adds no label column", {
  out  <- tidy_tam(fit)  # default label = 'model'
  expect_false("model" %in% names(out$obs_pred$catch))
  expect_false("model" %in% names(out$obs_pred$index))
  expect_false("model" %in% names(out$pop$N))
  expect_false("model" %in% names(out$pop$ssb))
})

test_that("tidy_tam: multiple models via ... use object names as labels (and label is first column)", {
  out <- tidy_tam(fit_2023, fit_2024)  # label = 'model'
  # obs_pred
  expect_true("model" %in% names(out$obs_pred$catch))
  expect_identical(names(out$obs_pred$catch)[1], "model")
  expect_equal(sort(unique(out$obs_pred$catch$model)), c("fit_2023", "fit_2024"))
  # pop
  expect_true("model" %in% names(out$pop$N))
  expect_identical(names(out$pop$N)[1], "model")
  expect_equal(sort(unique(out$pop$N$model)), c("fit_2023", "fit_2024"))
})

test_that("tidy_tam: model_list must be named; names are used as labels (even length 1) and numeric-like names are coerced", {
  # Unnamed -> error
  expect_error(tidy_tam(model_list = list(fit_2024)))

  # Named (length 1) -> label present and equals the name, first column, numeric conversion
  out1 <- tidy_tam(model_list = list(`2024` = fit_2024))
  expect_true("model" %in% names(out1$pop$N))
  expect_identical(names(out1$pop$N)[1], "model")
  expect_type(out1$pop$N$model, "integer")   # type.convert coerced it
  expect_equal(unique(out1$pop$N$model), 2024)

  # Named (length >1) -> labels in order of names
  out2 <- tidy_tam(model_list = list(`2023` = fit_2023, `2024` = fit_2024))
  expect_true("model" %in% names(out2$obs_pred$catch))
  expect_equal(sort(unique(out2$obs_pred$catch$model)), c(2023, 2024))
  expect_type(out2$obs_pred$catch$model, "integer")
})

test_that("tidy_tam: custom label name is respected across obs_pred and pop", {
  out <- tidy_tam(fit_2024, fit_2023, label = "fold")
  expect_true("fold" %in% names(out$obs_pred$index))
  expect_identical(names(out$obs_pred$index)[1], "fold")
  expect_true("fold" %in% names(out$pop$ssb))
  expect_identical(names(out$pop$ssb)[1], "fold")
  expect_equal(sort(unique(out$pop$ssb$fold)), c("fit_2023", "fit_2024"))
})

test_that("tidy_tam: outputs are lists of data.frames and preserve is_proj", {
  out <- tidy_tam(fit_2024)
  # structure
  expect_true(is.list(out$obs_pred) && is.list(out$pop))
  expect_true(all(vapply(out$obs_pred, is.data.frame, logical(1))))
  expect_true(all(vapply(out$pop,      is.data.frame, logical(1))))
  # is_proj presence (comes from tidy_* internals)
  expect_true("is_proj" %in% names(out$obs_pred$catch))
  expect_true("is_proj" %in% names(out$pop$ssb))
  expect_type(out$obs_pred$catch$is_proj, "logical")
  expect_type(out$pop$ssb$is_proj, "logical")
})

test_that("tidy_tam: interval argument is forwarded to tidy_pop (lwr/upr present)", {
  out <- tidy_tam(fit_2024, interval = 0.80)
  # sdreport-derived series like 'ssb' should have lwr/upr
  expect_true(all(c("est", "lwr", "upr") %in% names(out$pop$ssb)))
})

test_that("tidy_tam: stacking uses intersection of components across models", {
  # Just a smoke check that all returned subtables are data.frames even if models differ
  out <- tidy_tam(fit_2022, fit_2023, fit_2024)
  expect_true(all(vapply(out$obs_pred, is.data.frame, logical(1))))
  expect_true(all(vapply(out$pop,      is.data.frame, logical(1))))
})

test_that("tidy_tam: includes fixed_par and random_par with correct structure and labeling", {
  out <- tidy_tam(fit_2023, fit_2024)

  # Ensure the new elements exist
  expect_true(all(c("fixed_par", "random_par") %in% names(out)))

  # fixed_par should be a data.frame
  expect_s3_class(out$fixed_par, "data.frame")
  expect_true(all(c("par", "est", "se", "lwr", "upr") %in% names(out$fixed_par)))

  # random_par should be a list of data.frames
  expect_type(out$random_par, "list")
  expect_true(all(vapply(out$random_par, is.data.frame, logical(1))))

  # Label column consistency
  expect_true("model" %in% names(out$fixed_par))
  expect_identical(names(out$fixed_par)[1], "model")

  # Each random block should include the label as well
  expect_true(all(vapply(out$random_par, function(df) "model" %in% names(df), logical(1))))

  # Model labels should match expected values
  expect_equal(sort(unique(out$fixed_par$model)), c("fit_2023", "fit_2024"))

  # Key parameter names should appear
  expect_true(any(grepl("sd_*", out$fixed_par$par)))
  expect_true(any(grepl("q", out$fixed_par$par)))

  # Each random block should contain numeric estimates
  expect_true(all(vapply(out$random_par, function(df) is.numeric(df$est), logical(1))))
})

test_that("tidy_tam: single model returns fixed_par and random_par without label column", {
  out <- tidy_tam(fit_2024)
  expect_true(is.data.frame(out$fixed_par))
  expect_false("model" %in% names(out$fixed_par))
  expect_true(is.list(out$random_par))
  expect_true(all(vapply(out$random_par, is.data.frame, logical(1))))
})

