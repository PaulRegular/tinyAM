population_fixture <- function() {
  data <- cod_obs
  data$pop_covariates <- data$weight[c("year", "age")]
  data$pop_covariates$M_assumption <- .3
  data$pop_covariates$group <- factor(ifelse(data$pop_covariates$age < 5, "young", "old"))
  data$pop_covariates$temp <- ifelse(data$pop_covariates$age == 2, sin(data$pop_covariates$year), NA)
  data
}

test_that("tam_data adds only a class and validates again after edits", {
  data <- population_fixture()
  attr(data, "provenance") <- "fixture"
  x <- as_tam_data(data)
  expect_s3_class(x, "tam_data")
  expect_identical(unclass(x), unclass(data))
  expect_identical(as_tam_data(x), x)
  expect_invisible(check_tam_data(x))
  expect_identical(check_obs(x), check_tam_data(x))
  a <- make_test_dat(data = data, years = 1983:1990, ages = 2:6)
  b <- make_test_dat(data = x, years = 1983:1990, ages = 2:6)
  expect_identical(make_par(a), make_par(b))
  expect_equal(nll_fun(make_par(a), a), nll_fun(make_par(b), b))
  x$pop_covariates$year[1] <- NA
  expect_error(make_test_dat(data = x), "pop_covariates")
})

test_that("population keys are complete and cannot extend the model range", {
  x <- population_fixture()
  x$pop_covariates <- x$pop_covariates[-1, ]
  expect_error(check_tam_data(x), "coverage is incomplete")
  x <- population_fixture()
  x$pop_covariates <- rbind(x$pop_covariates, x$pop_covariates[1, ])
  expect_error(check_tam_data(x), "duplicate")
  x <- population_fixture()
  x$pop_covariates$year[1] <- max(x$catch$year) + 1
  expect_error(check_tam_data(x), "year-age range")
  x <- population_fixture()
  x$pop_covariates$age <- as.character(x$pop_covariates$age)
  expect_error(check_tam_data(x), "finite numeric integer")
})

test_that("all population formulas use keyed covariates and active ages", {
  x <- population_fixture()
  args <- list(data = x, years = 1983:1990, ages = 2:6,
    N_settings = list(sd_form = ~ group, rec_form = ~ temp + iid(year)),
    F_settings = list(process = "iid", mu_form = ~ group, sd_form = ~ group),
    M_settings = list(mu_supplied = ~ M_assumption))
  a <- do.call(prepare_tam, args)
  args$data$pop_covariates <- x$pop_covariates[rev(seq_len(nrow(x$pop_covariates))), ]
  b <- do.call(prepare_tam, args)
  expect_equal(a$F_modmat, b$F_modmat)
  expect_equal(a$rec$matrix, b$rec$matrix)
  expect_identical(a$process_sd, b$process_sd)
  expect_equal(nll_fun(make_par(a), a), nll_fun(make_par(b), b))
  x$weight$missing_variable <- 1
  expect_error(make_test_dat(data = x, M_settings = list(mu_supplied = ~ missing_variable)),
    "M_settings\\$mu_supplied.*data\\$pop_covariates")
  expect_error(make_test_dat(data = x, F_settings = list(process = "iid", mu_form = ~ missing_variable)), "F_settings")
  expect_error(make_test_dat(data = x, N_settings = list(rec_form = ~ missing_variable + iid(year))), "N_settings")
})

test_that("plus ages use explicit covariates and projections hold terminal rows", {
  x <- population_fixture()
  x$pop_covariates$M_assumption <- .1 + x$pop_covariates$age / 100
  before <- x
  dat <- make_test_dat(data = x, years = 1983:1990, ages = 2:6,
    M_settings = list(mu_supplied = ~ M_assumption), proj_settings = list(n_proj = 2))
  expect_equal(exp(matrix(dat$log_mu_supplied_m, 10, 5)[, 5]), rep(.16, 10))
  expect_false("obs" %in% names(dat$obs$pop_covariates))
  terminal <- subset(dat$pop_data, year == 1990)
  future <- subset(dat$pop_data, year == 1991)
  expect_identical(terminal$M_assumption, future$M_assumption)
  expect_identical(levels(terminal$group), levels(future$group))
  expect_identical(x, before)
  expect_identical(tinyAM:::.plus_fun(x, 6)$pop_covariates, x$pop_covariates)
  x$pop_covariates <- NULL
  dat <- make_test_dat(data = x, years = 1983:1990, ages = 2:6)
  expect_equal(names(dat$pop_data), c("year", "age"))
  expect_true(is.finite(nll_fun(make_par(dat), dat)))
})

test_that("printing shows compact table previews", {
  x <- as_tam_data(population_fixture())
  messages <- capture.output(text <- capture.output(print(x, n = 2)), type = "message")
  expect_true(any(grepl("pop_covariates", messages)))
  expect_lt(length(text), 65)
  expect_invisible(print(x, n = 1))
  expect_error(print(x, n = 0), "positive integer")
})
