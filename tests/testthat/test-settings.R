test_that("partial preparation uses complete defaults without changing mathematics", {
  partial <- make_test_dat(N_settings = list(init = "free"), F_settings = list(mean_ages = 3:8))
  explicit <- make_test_dat(N_settings = list(process = "iid", init = "free", sd_form = ~ 1), rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "rw", mu_form = NULL, sd_form = ~ 1, mean_ages = 3:8))
  expect_identical(partial$N_settings$process, explicit$N_settings$process)
  expect_identical(partial$N_settings$init, explicit$N_settings$init)
  expect_equal(partial$rec$matrix, explicit$rec$matrix)
  expect_identical(make_par(partial), make_par(explicit))
  expect_identical(nll_fun(make_par(partial), partial), nll_fun(make_par(explicit), explicit))
  expect_identical(make_par(make_test_dat(N_settings = list())), default_par)
})

test_that("settings merge preserves NULL and replaces whole formulas and vectors", {
  base <- tinyAM:::.settings_defaults("M_settings")
  base$mu_form <- ~ factor(age)
  base$mean_ages <- 3:8
  merged <- tinyAM:::.merge_settings(base, list(mu_form = NULL, mean_ages = 4:6), "M_settings")
  expect_null(merged$mu_form)
  expect_identical(merged$mean_ages, 4:6)
  expect_identical(merged$mu_supplied, base$mu_supplied)
  expect_identical(tinyAM:::.merge_settings(base, list(), "M_settings"), base)
  expect_null(tinyAM:::.merge_settings(base, list(mu_supplied = NULL), "M_settings")$mu_supplied)
})

test_that("all settings lists reject malformed entries and missing required values", {
  for (name in tinyAM:::.settings_names) {
    base <- if (name == "proj_settings") NULL else tinyAM:::.settings_defaults(name)
    expect_error(tinyAM:::.merge_settings(base, list(1), name), "unique, nonempty")
    expect_error(tinyAM:::.merge_settings(base, list(typo = 1), name), "Unknown")
    expect_error(tinyAM:::.merge_settings(base, setNames(list(1, 2), c("x", "x")), name), "unique, nonempty")
    expect_error(tinyAM:::.merge_settings(base, setNames(list(1), NA_character_), name), "unique, nonempty")
  }
  expect_error(make_test_dat(N_settings = list(process = NULL)), "non-NULL")
  expect_error(make_test_dat(index_settings = list(q_form = NULL)), "non-NULL")
  expect_error(make_test_dat(F_settings = list(process = "iid")), "mu_form is required")
})

test_that("updates merge stored inputs and return complete effective calls", {
  before <- default_fit$refit_args
  cl <- update(default_fit, N_settings = list(init = "free"), evaluate = FALSE)
  expect_identical(cl$N_settings$process, "iid")
  expect_identical(cl$N_settings$init, "free")
  expect_identical(default_fit$refit_args$rec_settings, before$rec_settings)
  expect_identical(default_fit$refit_args, before)
  local_mocked_bindings(fit_tam = function(...) list(refit_args = list(...), call = NULL))
  result <- update(default_fit, N_settings = list(init = "free"), index_settings = list(q_link = "logit"))
  expect_identical(result$refit_args$index_settings$q_form, before$index_settings$q_form)
  expect_identical(result$refit_args$index_settings$q_link, "logit")
  expect_identical(result$call$N_settings, result$refit_args$N_settings)
  class(result) <- c("tam_fit", "list")
  next_fit <- update(result, rec_settings = list(process = "iid", mu_form = ~ 1))
  expect_identical(next_fit$refit_args$N_settings$init, "free")
  expect_identical(next_fit$refit_args$rec_settings$mu_form, ~ 1)
  expect_identical(default_fit$refit_args, before)
})

test_that("updates retain requested formulas rather than preparation-adjusted formulas", {
  object <- default_fit
  object$refit_args$index_settings$sd_form <- ~ 1
  object$refit_args$index_settings$sd_supplied <- ~ I(.2)
  object$dat$index_settings$sd_form <- ~ 0
  cl <- update(object, index_settings = list(fill_missing = FALSE), evaluate = FALSE)
  expect_identical(cl$index_settings$sd_form, ~ 1)
  expect_identical(cl$index_settings$sd_supplied, ~ I(.2))
})

test_that("projection settings can be enabled, partially changed and disabled", {
  cl <- update(default_fit, proj_settings = list(n_proj = 2), evaluate = FALSE)
  expect_identical(cl$proj_settings, list(n_proj = 2, n_mean = 3, F_mult = 1))
  object <- default_fit
  object$refit_args$proj_settings <- list(n_proj = 4, n_mean = 2, F_mult = .8)
  cl <- update(object, proj_settings = list(n_proj = 1), evaluate = FALSE)
  expect_identical(cl$proj_settings, list(n_proj = 1, n_mean = 2, F_mult = .8))
  expect_null(update(object, proj_settings = NULL, evaluate = FALSE)$proj_settings)
  expect_equal(make_test_dat(proj_settings = list())$proj_settings$n_proj, 0)
  expect_equal(sum(make_test_dat(proj_settings = list(n_proj = 2))$is_proj), 2)
  expect_error(make_test_dat(proj_settings = list(n_mean = 2)), "n_proj")
  expect_error(make_test_dat(proj_settings = list(n_proj = NA)), "integer")
})

test_that("Usage shows every supported settings entry and optional projections", {
  expected <- list(
    N_settings = c("process", "init", "sd_form"),
    rec_settings = c("process", "mu_form", "sd", "phi"),
    F_settings = c("process", "mu_form", "sd_form", "mean_ages"),
    M_settings = c("process", "mu_form", "mu_supplied", "sd_form", "age_breaks", "first_dev_year", "mean_ages"),
    catch_settings = c("sd_form", "sd_supplied", "fill_missing"),
    index_settings = c("sd_form", "sd_supplied", "q_form", "q_link", "fill_missing"))
  for (fun in c("fit_tam", "prepare_tam")) {
    usage <- tinyAM:::.settings_usage(fun)
    call <- parse(text = sub("^@usage\n", "", usage))[[1L]]
    for (name in names(expected)) {
      entries <- eval(call[[name]], envir = asNamespace("tinyAM"))
      expect_identical(names(entries), expected[[name]])
      expect_equal(entries, tinyAM:::.settings_defaults(name))
    }
    expect_null(call$proj_settings)
    expect_match(usage, "# proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)", fixed = TRUE)
  }
})

test_that("recruitment settings validate and merge independently of survival", {
  expect_error(make_test_dat(N_settings = list(rec_form = ~ rw(year))), "has moved")
  expect_error(make_test_dat(rec_settings = list(process = "iid")), "mu_form is required")
  expect_error(make_test_dat(rec_settings = list(process = "ar1")), "mu_form is required")
  expect_error(make_test_dat(rec_settings = list(process = "off")), "must be")
  expect_error(make_test_dat(rec_settings = list(sd = 0)), "positive")
  expect_error(make_test_dat(rec_settings = list(sd = NA_real_)), "finite")
  expect_error(make_test_dat(rec_settings = list(phi = .5)), "only used")
  expect_error(make_test_dat(rec_settings = list(process = "ar1", mu_form = ~ 1, phi = 1)), "below")
  expect_error(make_test_dat(rec_settings = list(mu_form = ~ iid(year))), "residual process terms")
  before <- default_fit$refit_args
  call <- update(default_fit, rec_settings = list(process = "iid", mu_form = ~ 1), evaluate = FALSE)
  expect_identical(call$rec_settings, list(process = "iid", mu_form = ~ 1, sd = NULL, phi = NULL))
  expect_identical(default_fit$refit_args, before)
  expect_identical(make_test_dat(rec_settings = list())$rec_settings,
                   tinyAM:::.settings_defaults("rec_settings"))
  base <- list(process = "ar1", mu_form = ~ bh(ssb), sd = .2, phi = .5)
  changed <- tinyAM:::.merge_settings(base, list(mu_form = ~ 1, sd = NULL), "rec_settings")
  expect_identical(changed, list(process = "ar1", mu_form = ~ 1, sd = NULL, phi = .5))
})
