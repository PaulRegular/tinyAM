f_curve_dat <- function(form = ~ logistic(age), process = "iid", ...) {
  prepare_tam(cod_obs,
    years = 1983:1998, ages = 2:8,
    N_settings = list(process = "off"),
    F_settings = list(process = process, mu_form = form), ...
  )
}

test_that("F curves define the log mean without changing residual densities", {
  for (process in c("iid", "ar1")) {
    for (curve in c("mono", "logistic")) {
      d <- f_curve_dat(if (curve == "mono") ~ mono(age) else ~ logistic(age), process)
      p <- make_par(d)
      p$log_mu_f[] <- log(.4)
      if (curve == "mono") {
        p$dmu_f[] <- c(.3, .2, 0, 0, .1, 0)
        mean <- .4 * exp(c(0, cumsum(p$dmu_f)))
      } else {
        term <- d$F_terms[[1L]]
        expect_null(term$parameter)
        p[[term$a50_parameter]][] <- 4
        p[[term$slope_parameter]][] <- log(1.2)
        mean <- .4 * plogis(1.2 * (d$ages - 4))
      }
      effect <- matrix(tinyAM:::.mean_effects(p, d, "F")$contribution, length(d$years))
      expect_equal(as.numeric(exp(p$log_mu_f + effect[1, ])), as.numeric(mean))
      expect_true(all(diff(mean) >= 0))
      p$log_f[] <- rep(log(mean), each = nrow(p$log_f))
      obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
      expect_true(all(is.finite(obj$gr(obj$par))))
      expect_equal(as.numeric(obj$report()$mu_F[1, ]), as.numeric(mean))
      ordinary <- f_curve_dat(~ factor(age), process)
      equivalent <- make_par(ordinary)
      shared <- intersect(names(equivalent), names(p))
      equivalent[shared] <- p[shared]
      equivalent$log_mu_f <- setNames(
        c(log(mean[1L]), log(mean[-1L] / mean[1L])),
        colnames(ordinary$F_modmat)
      )
      expect_equal(as.numeric(nll_fun(p, d)), as.numeric(nll_fun(equivalent, ordinary)), tolerance = 1e-10)
      initial <- make_par(d)
      initial_mean <- matrix(tinyAM:::.mean_effects(initial, d, "F")$contribution, length(d$years))
      expect_equal(unname(initial$log_f), unname(initial_mean))
    }
  }
})

test_that("F curves and q curves have independent parameter blocks", {
  d <- f_curve_dat(index_settings = list(q_form = ~ logistic(age)))
  p <- make_par(d)
  expect_equal(anyDuplicated(names(p)), 0L)
  expect_true(all(c(d$F_terms[[1]]$a50_parameter, d$q_terms[[1]]$a50_parameter) %in% names(p)))
  expect_false(d$F_terms[[1]]$a50_parameter %in% tinyAM:::.formula_random_parameters(d))
  expect_true(all(is.finite(RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)$gr())))
})

test_that("F curve grouping uses population columns and retains warm-start labels", {
  data <- cod_obs
  data$pop_covariates$shape <- factor(ifelse(data$pop_covariates$year < 1990, "early", "late"))
  data$pop_covariates$size <- data$pop_covariates$age * 2
  for (curve in c("mono", "logistic")) {
    form <- if (curve == "mono") ~ shape + mono(size, by = shape) else ~ shape + logistic(size, by = shape)
    d <- prepare_tam(data, years = 1983:1998, ages = 2:8, F_settings = list(process = "iid", mu_form = form))
    p <- make_par(d)
    if (curve == "mono") expect_equal(length(p$dmu_f), 12L) else expect_length(d$F_terms[[1]]$groups, 2L)
    shortened <- prepare_tam(data, years = 1985:1996, ages = 2:8, F_settings = list(process = "iid", mu_form = form))
    aligned <- tinyAM:::.merge_start_par(make_par(shortened), p)
    expect_equal(aligned$log_f, p$log_f[as.character(1985:1996), , drop = FALSE])
  }
})

test_that("F curve safeguards do not unlock Gaussian temporal overlaps", {
  expect_error(f_curve_dat(~ mono(age), "rw"), "cancel")
  expect_error(f_curve_dat(~ logistic(age), "cor_rw"), "cancel")
  expect_error(f_curve_dat(~ mono(age) + logistic(age)), "one mono")
  expect_error(f_curve_dat(~ logistic(year)), "age-derived")
  expect_error(f_curve_dat(~ 0 + logistic(age)), "baseline")
  expect_error(f_curve_dat(~ 0 + age + logistic(age)), "baseline")
  expect_error(f_curve_dat(~ factor(age) + logistic(age)), "duplicates")
  expect_error(f_curve_dat(~ factor(age) + mono(age)), "ordinary")
  expect_error(f_curve_dat(~ mono(age, by = year)), "categorical")
  expect_error(f_curve_dat(~ logistic(age) + rw(year), "ar1"), "IID residuals")
  expect_error(f_curve_dat(~ logistic(age), M_settings = list(mu_form = ~ logistic(age))), "F only")
})

test_that("F logistic support cautions remain separate from numerical diagnostics", {
  d <- f_curve_dat()
  term <- d$F_terms[[1L]]
  fit <- list(dat = d, sdrep = list(
    pdHess = TRUE,
    par.fixed = setNames(12, term$a50_parameter), cov.fixed = matrix(25)
  ))
  advice <- tinyAM:::.formula_uncertainty_advisories(fit)
  expect_equal(advice$issue, "logistic_support")
  expect_match(advice$detail, "baseline F")
  expect_match(advice$detail, "2 to 8")
  fit$sdrep$cov.fixed[] <- .01
  expect_equal(nrow(tinyAM:::.formula_uncertainty_advisories(fit)), 0L)
  groups <- tinyAM:::.fixed_parameter_groups(data.frame(
    par = c("dmu_f", "F_a50_F_logistic_age", "F_slope_F_logistic_age")
  ))
  expect_equal(nrow(groups[["Mortality means"]]), 3L)
})

test_that("F curve simulation and projection use the same mean surface", {
  for (curve in c("mono", "logistic")) {
    d <- f_curve_dat(if (curve == "mono") ~ mono(age) else ~ logistic(age), "ar1",
      proj_settings = list(n_proj = 2, n_mean = 1, F_mult = c(1, .5))
    )
    p <- make_par(d)
    p$log_mu_f[] <- log(.3)
    p$log_sd_f <- log(.15)
    p$log_r0 <- log(1e6)
    # Negligible observation noise isolates predictions from simulated states.
    p$log_sd_catch[] <- p$log_sd_index[] <- -30
    set.seed(720)
    sim <- nll_fun(p, d, simulate = TRUE)
    expect_equal(dim(sim$log_f), c(16L, 7L))
    p$log_f <- sim$log_f
    p$log_r <- sim$log_r
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    rep <- obj$report()
    expect_equal(sim$log_obs[d$is_observed], rep$log_pred[d$is_observed], tolerance = 1e-9)
    expect_true(all(is.na(sim$log_obs[d$is_missing])))
    expect_equal(rep$F[17, ], rep$F[16, ], ignore_attr = TRUE)
    expect_equal(rep$F[18, ], rep$F[16, ] * .5, ignore_attr = TRUE)
  }
})

test_that("F curve fits retain tidy scales, uncertainty and update paths", {
  grid <- expand.grid(year = 2000:2019, age = 1:8)
  data <- list(
    catch = transform(grid, obs = 1, relative_sd = .08),
    index = transform(grid, obs = 1, survey = "survey", samp_time = .5, relative_sd = .08),
    weight = transform(grid, obs = .2 * age^1.5), maturity = transform(grid, obs = plogis(age - 4))
  )
  settings <- list(
    N_settings = list(process = "off"),
    F_settings = list(process = "iid", mu_form = ~ mono(age)),
    catch_settings = list(sd_form = ~ 0, sd_supplied = ~ relative_sd),
    index_settings = list(q_form = ~ 1, sd_form = ~ 0, sd_supplied = ~ relative_sd))
  d <- do.call(prepare_tam, c(list(data = data), settings))
  p <- make_par(d)
  p$log_r0 <- log(1e6)
  p$log_sd_r <- log(.2)
  p$log_sd_f <- log(.15)
  p$log_mu_f[] <- log(.1)
  p$dmu_f[] <- c(.6, .4, .2, 0, 0, 0, 0)
  p$log_q[] <- log(.7)
  set.seed(552)
  sim <- nll_fun(p, d, simulate = TRUE)
  values <- split(exp(sim$log_obs), d$obs_map$type)
  data$catch$obs <- values$catch
  data$index$obs <- values$index
  fit <- do.call(fit_tam, c(list(
    data = data, silent = TRUE,
    start_par = list(log_r0 = log(8e5))
  ), settings))
  expect_true(fit$is_converged)
  expect_equal(fit$bounds$lower[names(fit$opt$par) == "dmu_f"], rep(0, 7))
  steps <- subset(tidy_par(fit)$fixed, par == "dmu_f")
  expect_true(all(steps$est >= 0 & is.finite(steps$se)))
  expect_equal(steps$se, as.numeric(as.list(fit$sdrep, "Std. Error")$dmu_f))
  expect_true(all(steps$se_scale == "reported"))
  expect_true(all(is.finite(tidy_pop(fit)$mu_F$se)))
  expect_equal(nrow(tidy_pop(fit)$mu_F), 20 * 8)
  before <- fit$refit_args
  call <- update(fit, F_settings = list(process = "ar1"), evaluate = FALSE)
  expect_identical(call$F_settings$mu_form, before$F_settings$mu_form)
  expect_identical(fit$refit_args, before)
  projected <- update(fit,
    proj_settings = list(n_proj = 1), silent = TRUE,
    start_par = as.list(fit$sdrep, "Estimate")
  )
  expect_equal(projected$dat$F_mono_steps, d$F_mono_steps)
  expect_true(all(is.finite(sim_tam(projected, n = 1, par_uncertainty = "none", progress = FALSE)$F$est)))
  retro <- fit_retro(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  expect_true(all(vapply(retro$fits, inherits, logical(1), "tam_fit")))
  expect_equal(nrow(retro$fits[[1]]$rep$F), 19L)
  hindcast <- fit_hindcast(fit, folds = 1, start_from_fit = TRUE, progress = FALSE)
  expect_true(all(vapply(hindcast$fits, inherits, logical(1), "tam_fit")))
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    path <- tempfile(fileext = ".html")
    reference <- fit
    class(reference) <- c("tam_ref", "list")
    reference$pop$mu_F[, c("est", "lwr", "upr", "se")] <- NA_real_
    reference$dat$F_mono_modmat <- NULL
    reference$dat$F_terms <- list(list(type = NA_character_))
    expect_no_error(vis_tam(
      model_list = list(Reference = reference, tinyAM = fit),
      output_file = path, open_file = FALSE,
      render_args = list(quiet = TRUE)
    ))
    con <- file(path, "rb")
    html <- paste(readLines(con, warn = FALSE), collapse = "\n")
    close(con)
    expect_false(grepl("<h3>Mean F</h3>", html, fixed = TRUE))
    expect_match(html, "Annual trends", fixed = TRUE)
    expect_match(html, "Age profiles", fixed = TRUE)
    expect_match(html, "Mortality means", fixed = TRUE)
    scripts <- regmatches(html, gregexpr(
      '(?s)<script type="application/json"[^>]*>.*?</script>', html,
      perl = TRUE
    ))[[1L]]
    widgets <- lapply(scripts, function(script) {
      jsonlite::fromJSON(
        sub("</script>$", "", sub("^<script[^>]*>", "", script)),
        simplifyVector = FALSE
      )
    })
    F_plots <- Filter(function(widget) any(unlist(widget$x$layout$yaxis$title) == "F"), widgets)
    expect_length(F_plots, 2L)
    expect_setequal(vapply(
      F_plots, function(widget) unlist(widget$x$layout$xaxis$title),
      character(1)
    ), c("Year", "Age"))
    for (plot in F_plots) {
      ribbons <- vapply(plot$x$data, function(trace) {
        !is.null(trace$fill) && trace$fill %in% c("toself", "tonexty")
      }, logical(1))
      expect_true(any(ribbons))
      mean_lines <- Filter(function(trace) {
        !is.null(trace$hovertemplate) &&
          any(grepl("Process mean:", trace$hovertemplate, fixed = TRUE))
      }, plot$x$data)
      expect_true(any(vapply(mean_lines, function(trace) {
        any(is.finite(unlist(trace$y)))
      }, logical(1))))
    }
  }
})
