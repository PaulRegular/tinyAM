test_that("model diagrams read completed settings without calculations or mutation", {
  testthat::local_mocked_bindings(
    fit_tam = function(...) stop("Unexpected refit"),
    tidy_tam = function(...) stop("Unexpected diagnostics"),
    nll_fun = function(...) stop("Unexpected objective evaluation"))
  fit <- default_fit
  original <- fit
  set.seed(251)
  seed <- .Random.seed
  view <- plot_model(fit)
  expect_s3_class(view, "shiny.tag")
  html <- htmltools::renderTags(view)$html
  expect_match(html, "Stochastic RW")
  expect_match(html, "Deterministic: process off")
  expect_match(html, "Older-age initialization: exp")
  expect_match(html, "Recruitment boundary: free")
  expect_match(html, "F is the sum of fleet F")
  expect_match(html, "tinyAM-model.html", fixed = TRUE)
  expect_match(html, "q formula:")
  expect_false(grepl("Parent SSB:", html, fixed = TRUE))
  expect_match(html, 'data-formula="survey"', fixed = TRUE)
  expect_match(html, 'data-input-component="survey"', fixed = TRUE)
  expect_match(html, 'index: q_block', fixed = TRUE)
  expect_identical(fit, original)
  expect_identical(.Random.seed, seed)
})

test_that("model diagrams show deterministic means and stock-recruit feedback", {
  dat <- make_test_dat(years = 2000:2007, ages = 2:8,
    N_settings = list(process = "off", init = "latent"),
    rec_settings = list(process = "off", mu_form = ~ bh(ssb)),
    F_settings = list(process = "off", mu_form = ~ factor(year) + factor(age)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3)))
  fit <- structure(list(dat = dat), class = "tam_fit")
  html <- htmltools::renderTags(plot_model(fit))$html
  expect_match(html, "Mean: ~factor\\(year\\) \\+ factor\\(age\\)")
  expect_match(html, "Parent SSB: bh, lag 2", fixed = TRUE)
  expect_match(html, "Parent-SSB relationship: bh ; lag 2", fixed = TRUE)
  expect_match(html, "Older-age initialization: latent", fixed = TRUE)
  expect_match(html, "Deterministic (process off)", fixed = TRUE)
  expect_false(grepl("Stochastic process:", html, fixed = TRUE))
  expect_match(html, 'data-input-component="fishing"', fixed = TRUE)
  expect_match(html, "catch: year, age", fixed = TRUE)
  spec <- .model_structure_spec(dat)
  expect_identical(unname(spec$inputs["recruitment"]), "")
  dat$rec_settings$mu_form <- ~ factor(year)
  expect_identical(unname(.model_structure_spec(dat)$inputs["recruitment"]), "")
  dat$pop_data$temperature <- seq_len(nrow(dat$pop_data))
  dat$N_settings$mu_form <- ~ temperature
  dat$M_settings$mu_form <- ~ temperature
  spec <- .model_structure_spec(dat)
  expect_equal(unname(spec$inputs[c("population", "natural")]),
    rep("pop_covariates: temperature", 2L))
  expect_equal(unname(spec$formula_labels["natural"]), "Mean effects ~temperature")
})

test_that("model diagrams retain fleet, survey and shared block definitions", {
  data <- cod_obs
  data$catch <- rbind(transform(data$catch, fleet = "Trawl", obs = .6 * obs),
    transform(data$catch, fleet = "Gillnet", obs = .4 * obs))
  data$index <- rbind(transform(data$index, survey = "RV"),
    transform(data$index, survey = "Sentinel"))
  dat <- make_test_dat(data = data, years = 2000:2007, ages = 2:8,
    F_settings = list(process = "off", mu_form = ~ 0 +
      for_fleet("Trawl", ~ factor(year) + logistic(age)) +
      for_fleet("Gillnet", ~ factor(year) + double_logistic(age))),
    M_settings = list(process = "ar1", age_breaks = c(2, 5, 9)))
  fit <- structure(list(dat = dat), class = "tam_fit")
  html <- htmltools::renderTags(plot_model(fit))$html
  expect_match(html, "Fleets: Gillnet, Trawl", fixed = TRUE)
  expect_match(html, "Surveys: RV, Sentinel", fixed = TRUE)
  expect_match(html, "for_fleet", fixed = TRUE)
  expect_match(html, "double_logistic", fixed = TRUE)
  expect_match(html, "Age blocks: 2-4, 5-8", fixed = TRUE)
  expect_match(html, "Shared age blocks share actual M", fixed = TRUE)
})

test_that("model selectors handle references and independent escaped containers", {
  reference <- structure(list(pop = list()), class = "tam_ref")
  fits <- list(Accepted = reference, `Fit <two> \"quoted\"` = default_fit)
  html <- htmltools::renderTags(plot_model(model_list = fits))$html
  expect_match(html, '<option value="2" selected="selected">', fixed = TRUE)
  expect_match(html, 'data-model-pane="1" style="display:none"', fixed = TRUE)
  expect_match(html, "unavailable for a reporting reference", fixed = TRUE)
  expect_match(html, "Fit &lt;two&gt; &quot;quoted&quot;", fixed = TRUE)
  expect_match(html, "closest(&#39;.tam-model-structure&#39;)", fixed = TRUE)
  expect_false(grepl("document.querySelector", html, fixed = TRUE))
  all_refs <- htmltools::renderTags(plot_model(reference))$html
  expect_match(all_refs, "unavailable", fixed = TRUE)
  expect_false(grepl("<svg", all_refs, fixed = TRUE))
  expect_error(plot_model(default_fit, model_list = fits), "either")
})

test_that("the dedicated Model structure page preserves the dashboard initial page", {
  template <- readLines(system.file("rmd", "vis_tam.Rmd", package = "tinyAM"))
  structure <- readLines(system.file("rmd", "model_structure.Rmd", package = "tinyAM"))
  expect_true(any(grepl('child = "model_structure.Rmd"', template, fixed = TRUE)))
  expect_true(any(grepl("plot_model(model_list = fits)", structure, fixed = TRUE)))
  expect_gt(which(grepl('child = "model_structure.Rmd"', template, fixed = TRUE)),
    which(grepl('child = "fisheries.Rmd"', template, fixed = TRUE)))
})
