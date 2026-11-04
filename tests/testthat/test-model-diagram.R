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
  expect_match(html, "document.currentScript.closest('.tam-model-structure')", fixed = TRUE)
  expect_false(grepl("document.querySelector", html, fixed = TRUE))
  all_refs <- htmltools::renderTags(plot_model(reference))$html
  expect_match(all_refs, "unavailable", fixed = TRUE)
  expect_false(grepl("<svg", all_refs, fixed = TRUE))
  expect_error(plot_model(default_fit, model_list = fits), "either")
})

test_that("Model structure opens first and groups its diagram with function calls", {
  template <- readLines(system.file("rmd", "vis_tam.Rmd", package = "tinyAM"))
  structure <- readLines(system.file("rmd", "model_structure.Rmd", package = "tinyAM"))
  expect_true(any(grepl('child = "model_structure.Rmd"', template, fixed = TRUE)))
  expect_true(any(grepl("plot_model(model_list = fits)", structure, fixed = TRUE)))
  expect_lt(which(grepl('child = "model_structure.Rmd"', template, fixed = TRUE)),
    which(grepl('child = "fisheries.Rmd"', template, fixed = TRUE)))
  expect_lt(which(grepl('child = "model_structure.Rmd"', template, fixed = TRUE)),
    which(grepl('child = if (is.null(background))', template, fixed = TRUE)))
  expect_true(any(grepl("Row {.tabset}", structure, fixed = TRUE)))
  expect_true(all(c("### Diagram", "### Function calls") %in% structure))
  tables <- readLines(system.file("rmd", "tables.Rmd", package = "tinyAM"))
  expect_false(any(grepl("### Function calls", tables, fixed = TRUE)))
  expect_false(any(grepl('Settings {data-navmenu="Tables"}', tables, fixed = TRUE)))
})

test_that("diagram nodes inherit typography and explain colours and arrows below", {
  html <- htmltools::renderTags(plot_model(default_fit))$html
  expect_match(html, '<section class="tam-model-node"', fixed = TRUE)
  expect_match(html, 'data-process="stochastic"', fixed = TRUE)
  expect_match(html, 'data-process="deterministic"', fixed = TRUE)
  expect_match(html, 'data-process="neutral"', fixed = TRUE)
  expect_false(grepl("font-family:Arial", html, fixed = TRUE))
  expect_false(grepl("<text ", html, fixed = TRUE))
  expect_false(grepl("Stored specification for", html, fixed = TRUE))
  expect_match(html, "Navy: annual stochastic process. Orange: no annual process deviations.", fixed = TRUE)
  expect_lt(regexpr("Initial parameters and observations are separate.", html, fixed = TRUE)[1],
    regexpr("Arrows show biological and observation relationships", html, fixed = TRUE)[1])
})

test_that("responsive connectors avoid nodes and follow selection and visibility changes", {
  node <- Sys.which("node"); skip_if(!nzchar(node))
  result <- system2(node, c(shQuote(test_path("model-structure-controls.js")),
    shQuote(system.file("js", "model-structure.js", package = "tinyAM"))), stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "Responsive model structure controls passed", fixed = TRUE)
})
