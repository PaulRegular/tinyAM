.process_view_fixture <- function() {
  d <- expand.grid(year = 2000:2002, age = 2:3, model = c("sam", "ncam"))
  d$est <- sin(seq_len(nrow(d)))
  m <- data.frame(year = rep(2000:2002, 5),
    age_block = rep(c("2-4", "5-8", "2", "3-5", "6-8"), each = 3),
    model = rep(c("sam", "ncam"), c(6, 9)), est = seq(-2, 2, length.out = 15))
  rec <- data.frame(year = rep(2000:2002, 2), model = rep(c("sam", "ncam"), each = 3),
    est = c(-.3, .2, -.1, .4, -.2, 0))
  rec$lwr <- rec$est - .1; rec$upr <- rec$est + .1
  status <- data.frame(model = c("scaa", "sam", "ncam"), status = c("deterministic", "ok", "ok"),
    detail = c("No stochastic processes.", "Calculated.", "Calculated."))
  structure(list(recruitment = list(residual = rec), deviations = list(N = d, F = d, M = m),
    process_residuals = list(recruitment = rec, N = d, F = d, M = m, status = status)), label = "model")
}

test_that("process helpers return individual widgets without changing inputs or RNG", {
  out <- .process_view_fixture(); original <- serialize(out, NULL)
  set.seed(128); seed <- .Random.seed
  expect_equal(.process_models(out), c("scaa", "sam", "ncam"))
  for (viewer in list(plot_process_deviations, plot_process_residuals)) {
    for (component in c("Recruitment", "NFM")) {
      p <- viewer(out, component = component)
      expect_s3_class(p, "plotly")
      expect_gt(length(plotly::plotly_build(p)$x$data), 0L)
    }
    expect_error(viewer(out, component = "invalid"), "arg")
  }
  expect_identical(serialize(out, NULL), original)
  expect_identical(.Random.seed, seed)
  expect_error(plot_process_deviations(list()), "tidy_tam")
})

test_that("empty recruitment is a centred Plotly annotation and NFM retains empty models", {
  out <- structure(list(recruitment = list(), deviations = list(),
    process_residuals = list(status = data.frame(status = "deterministic", detail = "No stochastic processes."))), label = NULL)
  for (viewer in list(plot_process_deviations, plot_process_residuals)) {
    p <- plotly::plotly_build(viewer(out, component = "Recruitment"))
    expect_identical(p$x$layout$annotations[[1]]$text, "Unavailable for this model component")
    expect_equal(p$x$layout$annotations[[1]]$x, .5)
    expect_equal(p$x$layout$annotations[[1]]$y, .5)
    p <- viewer(out)
    spec <- p$jsHooks$render[[1L]]$data
    expect_equal(spec$models, "Model")
    if (identical(viewer, plot_process_residuals)) {
      expect_false(any(vapply(spec$components[[1]]$cells, `[[`, logical(1), "available")))
    } else expect_true(all(is.na(unlist(spec$metrics[[1]]$cells[[1]]$linear))))
  }
})

test_that("recruitment uses supplied signed estimates with intervals only for deviations", {
  out <- .process_view_fixture()
  fitted <- plotly::plotly_build(plot_process_deviations(out, "Recruitment"))
  sampled <- plotly::plotly_build(plot_process_residuals(out, "Recruitment"))
  expect_true(any(vapply(fitted$x$data, function(x) identical(x$fill, "toself"), logical(1))))
  expect_false(any(vapply(sampled$x$data, function(x) identical(x$fill, "toself"), logical(1))))
  expect_equal(vapply(sampled$x$layout$shapes, `[[`, numeric(1), "y0"), c(-2, 0, 2))
  expect_equal(unlist(lapply(sampled$x$data, `[[`, "y"), use.names = FALSE), out$process_residuals$recruitment$est)
  expect_null(sampled$x$layout$updatemenus)
})

test_that("one matrix provides NFM metrics, shared limits and model-specific M blocks", {
  out <- .process_view_fixture()
  out$deviations$N$age[out$deviations$N$age == 3] <- 7
  p <- plot_process_deviations(out)
  spec <- p$jsHooks$render[[1L]]$data
  expect_equal(spec$models, c("scaa", "sam", "ncam"))
  expect_equal(vapply(spec$metrics, `[[`, "", "name"), c("N", "F", "M"))
  expect_true(spec$signed)
  expect_equal(spec$metrics[[1]]$cells[[2]]$ages, c(2, 7))
  expect_equal(spec$metrics[[3]]$ylab, "M state block")
  expect_equal(spec$metrics[[3]]$linear, c(-2, 2))
  expect_equal(spec$metrics[[3]]$cells[[2]]$ticktext, c("2-4", "5-8"))
  expect_equal(spec$metrics[[3]]$cells[[3]]$ticktext, c("2", "3-5", "6-8"))
  expect_equal(lengths(lapply(spec$metrics[[3]]$cells[-1], `[[`, "linear")), c(2L, 3L))
  expect_equal(vapply(spec$colorscale, `[[`, "", 2L), c(tail(tam_cols(), 1L), "#FFFFFF", tam_cols()[1]))
  expect_equal(sort(unlist(spec$metrics[[1]]$cells[[2]]$linear)), sort(out$deviations$N$est[out$deviations$N$model == "sam"]))
  out$deviations$N <- NULL
  spec <- plot_process_deviations(out)$jsHooks$render[[1L]]$data
  expect_equal(spec$initial_metric, 1L)
})

test_that("bubble frames use supplied values, stable sign traces and one size scale", {
  out <- .process_view_fixture()
  spec <- plot_process_residuals(out)$jsHooks$render[[1L]]$data
  expect_equal(spec$models, c("scaa", "sam", "ncam"))
  expect_equal(c(spec$component, spec$model), c(0L, 1L))
  scales <- numeric()
  for (i in seq_along(spec$components)) {
    component <- spec$components[[i]]
    expect_equal(component$name, c("N", "F", "M")[i])
    expect_equal(component$cells[[1]]$layout$annotations[[1]]$text, "Unavailable for this model component")
    expect_true(all(vapply(component$cells[[1]]$data, function(trace) length(trace$x) == 0L, logical(1))))
    for (j in 2:3) {
      cell <- component$cells[[j]]
      values <- out$process_residuals[[component$name]]
      z <- values[values$model == spec$models[j], ]
      for (k in 1:2) {
        trace <- cell$data[[k]]
        expect_identical(trace$name, c("-", "+")[k])
        expect_equal(as.numeric(trace$marker$size), abs(z$est[if (k == 1) z$est <= 0 else z$est > 0]))
        expect_equal(trace$marker$color, c(tail(tam_cols(), 1), tam_cols()[1])[k])
        scales <- c(scales, trace$marker$sizeref)
      }
      if (component$name == "M") expect_equal(cell$layout$yaxis$ticktext, unique(z$age_block))
    }
  }
  expect_equal(unique(scales), 4 / 900)
})

test_that("bubble component buttons and native model frames coordinate independently", {
  node <- Sys.which("node"); skip_if(!nzchar(node))
  config <- tempfile(fileext = ".json"); on.exit(unlink(config))
  jsonlite::write_json(plot_process_residuals(.process_view_fixture())$jsHooks$render[[1L]]$data,
    config, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
  result <- system2(node, c(shQuote(test_path("process-bubbles-controls.js")),
    shQuote(system.file("js", "process-bubbles.js", package = "tinyAM")), shQuote(config)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "Process bubble controls passed")
})

test_that("dashboard dropdown adaptation retains native panes and handles selection and keyboard", {
  node <- Sys.which("node"); skip_if(!nzchar(node))
  result <- system2(node, c(shQuote(test_path("dashboard-dropdowns-controls.js")),
    shQuote(system.file("js", "dashboard-dropdowns.js", package = "tinyAM"))), stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "Dashboard dropdown controls passed")
})

test_that("the Random template emits an N dropdown only when both N panels exist", {
  template <- readLines(system.file("rmd", "parameters_random.Rmd", package = "tinyAM"))
  headings <- template[startsWith(template, '`r if ("log_n')]
  expect_length(headings, 2L)
  expressions <- lapply(sub("`$", "", substring(headings, 4L)), function(x) parse(text = x))
  render <- function(parameters) {
    env <- list2env(list(random_par = stats::setNames(rep(list(data.frame()), length(parameters)), parameters)))
    unlist(lapply(expressions, eval, envir = env), use.names = FALSE)
  }
  expect_null(render(character()))
  expect_identical(render("log_n"), "### N")
  expect_identical(render("log_n0"), "### N")
  both <- render(c("log_n", "log_n0"))
  expect_true(all(grepl('data-tab-group="N"', both, fixed = TRUE)))
  expect_match(both[1], "States", fixed = TRUE)
  expect_match(both[2], "Initial abundance", fixed = TRUE)
})
