residual_plot_data <- function() {
  d <- expand.grid(year = 2000:2002, age = 1:2,
    model = factor(c("A", "B"), levels = c("A", "B")))
  d$obs <- d$pred <- 1
  d$std_res <- seq_len(nrow(d)) / 10
  d$osa_res <- NA_real_
  d$osa_res[1:2] <- c(-.4, .7)
  d$osa_method <- "oneStepGeneric"
  d$is_proj <- FALSE
  d
}

residual_config <- function(p) p$jsHooks$render[[1L]]$data
residual_points <- function(cell) sum(vapply(cell$data, function(trace) length(trace$x), integer(1)))

test_that("initial resolution selects retained residuals per model and preserves OSA gaps", {
  d <- residual_plot_data(); before <- d
  auto <- .observation_residual_data(d, "auto")
  expect_identical(auto$status$residual_type, c("osa", "standardized"))
  expect_equal(auto$data$.residual_value[1:6], d$osa_res[1:6])
  expect_equal(auto$data$.residual_value[7:12], d$std_res[7:12])
  expect_equal(auto$status$n, c(2L, 6L))
  expect_identical(auto$status$osa_method[1L], "oneStepGeneric")
  expect_equal(.observation_residual_data(d, "standardized")$data$.residual_value,d$std_res)
  expect_equal(.observation_residual_data(d, "osa")$data$.residual_value,d$osa_res)
  d$obs[1] <- 0; d$is_proj[2] <- TRUE
  expect_true(all(is.na(.observation_residual_data(d, "osa")$data$.residual_value)))
  expect_identical(.observation_residual_data(d, "auto")$status$residual_type,
    c("standardized", "standardized"))
  expect_identical(before, residual_plot_data())
})

test_that("public residual plots retain exact model traces and unavailable OSA views", {
  d <- residual_plot_data(); before <- d
  local_mocked_bindings(.fitted_objective = function(...) stop("Must not compute OSA"))
  for (fun in list(plot_bubbles, plot_resid)) {
    p <- fun(d, frame = ~model)
    expect_s3_class(p, "plotly")
    config <- residual_config(p)
    expect_identical(config$frames, c("A", "B"))
    expect_identical(config$frameLabel, "Model")
    expect_null(p$x$frames)
    expect_equal(residual_points(list(data = plotly::plotly_build(p)$x$data)), 2L)
    expect_identical(config$view, 1L)
    expect_identical(config$defaultViews, c(1L, 0L))
    expect_true(config$showSelector)
    expect_identical(vapply(config$views, `[[`, character(1), "type"), c("standardized", "osa"))
    expect_equal(vapply(config$views[[1]]$cells, residual_points, integer(1)), c(6L, 6L))
    expect_equal(vapply(config$views[[2]]$cells, residual_points, integer(1)), c(2L, 0L))
    note <- config$views[[2]]$cells[[2]]$layout$annotations[[1]]$text
    expect_match(note, "OSA residuals unavailable", fixed = TRUE)
    osa_cell <- config$views[[2]]$cells[[1]]
    expect_length(osa_cell$layout$annotations, 0L)
    expect_identical(osa_cell$status$osa_method, "oneStepGeneric")
    expect_match(paste(unlist(lapply(osa_cell$data, `[[`, "text")), collapse = " "),
      "oneStepGeneric", fixed = TRUE)
    expect_identical(config$combineQuantile, c(TRUE, TRUE))
    expect_length(config$views[[2]]$cells[[2]]$data, 0L)
    expect_identical(residual_config(fun(d, residual_type = "osa"))$view, 1L)
  }
  expect_identical(d, before)
  expect_error(plot_resid(d, residual_type = "unknown"), "arg")
})

test_that("subsets choose one residual type across the whole panel", {
  d <- residual_plot_data()
  p <- plot_resid(d, frame = ~age, color = ~model, colors = c("#001122", "#ee8800"))
  config <- residual_config(p)
  expect_identical(config$frameLabel, "Age")
  status <- config$views[[1]]$cells[[2]]$status
  expect_identical(status$residual_type, c("standardized", "standardized"))
  expect_equal(status$n, c(3L, 3L))
  expect_identical(config$defaultViews, c(0L, 0L))
  expect_equal(residual_points(config$views[[1]]$cells[[2]]), 6L)
  p <- plot_resid(d, frame = ~model, color = ~model, colors = c("#001122", "#ee8800"))
  config <- residual_config(p)
  for (i in 1:2) {
    nonempty <- Filter(function(trace) length(trace$x), config$views[[1]]$cells[[i]]$data)
    expect_identical(nonempty[[1]]$marker$color, c("#001122", "#ee8800")[i])
  }
})

test_that("residual views preserve custom hover, scalar inputs and process bubbles", {
  d <- residual_plot_data()[1L, ]
  d$osa_method <- NULL
  p <- plot_resid(d, text = ~paste("custom", year), hovertemplate = "Custom: %{text}",
    residual_type = "standardized")
  expect_length(residual_config(p)$frames, 1L)
  expect_identical(residual_config(p)$view, 0L)
  traces <- Filter(function(trace) length(trace$x), residual_config(p)$views[[1]]$cells[[1]]$data)
  expect_identical(traces[[1]]$hovertemplate, "Custom: %{text}")
  expect_match(traces[[1]]$text, "custom", fixed = TRUE)
  expect_true(residual_config(plot_bubbles(d))$compactHeader)
  titled <- residual_config(plot_bubbles(d, title = "Catch residuals"))
  expect_false(titled$compactHeader)
  expect_identical(titled$views[[2]]$cells[[1]]$layout$title, "Catch residuals")
  expect_identical(residual_config(plot_bubbles(d))$views[[2]]$cells[[1]]$status$osa_method,
    "method unspecified")
  d$osa_res <- NA_real_
  expect_identical(residual_config(plot_resid(d))$views[[1]]$cells[[1]]$status$residual_type, "standardized")
  expect_false(residual_config(plot_resid(d))$showSelector)
  d$std_res <- NA_real_
  expect_equal(residual_points(residual_config(plot_resid(d))$views[[1]]$cells[[1]]), 0L)
  process <- data.frame(year = 2000:2001, age_block = c("2-4", "5-8"), est = c(-1, 2))
  p <- plot_bubbles(process, res = ~est, y = ~age_block, title = "Process")
  expect_s3_class(p, "plotly")
  expect_null(p$jsHooks$render)
  expect_equal(sum(vapply(plotly::plotly_build(p)$x$data, function(trace) length(trace$x), integer(1))), 2L)
  expect_true(all(vapply(plotly::plotly_build(p)$x$data, function(trace) identical(trace$marker$line$width, 0), logical(1))))
})

test_that("native residual selectors replace complete cells in independent widgets", {
  node <- Sys.which("node"); skip_if(!nzchar(node))
  config <- tempfile(fileext = ".json"); on.exit(unlink(config))
  d <- residual_plot_data()
  osa <- residual_config(plot_bubbles(d, frame = ~model))
  d$osa_res <- NA_real_
  d$family <- ifelse(d$model == "B", "pld_lognormal", "lognormal")
  d$obs[7L] <- 0
  d$quantile_res <- d$std_res + .5
  pld <- residual_config(plot_bubbles(d, frame = ~model))
  jsonlite::write_json(list(osa = osa, pld = pld),
    config, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
  result <- system2(node, c(shQuote(test_path("residual-selection-controls.js")),
    shQuote(system.file("js", "residual-selection.js", package = "tinyAM")), shQuote(config)),
    stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "Residual controls passed", fixed = TRUE)
})


test_that("overlay defaults require eligible OSA in every displayed model", {
  d <- residual_plot_data()
  p <- plot_resid(d)
  config <- residual_config(p)
  expect_identical(config$view, 0L)
  expect_equal(residual_points(list(data=p$x$data)), 12L)
  d$osa_res[7L] <- -.5
  p <- plot_resid(d, frame=~age)
  expect_identical(residual_config(p)$defaultViews, c(1L, 1L))
  expect_equal(residual_points(residual_config(p)$views[[2]]$cells[[2]]), 0L)
  traces <- do.call(c, lapply(residual_config(plot_bubbles(d))$views, function(view)
    do.call(c, lapply(view$cells, `[[`, "data"))))
  expect_true(all(vapply(traces, function(trace) identical(trace$marker$line$width, 0), logical(1))))
})

test_that("Gamma residuals select quantiles and only PLD includes sampled zeros", {
  d <- residual_plot_data()
  d$osa_res <- NA_real_
  d$family <- ifelse(d$model == "B", "gamma", "lognormal")
  d$quantile_res <- d$std_res + 1
  d$is_structural <- FALSE
  d$obs[7L] <- 0
  config <- residual_config(plot_bubbles(d, frame = ~model))
  expect_identical(config$defaultViews, c(0L, 2L))
  expect_identical(config$combineQuantile, c(TRUE, FALSE))
  expect_equal(residual_points(config$views[[3L]]$cells[[2L]]), 5L)
  d$family[d$model == "B"] <- "pld_gengamma"
  d$is_observed <- TRUE
  d$p_zero <- 0.2
  config <- residual_config(plot_bubbles(d, frame = ~model))
  expect_equal(residual_points(config$views[[3L]]$cells[[2L]]), 6L)
  frequency <- .zero_frequency_data(d)
  expect_true(all(frequency$model == "B"))
  expect_equal(sum(frequency$n), 6L)
})

test_that("explicit hover vectors follow model frames and partial OSA rows", {
  d <- residual_plot_data()
  labels <- paste0("row-", seq_len(nrow(d)))
  p <- plot_resid(d, frame = ~model, text = labels,
    customdata = seq_len(nrow(d)), hovertemplate = "%{text}: %{customdata}")
  osa <- residual_config(p)$views[[2]]$cells[[1]]
  traces <- Filter(function(trace) length(trace$x), osa$data)
  expect_identical(as.character(traces[[1]]$text), labels[1:2])
  expect_equal(as.numeric(traces[[1]]$customdata), 1:2)
  expect_true(all(as.character(traces[[1]]$hovertemplate) == "%{text}: %{customdata}"))
  standard <- residual_config(p)$views[[1]]$cells[[2]]
  traces <- Filter(function(trace) length(trace$x), standard$data)
  expect_identical(as.character(traces[[1]]$text), labels[7:12])
  expect_equal(as.numeric(traces[[1]]$customdata), 7:12)
})
