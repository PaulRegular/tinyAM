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

test_that("Auto selects retained residuals per model and preserves OSA gaps", {
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

test_that("public residual plots retain models, labels and unavailable OSA frames", {
  d <- residual_plot_data(); before <- d
  local_mocked_bindings(.fitted_objective = function(...) stop("Must not compute OSA"))
  for (fun in list(plot_bubbles, plot_resid)) {
    p <- fun(d, frame = ~model)
    expect_s3_class(p, "plotly")
    config <- residual_config(p)
    expect_identical(config$frames, c("A", "B"))
    expect_identical(vapply(p$x$frames, `[[`, character(1), "name"), c("A", "B"))
    expect_equal(vapply(p$x$frames, residual_points, integer(1)), c(2L, 6L))
    expect_true(all(vapply(p$x$frames, function(frame) length(frame$traces) > 0L, logical(1))))
    expect_identical(plotly::plotly_build(p)$x$frames, p$x$frames)
    expect_identical(config$view, 0L)
    expect_identical(vapply(config$views, `[[`, character(1), "type"), c("auto", "standardized", "osa"))
    expect_equal(vapply(config$views[[1]]$cells, residual_points, integer(1)), c(2L, 6L))
    expect_equal(vapply(config$views[[2]]$cells, residual_points, integer(1)), c(6L, 6L))
    expect_equal(vapply(config$views[[3]]$cells, residual_points, integer(1)), c(2L, 0L))
    note <- config$views[[3]]$cells[[2]]$layout$annotations[[1]]$text
    expect_match(note, "OSA residuals unavailable", fixed = TRUE)
    expect_match(config$views[[1]]$cells[[1]]$layout$annotations[[1]]$text,
      "oneStepGeneric", fixed = TRUE)
    slots <- unlist(lapply(config$views, function(view) vapply(view$cells, function(cell) length(cell$data), integer(1))))
    expect_length(unique(slots), 1L)
    expect_identical(residual_config(fun(d, residual_type = "osa"))$view, 2L)
  }
  expect_identical(d, before)
  expect_error(plot_resid(d, residual_type = "unknown"), "arg")
})

test_that("Auto is chosen for the whole panel rather than separately for animation ages", {
  d <- residual_plot_data()
  p <- plot_resid(d, frame = ~age, color = ~model, colors = c("#001122", "#ee8800"))
  config <- residual_config(p)
  status <- config$views[[1]]$cells[[2]]$status
  expect_identical(status$residual_type, c("osa", "standardized"))
  expect_equal(status$n, c(0L, 3L))
  expect_equal(residual_points(config$views[[1]]$cells[[2]]), 3L)
  p <- plot_resid(d, frame = ~model, color = ~model, colors = c("#001122", "#ee8800"))
  config <- residual_config(p)
  for (i in 1:2) {
    nonempty <- Filter(function(trace) length(trace$x), config$views[[2]]$cells[[i]]$data)
    expect_identical(nonempty[[1]]$marker$color, c("#001122", "#ee8800")[i])
  }
})

test_that("residual views preserve custom hover, scalar inputs and process bubbles", {
  d <- residual_plot_data()[1L, ]
  d$osa_method <- NULL
  p <- plot_resid(d, text = ~paste("custom", year), hovertemplate = "Custom: %{text}",
    residual_type = "standardized")
  expect_false(residual_config(p)$animated)
  expect_identical(residual_config(p)$view, 1L)
  traces <- Filter(function(trace) length(trace$x), residual_config(p)$views[[2]]$cells[[1]]$data)
  expect_identical(traces[[1]]$hovertemplate, "Custom: %{text}")
  expect_match(traces[[1]]$text, "custom", fixed = TRUE)
  expect_match(residual_config(plot_bubbles(d))$views[[1]]$cells[[1]]$layout$annotations[[1]]$text,
    "method unspecified", fixed = TRUE)
  d$osa_res <- NA_real_
  expect_identical(residual_config(plot_resid(d))$views[[1]]$cells[[1]]$status$residual_type, "standardized")
  d$std_res <- NA_real_
  expect_equal(residual_points(residual_config(plot_resid(d))$views[[1]]$cells[[1]]), 0L)
  process <- data.frame(year = 2000:2001, age_block = c("2-4", "5-8"), est = c(-1, 2))
  p <- plot_bubbles(process, res = ~est, y = ~age_block, title = "Process")
  expect_s3_class(p, "plotly")
  expect_null(p$jsHooks$render)
  expect_equal(sum(vapply(plotly::plotly_build(p)$x$data, function(trace) length(trace$x), integer(1))), 2L)
})

test_that("native residual dropdown and frames coordinate in independent widgets", {
  node <- Sys.which("node"); skip_if(!nzchar(node))
  config <- tempfile(fileext = ".json"); on.exit(unlink(config))
  jsonlite::write_json(residual_config(plot_bubbles(residual_plot_data(), frame = ~model)),
    config, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
  result <- system2(node, c(shQuote(test_path("residual-selection-controls.js")),
    shQuote(system.file("js", "residual-selection.js", package = "tinyAM")), shQuote(config)),
    stdout = TRUE, stderr = TRUE)
  expect_null(attr(result, "status"))
  expect_match(paste(result, collapse = "\n"), "Residual controls passed", fixed = TRUE)
})
