
## Helpers ----

`%||%` <- function(x, y) if (is.null(x)) y else x
plotly_traces <- function(p) {
  # Materialize traces so we can inspect type/mode/fill/etc.
  plotly::plotly_build(p)$x$data
}
trace_types <- function(p) {
  trs <- plotly_traces(p)
  vapply(trs, function(tr) tr[["type"]] %||% "", character(1))
}
trace_modes <- function(p) {
  trs <- plotly_traces(p)
  vapply(trs, function(tr) tr[["mode"]] %||% "", character(1))
}
trace_fills <- function(p) {
  trs <- plotly_traces(p)
  vapply(trs, function(tr) tr[["fill"]] %||% "", character(1))
}


## Globals ----

N_dev <- fit_tam(
  cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  M_settings = list(process = "off", assumption = ~ I(0.3)),
  proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1),
  silent = TRUE
)
M_dev <- update(
  N_dev,
  N_settings = list(process = "off", init_N0 = TRUE),
  M_settings = list(process = "ar1", assumption = ~ I(0.3),
                    age_breaks = c(2, 14))
)

tabs <- tidy_tam(N_dev, M_dev)


## plot_trend ----

test_that("plot_trend produces a valid plotly object", {
  p1 <- plot_trend(N_dev$pop$biomass)
  p2 <- plot_trend(tabs$pop$biomass, color = ~model, ylab = "Biomass")

  expect_s3_class(p1, "plotly")
  expect_s3_class(p2, "plotly")

  # After build, there should be a scatter trace with mode 'lines'
  expect_true(any(grepl("lines", trace_modes(p1), fixed = TRUE)))
  expect_true(any(grepl("lines", trace_modes(p2), fixed = TRUE)))
})

test_that("plot_trend handles confidence intervals and log buttons", {
  p <- plot_trend(N_dev$pop$biomass, add_intervals = TRUE, add_buttons = TRUE)
  expect_s3_class(p, "plotly")

  # ribbons materialize as scatter traces with fill='toself'
  expect_true(any(trace_fills(p) == "toself"))

  # updatemenus should include yaxis.type toggle
  um <- p$x$layoutAttrs[[1]]$updatemenus[[1]]$buttons
  btn_args <- unlist(lapply(um, `[[`, "args"), recursive = TRUE, use.names = FALSE)
  expect_true(any(grepl("yaxis.type", as.character(btn_args), fixed = TRUE)))
})

## plot_heatmap ----

test_that("plot_heatmap works for single dataset", {
  p <- plot_heatmap(N_dev$pop$F, zlab = "F")
  expect_s3_class(p, "plotly")
  expect_true(any(trace_types(p) == "heatmap"))
})

test_that("plot_heatmap aborts when frame is supplied", {
  expect_error(
    plot_heatmap(tabs$pop$F, zlab = "F", frame = ~model),
    class = "rlang_error"
  )
})

## plot_obs_pred ----

test_that("plot_obs_pred includes both observed (markers) and predicted (lines)", {
  p <- plot_obs_pred(N_dev$obs_pred$catch, frame = ~age)
  expect_s3_class(p, "plotly")

  modes <- suppressWarnings(trace_modes(p))
  expect_true(any(grepl("markers", modes, fixed = TRUE)))
  expect_true(any(grepl("lines",   modes, fixed = TRUE)))
})

## plot_bubbles ----

test_that("plot_bubbles produces marker traces", {
  p <- plot_bubbles(N_dev$obs_pred$catch)
  expect_s3_class(p, "plotly")
  expect_true(any(grepl("markers", trace_modes(p), fixed = TRUE)))
})

test_that("plot_bubbles supports model comparisons via frame", {
  p <- plot_bubbles(tabs$obs_pred$index, frame = ~model)
  expect_s3_class(p, "plotly")
  # still at least one marker trace after build
  expect_true(any(grepl("markers", trace_modes(p), fixed = TRUE)))
})

## plot_resid ----

test_that("plot_resid works with default and custom x", {
  p1 <- plot_resid(N_dev$obs_pred$catch)
  p2 <- plot_resid(N_dev$obs_pred$catch, x = ~age, xlab = "Age")

  expect_s3_class(p1, "plotly")
  expect_s3_class(p2, "plotly")

  # after build, there should be marker traces
  expect_true(any(grepl("markers", trace_modes(p1), fixed = TRUE)))
  expect_true(any(grepl("markers", trace_modes(p2), fixed = TRUE)))
})
