
.tam_cols <- c(
  "#003e56",
  "#007187",
  "#028f8d",
  "#2aaa85",
  "#93c25c",
  "#fec124",
  "#ffa31d",
  "#fe8427",
  "#f05c28"
)

.buttons <- list(
  list(
    type = "buttons",
    y = 1, x = 0, pad = list(r = 75), xanchor = "right", yanchor = "top",
    buttons = list(
      list(method = "relayout", args = list("yaxis.type", "linear"), label = "linear"),
      list(method = "relayout", args = list("yaxis.type", "log"),    label = "log")
    )
  )
)

#' Plot utilities for `tidy_tam()` outputs
#'
#' These functions provide minimal, Plotly-based visualizations for outputs
#' from [`tidy_tam()`] and related TAM model objects:
#'
#' - `plot_trend()`: time‐series trends with optional confidence ribbons.
#' - `plot_heatmap()`: age × year heatmaps for matrix quantities (e.g., `F`, `M`, `N`).
#' - `plot_obs_pred()`: observed vs predicted fits by year.
#' - `plot_bubbles()`: residuals (OSA or standardized) as bubbles by age × year.
#' - `plot_resid()`: residual diagnostics vs year, cohort, age, and expected values.
#'
#' @param data A data frame containing at least `year` and `est` columns.
#' @param ylab Label for the y-axis.
#' @param zlab Label for the colorbar in heatmaps.
#' @param x,xlab For residual plots, specify the x variable and its label.
#' @param title Plot title.
#' @param add_intervals Logical; if `TRUE`, show confidence intervals when `lwr` and `upr` are present.
#' @param add_buttons Logical; if `TRUE`, include buttons to toggle y-axis between linear/log.
#' @param ... Additional arguments to pass to [plotly::plot_ly()] such as `color`, `colors`, etc.
#'
#' @return A Plotly object.
#'
#' @examples
#' ## Fit two models
#' N_dev <- fit_tam(
#'   cod_obs, years = 1983:2024, ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   M_settings = list(process = "off", assumption = ~ I(0.3))
#' )
#' M_dev <- update(
#'   N_dev,
#'   N_settings = list(process = "off", init_N0 = TRUE),
#'   M_settings = list(process = "ar1", assumption = ~ I(0.3),
#'                     age_breaks = seq(2, 14, by = 6))
#' )
#' tabs <- tidy_tam(N_dev, M_dev)
#'
#' ## Trend plots
#' plot_trend(N_dev$pop$biomass, ylab = "Biomass")
#' plot_trend(tabs$pop$biomass, color = ~model, ylab = "Biomass")
#'
#' N_dev$pop$N |>
#'   dplyr::mutate(age = factor(age, levels = sort(unique(age)))) |>
#'   plot_trend(ylab = "N", color = ~age, colors = viridis::viridis(14))
#'
#' ## Age-Year Heatmap
#' plot_heatmap(N_dev$pop$F, zlab = "F")
#' plot_heatmap(M_dev$pop$F, zlab = "F")
#'
#' ## Observed vs Predicted
#' plot_obs_pred(N_dev$obs_pred$catch, frame = ~age, ylab = "Catch")
#' plot_obs_pred(tabs$obs_pred$index,
#'               color = ~model, legendgroup = ~model,
#'               frame = ~age, ylab = "Index")
#'
#' ## Residual Bubbles
#' plot_bubbles(N_dev$obs_pred$catch)
#' plot_bubbles(tabs$obs_pred$index, frame = ~model)
#'
#' ## Residual Scatter
#' plot_resid(N_dev$obs_pred$catch, x = ~age, xlab = "Age")
#' plot_resid(tabs$obs_pred$index, frame = ~model,
#'             x = ~year, xlab = "Year", showlegend = FALSE)
#'
#' @name plot_tam
#' @import plotly
#' @export
plot_trend <- function(
    data,
    ylab = "",
    title = "",
    add_intervals = TRUE,
    add_buttons = TRUE,
    ...
) {
  args <- list(...)
  if (!"color" %in% names(args)) {
    args$color <- I('#1f77b4')
    args$showlegend <- FALSE
  }
  p <- do.call(plot_ly, c(list(data = data, x = ~year), args))

  if (add_intervals && all(c("lwr", "upr") %in% names(data))) {
    p <- p |> add_ribbons(ymin = ~lwr, ymax = ~upr,
                          opacity = 0.3, line = list(width = 0),
                          showlegend = FALSE)
  }

  buttons <- if (add_buttons) .buttons else NULL

  p |>
    add_lines(y = ~est, line = list(width = 2)) |>
    layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = ylab),
      updatemenus = buttons
    )
}

#' @rdname plot_tam
#' @export
plot_heatmap <- function(
    data,
    zlab = "Estimate",
    title = "",
    ...
) {
  args <- list(...)
  if ("frame" %in% names(args) && !is.null(args$frame)) {
    cli::cli_abort(
      "Argument {.arg frame} is not supported for heatmaps.
       Supply a single frame at a time, or split/filter your data by the grouping
       variable and call {.fn plot_heatmap} per group."
    )
  }

  plot_ly(data, x = ~year, y = ~age, z = ~est,
          colorbar = list(title = zlab), ...) |>
    add_heatmap() |>
    layout(title = title,
           xaxis = list(title = "Year"),
           yaxis = list(title = "Age"))
}


#' @rdname plot_tam
#' @export
plot_obs_pred <- function(
    data,
    ylab = "",
    title = "",
    add_buttons = TRUE,
    ...
) {
  args <- list(...)
  if (!"color" %in% names(args)) {
    args$color <- I('#1f77b4')
    args$showlegend <- FALSE
    args$name <- "Predicted"
  }

  p <- do.call(plot_ly, c(list(data = data, x = ~year), args))

  buttons <- if (add_buttons) .buttons else NULL

  p |>
    add_markers(
      y = ~obs,
      color = ~I("darkgrey"),
      name = "Observed",
      showlegend = FALSE,
      opacity = 0.6
    ) |>
    add_lines(y = ~pred) |>
    layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = ylab),
      updatemenus = buttons
    )
}


#' @rdname plot_tam
#' @export
plot_bubbles <- function(data, ...) {
  if ("osa_res" %in% names(data)) {
    title <- "OSA residuals"
    data$res <- data$osa_res
  } else {
    title <- "Standardized residuals"
    data$res <- data$std_res
  }
  data$sign <- ifelse(data$res > 0, "+", "-")
  data$abs_res <- jitter(abs(data$res), factor = 1e-12) # add a small amount of noise, otherwise plotly does not plot values with the same number
  data <- data[!is.na(data$res), ]

  plot_ly(
    data,
    x = ~year,
    y = ~age,
    marker = list(
      size = ~abs_res,
      sizemin = 1,
      sizeref = 0.01,
      sizemode = "area"
    ),
    color = ~sign,
    colors = c("#377EB8", "#E41A1C"),
    text = ~round(res, 2),
    ...
  ) |>
    add_markers() |>
    layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = "Age")
    )
}


#' @rdname plot_tam
#' @export
plot_resid <- function(
    data,
    x = ~year,
    xlab = "Year",
    ...
) {
  if ("osa_res" %in% names(data)) {
    title <- "OSA residuals"
    data$res <- data$osa_res
  } else {
    title <- "Standardized residuals"
    data$res <- data$std_res
  }
  data <- data[!is.na(data$res), ]

  plot_ly(data, x = x, y = ~res, ...) |>
    add_markers(alpha = 0.8) |>
    layout(xaxis = list(title = xlab),
           yaxis = list(title = title))

}

