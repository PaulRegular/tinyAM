
.log_linear_buttons <- list(
  list(
    type = "buttons",
    y = 1, x = 0, pad = list(r = 60), xanchor = "right", yanchor = "top",
    buttons = list(
      list(method = "relayout", args = list("yaxis.type", "linear"), label = "Linear"),
      list(method = "relayout", args = list("yaxis.type", "log"),    label = "Log")
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
#' @param data A data frame of observations or estimates.
#' @param x Formula specifying x variable
#' @param y Formula specifying y variable
#' @param ylab Label for the y-axis.
#' @param xlab Label for the x-axis.
#' @param zlab Label for the colorbar in heatmaps.
#' @param title Plot title.
#' @param add_intervals Logical; if `TRUE`, show confidence intervals when `lwr` and `upr` are present.
#' @param add_buttons Logical; if `TRUE`, include buttons to toggle y-axis between Linear / Log and, if
#'   applicable, Show CI / Hide CI.
#' @param ... Additional arguments to pass to [plotly::plot_ly()] such as `color`, `colors`, etc.
#'
#' @return A Plotly object.
#'
#' @example inst/examples/example_fits.R
#' @examples
#' ## Stack outputs ----
#' tabs <- tidy_tam(model_list = fits)
#'
#' ## Trend plots ----
#' plot_trend(N_dev$pop$biomass, ylab = "Biomass")
#' plot_trend(tabs$pop$biomass, color = ~model, ylab = "Biomass")
#'
#' N_dev$pop$N |>
#'   subset(age %in% 2:8) |>
#'   plot_trend(ylab = "N", color = ~as.character(age))
#'
#' ## Age-Year Heatmap ----
#' plot_heatmap(N_dev$pop$F, zlab = "F")
#' plot_heatmap(M_dev$pop$F, zlab = "F")
#'
#' ## Observed vs Predicted ----
#' plot_obs_pred(N_dev$obs_pred$catch, frame = ~age, ylab = "Catch")
#' plot_obs_pred(tabs$obs_pred$index,
#'               color = ~model, legendgroup = ~model,
#'               frame = ~age, ylab = "Index")
#'
#' ## Residual Bubbles ----
#' plot_bubbles(N_dev$obs_pred$catch)
#' plot_bubbles(tabs$obs_pred$index, frame = ~model)
#'
#' ## Residual Scatter ----
#' plot_resid(N_dev$obs_pred$catch, x = ~age, xlab = "Age")
#' plot_resid(tabs$obs_pred$index, frame = ~model,
#'             x = ~year, xlab = "Year", showlegend = FALSE)
#'
#' ## Parameter Estimates ----
#' plot_par(N_dev$fixed_par)
#' plot_par(tabs$fixed_par, color = ~model)
#'
#' @name plot_tam
#' @import plotly
#' @export
plot_trend <- function(
    data,
    x = ~year,
    y = ~est,
    xlab = "Year",
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

  p <- do.call(plotly::plot_ly, c(list(data = data, x = x), args))
  max_y <- max(model.frame(y, data = data), na.rm = TRUE) * 1.05

  has_ci <- add_intervals && all(c("lwr", "upr") %in% names(data))
  if (has_ci) {
    max_y <- max(data$upr, na.rm = TRUE) * 1.05
    p <- p |>
      plotly::add_ribbons(
        ymin = ~lwr, ymax = ~upr,
        opacity = 0.3, line = list(width = 0),
        showlegend = FALSE
      )
  }

  shapes <- NULL
  if ("is_proj" %in% names(data) && any(data$is_proj)) {
    max_obs_year <- suppressWarnings(max(data$year[!data$is_proj], na.rm = TRUE))
    if (is.finite(max_obs_year)) {
      shapes <- list(list(
        type  = "line",
        x0    = max_obs_year, x1 = max_obs_year,
        y0    = 0, y1 = 1, yref = "paper",
        line  = list(dash = "dot", width = 1, color = "black")
      ))
    }
  }

  p <- p |>
    plotly::add_lines(y = y) |>
    plotly::layout(
      title = title,
      xaxis = list(title = xlab),
      yaxis = list(title = ylab, range = c(0, max_y)),
      shapes = shapes
    )

  if (add_buttons) {
    pb <- plotly::plotly_build(p)
    data_traces <- pb$x$data
    if (is.null(data_traces)) {
      data_traces <- list()
    }

    # find ribbon traces (CI)
    rib_idx <- which(vapply(
      data_traces,
      function(tr) identical(tr$type, "scatter") && !is.null(tr$fill) && nzchar(tr$fill),
      logical(1)
    ))

    menus <- .log_linear_buttons

    if (length(rib_idx)) {
      idx0 <- as.integer(rib_idx - 1L)
      ci_menu <- list(
        type = "buttons",
        y = 1, x = 0, pad = list(r = 60, t = 75), xanchor = "right", yanchor = "top",
        buttons = list(
          list(
            label = "Show CI",
            method = "restyle",
            args = list(list(visible = rep(TRUE, length(idx0))), idx0)
          ),
          list(
            label = "Hide CI",
            method = "restyle",
            args = list(list(visible = rep(FALSE, length(idx0))), idx0)
          )
        )
      )
      menus <- c(menus, list(ci_menu))
    }

    if (is.null(p$x$layout)) {
      p$x$layout <- list()
    }
    p$x$layout$updatemenus <- menus
    if (!is.null(pb$x$frames)) {
      p$x$frames <- pb$x$frames
    }
  }

  p
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
  max_y <- max(c(data$obs, data$pred), na.rm = TRUE) * 1.05

  shapes <- NULL
  if ("is_proj" %in% names(data) && any(data$is_proj)) {
    max_obs_year <- suppressWarnings(max(data$year[!data$is_proj], na.rm = TRUE))
    if (is.finite(max_obs_year)) {
      shapes <- list(list(
        type  = "line",
        x0    = max_obs_year,
        x1    = max_obs_year,
        y0    = 0,
        y1    = 1,
        yref  = "paper",
        line  = list(dash = "dot", width = 1, color = "black")
      ))
    }
  }

  p <- do.call(plot_ly, c(list(data = data, x = ~year), args))

  buttons <- if (add_buttons) .log_linear_buttons else NULL

  p |>
    add_markers(
      y = ~obs,
      color = ~I("grey40"),
      name = "Observed",
      showlegend = FALSE,
      opacity = 0.5
    ) |>
    add_lines(y = ~pred) |>
    layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = ylab, range = c(0, max_y)),
      updatemenus = buttons,
      shapes = shapes
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
  data$abs_res <- jitter(abs(data$res), factor = 1e-12) # tiny jitter to avoid identical-size rendering issues
  data <- data[!is.na(data$res), ]

  args <- list(...)
  if (is.null(args$colors)) {
    args$colors <- c("#377EB8", "#E41A1C")
  }

  p <- do.call(
    plotly::plot_ly,
    c(
      list(
        data  = data,
        x     = ~year,
        y     = ~age,
        marker = list(
          size     = ~abs_res,
          sizemin  = 1,
          sizeref  = 0.01,
          sizemode = "area"
        ),
        color = ~sign,
        text  = ~round(res, 2)
      ),
      args
    )
  )

  p |>
    plotly::add_markers() |>
    plotly::layout(
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

#' @rdname plot_tam
#' @export
plot_par <- function(data, ...) {
  data$coef <- factor(ifelse(is.na(data$coef), data$par, paste0(data$par, ": ", data$coef)))
  max_y <- max(data$upr[grepl("q|sd", data$par)]) * 1.05

  args <- list(...)
  if (is.null(args$color)) {
    args$color <- I("#377EB8")
  }
  p <- do.call(plot_ly, c(list(data = data, x = ~est, y = ~coef), args))

  p |>
    add_segments(
      x = ~lwr, xend = ~upr,
      y = ~coef, yend = ~coef,
      showlegend = FALSE
    ) |>
    add_markers(showlegend = TRUE) |>
    layout(
      xaxis = list(
        title = "Estimate",
        range = c(0, max_y)
      ),
      yaxis = list(
        title = ""
      )
    )
}
