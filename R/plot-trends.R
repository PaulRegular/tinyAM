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
#' - `plot_par()`: parameter estimates and confidence limits, including signed effects.
#' - `add_plot_buttons()`: finish Linear / Log and CI controls after adding layers.
#'
#' @details
#' These functions display the supplied tables without refitting the model or
#' recalculating uncertainty. Trend ribbons use `lwr` and `upr` confidence limits;
#' they are not prediction intervals for future observations. Catch and survey
#' predictions are conditional medians. Residuals come from [tidy_obs_pred()].
#' Subset to a single survey or supply a grouping aesthetic when combining surveys
#' so distinct series are not connected. A log axis requires positive values.
#' Both trend and observed/predicted plots accept an `x` formula and `xlab`;
#' use `x = ~ age, xlab = "Age", frame = ~ year` for annual age profiles.
#' The historical/projection separator is shown only on a year axis.
#' Tables with `fold` (forecast origin) and `horizon` columns add forecast
#' status to default hover text. Different origins suppress a shared separator.
#' Explicit `text` or `hovertemplate` arguments override these hover defaults.
#' [plot_age_year()] assembles both orientations with view-switching buttons.
#' For custom overlays, call `plot_trend(add_buttons = FALSE)`, add the layers
#' with Plotly, then call `add_plot_buttons()` so CI controls include all ribbons.
#' See [tinyAM-model] for definitions and units of population quantities.
#' Optional `unit` columns supply default y-axis labels and hover units for
#' natural estimates; an explicit `ylab` or `hovertemplate` takes precedence.
#' Known inconsistent units must be converted explicitly
#' before overlaying models. Single-coordinate trends show a point and a
#' vertical confidence interval instead of an invisible line segment.
#' Missing intervals split confidence ribbons into separate runs. Signed values
#' and intervals expand the default axis below zero; nonnegative trends retain
#' a zero baseline.
#'
#' @param data A data frame of observations or estimates.
#' @param x Formula specifying x variable
#' @param y Formula specifying y variable
#' @param res Optional one-column formula for signed bubble values, e.g.
#'   `~est` for process deviations. Omitted values use OSA or standardized
#'   observation residuals. Use `y = ~age_block` for shared M state blocks.
#' @param residual_type Initial observation residual view: `"auto"` (default),
#'   `"standardized"` or `"osa"`. The default uses retained OSA residuals
#'   when available, otherwise standardized residuals. Model frames choose
#'   separately; overlays use OSA only when all models have them. A dropdown
#'   offers OSA and Standardized when OSAs are available. Partial gaps stay missing;
#'   unavailable OSA values are identified without computing them. A supplied
#'   `res` formula in `plot_bubbles()` bypasses these observation controls.
#' @param ylab Label for the y-axis.
#' @param xlab Label for the x-axis.
#' @param zlab Label for the colorbar in heatmaps.
#' @param title Plot title.
#' @param p A Plotly object, with all desired layers already added.
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
#'   plot_trend(ylab = "N", color = ~ as.character(age))
#'
#' ## Age-Year Heatmap ----
#' plot_heatmap(N_dev$pop$F, zlab = "F")
#' plot_heatmap(M_dev$pop$F, zlab = "F")
#'
#' ## Observed vs Predicted ----
#' plot_obs_pred(N_dev$obs_pred$catch, frame = ~age, ylab = "Catch")
#' plot_obs_pred(tabs$obs_pred$index,
#'   color = ~model, legendgroup = ~model,
#'   frame = ~age, ylab = "Index"
#' )
#'
#' ## Residual Bubbles ----
#' plot_bubbles(N_dev$obs_pred$catch)
#' plot_bubbles(tabs$obs_pred$index, frame = ~model)
#'
#' ## Residual Scatter ----
#' plot_resid(N_dev$obs_pred$catch, x = ~age, xlab = "Age")
#' plot_resid(tabs$obs_pred$index,
#'   frame = ~model,
#'   x = ~year, xlab = "Year", showlegend = FALSE
#' )
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
  unit <- if (identical(y[[2L]], quote(est))) .unit_label(data, "") else ""
  unit_label <- if (nzchar(unit)) paste0("Estimate (", unit, ")") else ""
  if (missing(ylab)) ylab <- unit_label
  if (!is.data.frame(data) || !nrow(data)) {
    return(plotly::plot_ly(x = numeric(), y = numeric(), type = "scatter", mode = "lines"))
  }
  args <- .forecast_plot_args(data, list(...))
  if (nzchar(unit) && is.null(args$hovertemplate)) {
    args$hovertemplate <- paste0("%{x}<br>%{y} ", unit,
      if (!is.null(args$text)) "<br>%{text}" else "", "<extra>%{fullData.name}</extra>")
  }
  if (!"color" %in% names(args)) {
    args$color <- I("#1f77b4")
    args$showlegend <- FALSE
  }
  if (inherits(args$color, "formula") && is.null(args$legendgroup)) args$legendgroup <- args$color

  values <- unlist(model.frame(y, data = data), use.names = FALSE)
  if (!any(is.finite(values))) {
    return(plotly::plot_ly(x = numeric(), y = numeric(), type = "scatter", mode = "lines"))
  }
  limits <- values[is.finite(values)]
  if (!is.null(args$frame) && identical(all.vars(y), "est")) {
    # Retain empty model traces in every animation frame.
    data$est[is.na(data$est)] <- Inf
  }
  p <- do.call(plotly::plot_ly, c(list(data = data, x = x), args))

  has_ci <- add_intervals && all(c("lwr", "upr") %in% names(data)) && any(is.finite(data$upr))
  single_x <- length(unique(unlist(stats::model.frame(x, data = data)))) == 1L
  if (has_ci) {
    limits <- c(limits, data$lwr[is.finite(data$lwr)], data$upr[is.finite(data$upr)])
    if (!single_x) {
      # Split before making polygons: NA vertices at the closing edge create
      # diagonal fills, while dropping interior gaps would invent intervals.
      groups <- unique(unlist(lapply(args[intersect(names(args),
        c("color", "split", "group", "legendgroup", "frame"))], function(z) {
          if (inherits(z, "formula")) all.vars(z) else character()
        })))
      groups <- intersect(groups, names(data))
      coordinate <- all.vars(x)
      ci <- data[do.call(order, data[c(groups, coordinate)]), , drop = FALSE]
      key <- if (length(groups)) interaction(ci[groups], drop = TRUE) else rep(1L, nrow(ci))
      valid <- is.finite(ci$lwr) & is.finite(ci$upr)
      frame_columns <- if (inherits(args$frame, "formula")) all.vars(args$frame) else character()
      series <- setdiff(groups, frame_columns)
      series <- if (length(series)) interaction(ci[series], drop = TRUE) else rep(1L, nrow(ci))
      runs <- stats::ave(valid, key, FUN = function(v) cumsum(v & !c(FALSE, utils::head(v, -1L))))
      ci$.ci_run <- interaction(series, runs, drop = TRUE)
      ci <- ci[valid, , drop = FALSE]
      if (nrow(ci) && length(frame_columns)) {
        series_columns <- setdiff(groups, frame_columns)
        slots <- merge(unique(data[c(series_columns, frame_columns)]),
          unique(ci[c(series_columns, ".ci_run")]), by = series_columns, sort = FALSE)
        key_columns <- c(frame_columns, ".ci_run")
        row_key <- function(x) do.call(paste, c(x[key_columns], sep = "\r"))
        absent <- slots[!row_key(slots) %in% row_key(ci), , drop = FALSE]
        if (nrow(absent)) {
          # Preserve empty ribbon traces across frames without drawing intervals.
          empty <- ci[match(absent$.ci_run, ci$.ci_run), , drop = FALSE]
          empty[frame_columns] <- absent[frame_columns]
          empty$lwr <- empty$upr <- Inf
          ci <- rbind(ci, empty)
        }
      }
      if (nrow(ci)) p <- p |> plotly::add_ribbons(
        data = ci, split = ~.ci_run,
        ymin = ~lwr, ymax = ~upr,
        opacity = 0.3, line = list(width = 0),
        showlegend = FALSE, hoverinfo = "skip"
      )
    }
    if (single_x) p <- p |> plotly::add_segments(
      y = ~lwr, yend = ~upr, xend = x, showlegend = FALSE,
      meta = list(tam_ci = TRUE), hoverinfo = "skip"
    )
  }

  shapes <- NULL
  if (identical(all.vars(x), "year") && "is_proj" %in% names(data) && any(data$is_proj)) {
    max_obs_year <- .projection_boundary(data)
    if (is.finite(max_obs_year)) {
      shapes <- list(list(
        type = "line",
        x0 = max_obs_year, x1 = max_obs_year,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(dash = "dot", width = 1, color = "black")
      ))
    }
  }

  limits <- range(c(0, limits))
  padding <- max(diff(limits), abs(limits), 1e-8) * .05
  limits <- c(if (limits[1] < 0) limits[1] - padding else 0, limits[2] + padding)
  p <- if (single_x) plotly::add_markers(p, data = data, y = y) else plotly::add_lines(p, data = data, y = y)
  p <- p |> plotly::layout(
      title = title,
      xaxis = list(title = xlab),
      yaxis = list(title = ylab, range = limits),
      shapes = shapes
    )

  if (add_buttons) p <- add_plot_buttons(p)
  p
}

#' @name plot_tam
#' @rdname plot_tam
#' @export
plot_obs_pred <- function(
  data,
  ylab = "",
  title = "",
  add_buttons = TRUE,
  x = ~year,
  xlab = "Year",
  ...
) {
  unit <- .unit_label(data, "")
  unit_label <- if (nzchar(unit)) paste0("Observation (", unit, ")") else ""
  if (missing(ylab)) ylab <- unit_label
  if (!is.data.frame(data) || !nrow(data) ||
    !any(is.finite(c(data$obs, data$pred)))) {
    return(plotly::plot_ly(x = numeric(), y = numeric(), type = "scatter", mode = "lines"))
  }
  args <- .forecast_plot_args(data, list(...))
  if (nzchar(unit) && is.null(args$hovertemplate)) {
    args$hovertemplate <- paste0("%{x}<br>%{y} ", unit,
      if (!is.null(args$text)) "<br>%{text}" else "", "<extra>%{fullData.name}</extra>")
  }
  if (!"color" %in% names(args)) {
    args$color <- I("#1f77b4")
    args$showlegend <- FALSE
    args$name <- "Predicted"
  }
  max_y <- max(c(data$obs, data$pred), na.rm = TRUE) * 1.05

  if (!is.null(args$frame)) {
    # Keep empty series in every frame so Plotly retains the same trace indices.
    # Inf becomes a non-rendered point; the original observations stay unchanged.
    data$obs[is.na(data$obs)] <- Inf
    data$pred[is.na(data$pred)] <- Inf
  }

  shapes <- NULL
  if (identical(all.vars(x), "year") && "is_proj" %in% names(data) && any(data$is_proj)) {
    max_obs_year <- .projection_boundary(data)
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

  p <- do.call(plot_ly, c(list(data = data, x = x), args))

  buttons <- if (add_buttons) .log_linear_buttons else NULL

  p |>
    add_markers(
      y = ~obs,
      color = ~ I("grey40"),
      name = "Observed",
      showlegend = FALSE,
      opacity = 0.5
    ) |>
    add_lines(y = ~pred) |>
    layout(
      title = title,
      xaxis = list(title = xlab),
      yaxis = list(title = ylab, range = c(0, max_y)),
      updatemenus = buttons,
      shapes = shapes
    )
}


#' Switch between annual trends and age profiles
#'
#' @description
#' Display an age-year table with Annual trends and Age profiles buttons, using
#' the same plots and controls as the dashboard. Each view retains its selected
#' animation frame, axis scale, confidence intervals and legend selections.
#'
#' @param data Data frame containing numeric year and age columns, plus the
#'   estimates or observations required by the plotting function.
#' @param plot Function used to draw each view. Defaults to [plot_trend()]; use
#'   [plot_obs_pred()] for observations and predictions. A custom function must
#'   accept `data`, `x`, `xlab`, `frame` and additional plotting arguments, and
#'   return a Plotly object. Add custom layers inside that function.
#' @inheritParams plot_tam
#' @param frame Formula selecting the initial view: `~ age` for Annual trends
#'   (the default), or `~ year` with `x = ~ age` for Age profiles.
#'   Other frame/axis combinations produce a single ordinary plot.
#' @details
#' The two views use separate widgets, without refitting or changing the input
#' table. Projection boundaries appear only on year axes; projection years are
#' identified in the age-profile selector and default hover text. Multiple
#' displays can be embedded in the same HTML document independently.
#' A single-age table instead returns an annual Plotly trend labeled with its
#' age, without animation or view buttons.
#' Use [plot_trend()] or [plot_obs_pred()] directly for a single view; their
#' returned Plotly objects remain available for adding layers or changing layout.
#' @return A browsable HTML container holding two Plotly widgets for age-year
#'   views. Other frame/axis combinations return the supplied plotting function's
#'   Plotly object, as do single-age tables. Modify individual plots through
#'   `plot` before assembly.
#' @examples
#' d <- expand.grid(year = 2000:2004, age = 1:3)
#' d$est <- exp(0.1 * (d$year - 2000) + 0.2 * d$age)
#' annual <- plot_age_year(d, ylab = "Abundance")
#' profiles <- plot_age_year(d, x = ~age, xlab = "Age", frame = ~year)
#' d$pred <- d$est
#' d$obs <- d$est * 1.1
#' catch_views <- plot_age_year(d, plot = plot_obs_pred, ylab = "Catch")
#' # Print a view to display it, or save it with htmltools::save_html().
#' @seealso [plot_trend()], [plot_obs_pred()], [add_plot_buttons()], [vis_tam()]
#' @export
plot_age_year <- function(data, plot = plot_trend, x = ~year, xlab = "Year",
                          frame = ~age, ...) {
  args <- list(...)
  eligible <- all(c("year", "age") %in% names(data)) &&
    ((identical(all.vars(x), "year") && identical(all.vars(frame), "age")) ||
      (identical(all.vars(x), "age") && identical(all.vars(frame), "year")))
  if (!eligible || !nrow(data)) {
    return(do.call(plot, c(list(data = data, x = x, xlab = xlab, frame = frame), args)))
  }
  data <- data[order(data$year, data$age), , drop = FALSE]
  if (length(unique(data$age)) == 1L) {
    if (is.null(args$title)) args$title <- paste("Age", unique(data$age))
    if (is.null(args$text) && is.null(args$hovertemplate)) {
      data$.view_hover <- paste0("Year: ", data$year, "<br>Age: ", data$age, .forecast_hover(data))
      args$text <- ~.view_hover
    }
    return(do.call(plot, c(list(data = data, x = ~year, xlab = "Year", frame = NULL), args)))
  }
  projection <- if ("is_proj" %in% names(data)) data$is_proj else rep(FALSE, nrow(data))
  years <- unique(data$year)
  suffix <- vapply(years, function(year) {
    flags <- projection[data$year == year]
    if (all(flags)) " (projection)" else if (any(flags)) " (includes projections)" else ""
  }, character(1))
  year_labels <- paste0(data$year, suffix[match(data$year, years)])
  data$.view_year <- factor(year_labels, levels = unique(year_labels))
  if (is.null(args$text) && is.null(args$hovertemplate)) {
    data$.view_hover <- paste0(
      "Year: ", data$year,
      ifelse(projection, " (projection)", ""), "<br>Age: ", data$age,
      .forecast_hover(data)
    )
    args$text <- ~.view_hover
  }
  make_view <- function(x, xlab, frame, prefix, view) {
    animated <- length(unique(data[[all.vars(frame)]])) > 1L
    if (!animated) frame <- NULL
    p <- do.call(plot, c(list(data = data, x = x, xlab = xlab, frame = frame), args))
    p <- .age_year_menu(p, view)
    if (animated && length(p$x$frames)) {
      p <- plotly::animation_slider(p,
        currentvalue = list(prefix = prefix), pad = list(t = 45)
      )
    }
    p$width <- "100%"
    p$height <- "100%"
    p
  }
  annual <- make_view(~year, "Year", ~age, "Age: ", "annual")
  profile <- make_view(~age, "Age", ~.view_year, "Year: ", "profile")
  initial <- if (identical(all.vars(frame), "age")) "annual" else "profile"
  panel <- function(view, widget) {
    htmltools::tags$div(
      class = "tam-age-year-panel", `data-view` = view,
      style = paste0("flex:1;min-height:0;", if (view != initial) "display:none;"), widget
    )
  }
  script <- paste(readLines(system.file("js", "age-year.js", package = "tinyAM"),
    warn = FALSE
  ), collapse = "\n")
  htmltools::browsable(htmltools::tags$div(
    class = "tam-age-year", `data-default-view` = initial,
    style = "display:flex;flex-direction:column;width:100%;height:100%;min-height:400px;",
    panel("annual", annual), panel("profile", profile),
    htmltools::tags$script(htmltools::HTML(script))
  ))
}

.age_year_menu <- function(p, view) {
  p <- plotly::plotly_build(p)
  menus <- p$x$layout$updatemenus
  controls <- Filter(function(menu) {
    any(unlist(lapply(menu$buttons, `[[`, "label")) %in%
      c("Linear", "Log", "Show CI", "Hide CI"))
  }, menus)
  menu <- if (length(controls)) controls[[length(controls)]] else .log_linear_buttons[[1L]]
  menu$name <- "tam-age-year-view"
  menu$pad$t <- if (is.null(menu$pad$t)) 75 else menu$pad$t + 75
  menu$active <- match(view, c("annual", "profile")) - 1L
  menu$buttons <- list(
    list(method = "skip", execute = FALSE, args = list("annual"), label = "Annual trends"),
    list(method = "skip", execute = FALSE, args = list("profile"), label = "Age profiles")
  )
  if (is.null(p$x$layout)) p$x$layout <- list()
  p$x$layout$updatemenus <- c(menus, list(menu))
  p
}

