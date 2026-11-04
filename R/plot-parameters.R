#' @rdname plot_tam
#' @export
plot_par <- function(data, ...) {
  data$coef <- factor(ifelse(is.na(data$coef), data$par, paste0(data$par, ": ", data$coef)))

  args <- list(...)
  if (is.null(args$color)) {
    args$color <- I("#377EB8")
  }
  p <- do.call(plot_ly, c(list(data = data, x = ~est, y = ~coef), args))

  if (any(is.finite(data$lwr) & is.finite(data$upr))) {
    p <- p |>
      add_segments(
        x = ~lwr, xend = ~upr,
        y = ~coef, yend = ~coef,
        showlegend = FALSE
      )
  }
  p |>
    add_markers(showlegend = TRUE) |>
    layout(
      xaxis = list(
        title = "Estimate"
      ),
      yaxis = list(
        title = ""
      )
    )
}

# Signed formula effects use a linear axis and a zero reference. Keep the
# existing trend/ribbon implementation for temporal effects.
.plot_formula_effect <- function(data, points = FALSE, frame = NULL, ...) {
  if (!nrow(data)) {
    return(plotly::plot_ly())
  }
  limits <- range(c(0, data$est, data$lwr, data$upr), finite = TRUE)
  pad <- max(diff(limits) * .05, .01)
  if (points) {
    p <- plotly::plot_ly(data,
      x = ~coordinate, y = ~est, frame = frame,
      error_y = list(
        type = "data", array = ~ upr - est,
        arrayminus = ~ est - lwr, symmetric = FALSE
      ), ...
    ) |>
      plotly::add_markers()
  } else {
    p <- plot_trend(data,
      x = ~coordinate, xlab = unique(data$variable)[1L],
      frame = frame, add_buttons = FALSE, ...
    )
  }
  shapes <- list(list(
    type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0,
    line = list(color = "grey", width = 1)
  ))
  if (any(data$is_proj)) {
    boundary <- max(data$coordinate[!data$is_proj])
    shapes <- c(shapes, list(list(
      type = "line", x0 = boundary, x1 = boundary,
      y0 = 0, y1 = 1, yref = "paper", line = list(color = "black", dash = "dot", width = 1)
    )))
  }
  plotly::layout(p,
    xaxis = list(
      title = unique(data$variable)[1L],
      tickvals = unique(data$coordinate), ticktext = unique(data$level)
    ),
    yaxis = list(
      title = if ("component" %in% names(data) && unique(data$component) != "q") {
        if (startsWith(unique(data$component), "sd_")) {
          "Effect on log SD"
        } else {
          paste("Effect on log", unique(data$component))
        }
      } else {
        "Effect on catchability's link scale"
      }, type = "linear",
      range = limits + c(-pad, pad)
    ), shapes = shapes
  )
}

