#' Explore fitted process deviations and joint-sample residuals
#'
#' These plots consume a single [tidy_tam()] result. They never
#' refit a model or draw another sample. Recruitment uses signed line plots;
#' N/F/M use selectable matrices or diagnostic bubbles.
#'
#' @param out A result from [tidy_tam()].
#' @param component Either `"NFM"` for selectable N/F/M matrices or bubbles,
#'   or `"Recruitment"` for a recruitment line plot.
#' @details Fitted deviations retain their original log scale: departures from
#' cohort/mean predictions for IID/AR1 and annual increments for RW/cor_rw.
#' Recruitment retains the residual definition in [tidy_recruitment()]. Available
#' recruitment intervals describe fitted deviations, not diagnostic draws.
#'
#' Joint-sample residuals are standardized innovations, evaluated at the fitted
#' mode and perturbed using the AD Jacobian and one draw from the joint Gaussian
#' parameter uncertainty. This is a local Gaussian diagnostic approximation,
#' not a posterior simulation or an independent draw for each process. Zero
#' and +/-2 lines aid interpretation; no confidence ribbons are attached to
#' diagnostic draws. Missing/invalid joint uncertainty is reported explicitly.
#'
#' Each call returns one plot; the dashboard supplies its dropdown containers.
#' N/F/M buttons retain the selected model. Signed matrices share symmetric colour ranges across
#' models within a component and have no logarithmic scale. M uses unique,
#' ordered state blocks. Bubble areas use absolute residuals with one scale
#' across all components/models; orange is negative and navy is positive.
#'
#' @return A Plotly htmlwidget, usable alone or in a flexdashboard chart.
#' @examples
#' \dontrun{
#' out <- tidy_tam(scaa, sam, ncam, process_seed = 123456L)
#' plot_process_deviations(out)
#' plot_process_residuals(out)
#' plot_process_deviations(out, component = "Recruitment")
#' plot_process_residuals(out, component = "Recruitment")
#' out$process_residuals$status
#' }
#' @seealso [plot_matrix()], [vis_tam()]
#' @export
plot_process_deviations <- function(out, component = c("NFM", "Recruitment")) {
  .plot_process(out, match.arg(component), residuals = FALSE)
}

#' @rdname plot_process_deviations
#' @export
plot_process_residuals <- function(out, component = c("NFM", "Recruitment")) {
  .plot_process(out, match.arg(component), residuals = TRUE)
}

.process_unavailable <- function() {
  plotly::layout(plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(), y = numeric()),
    xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
    annotations = .process_empty_annotation(), margin = list(t = 20))
}

.process_empty_annotation <- function() list(list(
  text = "Unavailable for this model component", x = .5, y = .5,
  xref = "paper", yref = "paper", showarrow = FALSE))

.process_models <- function(out) {
  label <- attr(out, "label", exact = TRUE)
  tables <- c(out$process_residuals["status"], out$deviations, out$process_residuals, out$recruitment, out$random_par)
  if (is.null(label)) return("Model")
  values <- unlist(lapply(tables, function(d) {
    if (!is.data.frame(d) || !label %in% names(d)) return(NULL)
    x <- d[[label]]
    if (is.factor(x)) levels(x) else as.character(x)
  }), use.names = FALSE)
  unique(values)
}

.process_plot_data <- function(d, out) {
  if (!is.data.frame(d) || !nrow(d)) return(NULL)
  label <- attr(out, "label", exact = TRUE)
  d$model <- if (is.null(label)) "Model" else as.character(d[[label]])
  d$model <- factor(d$model, levels = .process_models(out))
  if ("is_proj" %in% names(d)) d <- d[!d$is_proj, , drop = FALSE]
  d
}

.process_lines <- function(d, models, residuals) {
  if (is.null(d) || !nrow(d)) return(.process_unavailable())
  p <- plotly::plot_ly()
  cols <- tam_pal(length(models))
  for (i in seq_along(models)) {
    z <- d[as.character(d$model) == models[i], , drop = FALSE]
    if (!nrow(z)) next
    z <- z[order(z$year), , drop = FALSE]
    if (!residuals && all(c("lwr", "upr") %in% names(z)) && any(is.finite(z$lwr) & is.finite(z$upr))) {
      p <- plotly::add_trace(p, x = c(z$year, rev(z$year)), y = c(z$upr, rev(z$lwr)),
        type = "scatter", mode = "lines", fill = "toself", fillcolor = .process_alpha(cols[i]),
        line = list(color = "transparent"), showlegend = FALSE, hoverinfo = "skip", legendgroup = models[i])
    }
    p <- plotly::add_trace(p, x = z$year, y = z$est, type = "scatter", mode = "lines+markers",
      name = models[i], legendgroup = models[i], line = list(color = cols[i]),
      marker = list(color = cols[i], size = 4), hovertemplate = "Year: %{x}<br>Value: %{y}<extra>%{fullData.name}</extra>")
  }
  refs <- if (residuals) c(-2, 0, 2) else 0
  plotly::layout(p, xaxis = list(title = "Year"),
    yaxis = list(title = if (residuals) "Joint-sample standardized residual" else "Log recruitment deviation", zeroline = FALSE),
    shapes = lapply(refs, function(y) list(type = "line", xref = "paper", x0 = 0, x1 = 1,
      y0 = y, y1 = y, line = list(color = "#888888", width = 1, dash = if (y == 0) "solid" else "dot"))),
    margin = list(l = 65, r = 20, t = 20, b = 50))
}

.process_alpha <- function(col) {
  rgb <- grDevices::col2rgb(col)
  sprintf("rgba(%d,%d,%d,0.15)", rgb[1L], rgb[2L], rgb[3L])
}

.process_nfm <- function(tables, out, models) {
  pieces <- lapply(c("N", "F", "M"), function(component) {
    d <- .process_plot_data(tables[[component]], out)
    if (is.null(d)) d <- data.frame(year = numeric(), est = numeric(), model = character())
    d$coordinate <- if (component == "M") integer(nrow(d)) else if (nrow(d)) d$age else numeric()
    if (component == "M" && nrow(d)) {
      d$coordinate <- integer(nrow(d))
      for (model in models) {
        i <- which(d$model == model)
        d$coordinate[i] <- match(as.character(d$age_block[i]), unique(as.character(d$age_block[i])))
      }
    }
    d$coordinate_label <- if (component == "M") as.character(d$age_block) else as.character(d$age)
    d$component <- rep(component, nrow(d))
    d <- d[intersect(c("year", "est", "model", "coordinate", "coordinate_label", "component", "definition", "notes"), names(d))]
    # NA-only rows keep every model/component selectable, without inventing values.
    for (model in setdiff(models, as.character(d$model))) {
      d <- .bind_rows(list(d, data.frame(year = 0, est = NA_real_, model = model,
        coordinate = 1, coordinate_label = "", component = component)))
    }
    d
  })
  d <- .bind_rows(pieces)
  years <- d$year[is.finite(d$est)]
  d$year[d$year == 0] <- if (length(years)) min(years) else 0
  d$model <- factor(d$model, levels = models)
  d$component <- factor(d$component, levels = c("N", "F", "M"))
  d
}

.process_bubbles <- function(d, models) {
  limit <- max(c(.01, abs(d$est[is.finite(d$est)])))
  colours <- c(tail(tam_cols(), 1L), tam_cols()[1L])
  span <- function(x) if (length(x)) range(x) + c(-.5, .5) else c(0, 1)
  components <- lapply(c("N", "F", "M"), function(component) {
    values <- d[as.character(d$component) == component, , drop = FALSE]
    xrange <- span(values$year[is.finite(values$est)])
    cells <- lapply(models, function(model) {
      z <- values[values$model == model & is.finite(values$est), , drop = FALSE]
      traces <- list()
      if (nrow(z)) {
        z$kind <- paste0("Model: ", model, "<br>", component, " residual<br>",
          if (component == "M") "M state block: " else "Age: ", z$coordinate_label,
          "<br>Year: ", z$year, "<br>Residual")
        p <- plotly::plotly_build(plot_bubbles(z, res = ~est, y = ~coordinate,
          title = "", colors = colours))
        traces <- p$x$data
      }
      # Two stable sign traces in every frame prevent stale markers on empty frames.
      traces <- lapply(c("-", "+"), function(sign) {
        found <- Filter(function(trace) identical(trace$name, sign), traces)
        trace <- if (length(found)) found[[1L]] else list(type = "scatter", mode = "markers",
          x = numeric(), y = numeric(), marker = list(), name = sign)
        trace$showlegend <- FALSE
        trace$marker$size <- abs(z$est[if (sign == "+") z$est > 0 else z$est <= 0])
        trace$marker$color <- colours[if (sign == "-") 1L else 2L]
        trace$marker$sizeref <- 2 * limit / 30^2
        trace$marker$sizemode <- "area"
        trace$marker$sizemin <- 0
        trace$frame <- NULL
        trace
      })
      ticks <- sort(unique(z$coordinate))
      axis <- list(title = if (component == "M") "M state block" else "Age",
        type = "linear", range = span(z$coordinate), autorange = FALSE,
        tickmode = if (component == "M") "array" else "auto",
        tickvals = ticks, ticktext = z$coordinate_label[match(ticks, z$coordinate)])
      list(data = traces, layout = list(xaxis = list(title = "Year", range = xrange, autorange = FALSE),
        yaxis = axis, annotations = if (nrow(z)) list() else .process_empty_annotation()), available = nrow(z) > 0L)
    })
    list(name = component, cells = cells)
  })
  initial_component <- which(vapply(components, function(x) any(vapply(x$cells, `[[`, logical(1), "available")), logical(1)))[1L]
  if (is.na(initial_component)) initial_component <- 1L
  initial_model <- which(vapply(components[[initial_component]]$cells, `[[`, logical(1), "available"))[1L]
  if (is.na(initial_model)) initial_model <- 1L
  config <- list(models = models, components = components, component = initial_component - 1L, model = initial_model - 1L)
  cell <- components[[initial_component]]$cells[[initial_model]]
  p <- plotly::plotly_build(plotly::plot_ly(type = "scatter", x = numeric(), y = numeric(), mode = "markers"))
  p$x$data <- cell$data
  p$x$layout <- utils::modifyList(p$x$layout, cell$layout)
  script <- paste(readLines(system.file("js", "process-bubbles.js", package = "tinyAM"), warn = FALSE), collapse = "\n")
  htmlwidgets::onRender(p, script, data = config)
}

.plot_process <- function(out, component, residuals) {
  if (!is.list(out) || !all(c("recruitment", "deviations") %in% names(out))) {
    cli::cli_abort("{.arg out} must be a result from tidy_tam().")
  }
  models <- .process_models(out)
  if (!length(models)) models <- "Model"
  tables <- if (residuals) out$process_residuals else out$deviations
  if (component == "Recruitment") {
    rec <- .process_plot_data(if (residuals) tables$recruitment else out$recruitment$residual, out)
    if (!is.null(rec)) rec <- rec[is.finite(rec$est), , drop = FALSE]
    return(.process_lines(rec, models, residuals))
  }
  d <- .process_nfm(tables, out, models)
  if (residuals) return(.process_bubbles(d, models))
  plot_matrix(d, metric = ~component, frame = ~model, signed = TRUE,
    y = ~coordinate, y_labels = ~coordinate_label,
    ylab = c(N = "Age", F = "Age", M = "M state block"))
}
