# Observation plots consume retained residuals; selecting a view never invokes
# residuals() or oneStepPredict(). Missing OSA values stay missing within a model.
.observation_residual_data <- function(data, type, selection = NULL) {
  values <- function(name) if (name %in% names(data)) data[[name]] else rep(NA_real_, nrow(data))
  standardized <- values("std_res")
  osa <- values("osa_res")
  eligible <- rep(TRUE, nrow(data))
  if ("obs" %in% names(data)) eligible <- eligible & is.finite(data$obs) & data$obs > 0
  if ("is_proj" %in% names(data)) eligible <- eligible & !is.na(data$is_proj) & !data$is_proj
  osa[!eligible] <- standardized[!eligible] <- NA_real_
  model <- if ("model" %in% names(data)) as.character(data$model) else rep("Model", nrow(data))
  models <- if (is.factor(data$model)) levels(data$model) else unique(model)
  models <- models[models %in% model & !is.na(models)]
  if (!length(models)) models <- "Model"
  method <- if ("osa_method" %in% names(data)) as.character(data$osa_method) else rep(NA_character_, nrow(data))
  method[is.na(method) | !nzchar(method)] <- "method unspecified"
  data$.residual_value <- rep(NA_real_, nrow(data))
  data$.residual_type <- rep(NA_character_, nrow(data))
  data$.residual_source <- rep(NA_character_, nrow(data))
  status <- lapply(models, function(name) {
    rows <- which(model == name)
    selected <- if (type == "auto" && !is.null(selection)) selection[[name]] else if (type == "auto") {
      if (any(is.finite(osa[rows]))) "osa" else "standardized"
    } else type
    data$.residual_value[rows] <<- if (selected == "osa") osa[rows] else standardized[rows]
    data$.residual_type[rows] <<- selected
    label <- if (selected == "osa") paste0("OSA residuals (", method[rows], ")") else "Standardized residuals"
    data$.residual_source[rows] <<- label
    finite <- rows[is.finite(data$.residual_value[rows])]
    data.frame(model = name, residual_type = selected,
      osa_method = if (selected == "osa" && length(finite)) paste(unique(method[finite]), collapse = ", ") else NA_character_,
      n = length(finite), stringsAsFactors = FALSE)
  })
  list(data = data, status = do.call(rbind, status))
}

.residual_axis_range <- function(data, formula) {
  value <- stats::model.frame(formula, data, na.action = stats::na.pass)[[1L]]
  if (!is.numeric(value)) return(NULL)
  value <- value[is.finite(value)]
  if (!length(value)) return(NULL)
  limits <- range(value)
  limits + c(-1, 1) * max(diff(limits) * .05, .1)
}

.residual_frame_values <- function(data, frame) {
  if (is.null(frame)) return(list(values = rep("Model", nrow(data)), levels = "Model", prefix = ""))
  value <- if (inherits(frame, "formula")) {
    tab <- stats::model.frame(frame, data, na.action = stats::na.pass)
    if (ncol(tab) != 1L) cli::cli_abort("Residual plot frames must select one column.")
    tab[[1L]]
  } else if (is.character(frame) && length(frame) == 1L && frame %in% names(data)) {
    data[[frame]]
  } else frame
  if (length(value) != nrow(data) || anyNA(value)) {
    cli::cli_abort("Residual plot frames must provide one nonmissing value per row.")
  }
  levels <- if (is.factor(value)) levels(value)[levels(value) %in% value] else unique(as.character(value))
  if (!length(levels)) levels <- "Model"
  label <- if (inherits(frame, "formula")) all.vars(frame)[1L] else "Frame"
  if (!length(label) || is.na(label)) label <- "Frame"
  list(values = as.character(value), levels = levels,
    prefix = paste0(toupper(substr(label, 1L, 1L)), substring(label, 2L), ": "))
}

.residual_plot_label <- function(status, type) {
  available <- status[status$n > 0L, , drop = FALSE]
  kinds <- unique(available$residual_type)
  if (identical(kinds, "osa")) "OSA residuals" else if (identical(kinds, "standardized")) {
    "Standardized residuals"
  } else if (type == "osa") "OSA residuals" else if (type == "standardized") {
    "Standardized residuals"
  } else "Residuals (Auto)"
}

.residual_plot_annotations <- function(status, type) {
  available <- status$n > 0L
  if (!any(available)) {
    text <- if (type == "osa" || all(status$residual_type == "osa")) "OSA residuals unavailable for this model/panel" else
      "Residuals unavailable for this model/panel"
    return(list(list(text = text, x = .5, y = .5, xref = "paper", yref = "paper", showarrow = FALSE)))
  }
  labels <- ifelse(status$residual_type == "osa", paste0("OSA (", status$osa_method, ")"), "Standardized")
  labels[!available] <- if (type == "osa") "OSA unavailable" else "Unavailable"
  if (nrow(status) > 1L) labels <- paste0(status$model, ": ", labels)
  text <- paste(labels, collapse = "; ")
  if (type == "auto") text <- paste0("Auto: ", text)
  list(list(text = as.character(htmltools::htmlEscape(text)), x = 0, y = 1,
    xref = "paper", yref = "paper", xanchor = "left", yanchor = "bottom",
    yshift = 8, showarrow = FALSE, font = list(size = 10)))
}

.residual_plot_cell <- function(data, type, kind, args, x, xlab, y, ylab, title, selection) {
  selected <- .observation_residual_data(data, type, selection)
  d <- selected$data
  label <- .residual_plot_label(selected$status, type)
  d <- d[is.finite(d$.residual_value), , drop = FALSE]
  if (!nrow(d)) {
    p <- plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(), y = numeric())
  } else if (kind == "bubbles") {
    d$kind <- as.character(htmltools::htmlEscape(d$.residual_source))
    p <- do.call(plot_bubbles, c(list(data = d, res = ~.residual_value,
      y = y, ylab = ylab, title = if (is.null(title)) label else title), args))
  } else {
    d$.residual_hover <- paste0(htmltools::htmlEscape(d$.residual_source), ": ", signif(d$.residual_value, 4))
    if (is.null(args$text)) args$text <- ~.residual_hover
    if (is.null(args$hovertemplate)) args$hovertemplate <- "%{x}<br>%{text}<extra>%{fullData.name}</extra>"
    p <- do.call(plotly::plot_ly, c(list(data = d, x = x, y = ~.residual_value), args)) |>
      plotly::add_markers(alpha = .8)
  }
  p <- plotly::plotly_build(p)
  if (!nrow(d)) for (j in seq_along(p$x$data)) {
    p$x$data[[j]]$x <- p$x$data[[j]]$y <- numeric()
    p$x$data[[j]]$showlegend <- FALSE
  }
  p$x$layout$xaxis$title <- if (kind == "bubbles") "Year" else xlab
  p$x$layout$yaxis$title <- if (kind == "bubbles") ylab else label
  if (kind == "bubbles") p$x$layout$title <- if (is.null(title)) label else title
  p$x$layout$annotations <- .residual_plot_annotations(selected$status, type)
  list(data = p$x$data, layout = p$x$layout, status = selected$status)
}

.residual_trace_key <- function(trace) {
  paste(trace$name, trace$legendgroup, trace$type, trace$mode, collapse = "\r")
}

.plot_observation_residuals <- function(data, type, kind, args,
                                      x = ~year, xlab = "Year", y = ~age,
                                      ylab = "Age", title = NULL) {
  if (!is.data.frame(data)) cli::cli_abort("Residual plots require a data frame.")
  frame <- .residual_frame_values(data, args$frame)
  animated <- !is.null(args$frame) && length(frame$levels) > 1L
  args$frame <- NULL
  types <- c("auto", "standardized", "osa")
  models <- if (is.factor(data$model)) levels(data$model) else unique(as.character(data$model))
  model_colors <- NULL
  if (kind == "scatter" && inherits(args$color, "formula") && identical(all.vars(args$color), "model")) {
    palette <- args$colors
    if (is.null(palette)) palette <- tam_pal(length(models))
    if (is.function(palette)) palette <- palette(length(models))
    if (is.character(palette) && length(palette) >= length(models)) model_colors <- stats::setNames(palette[seq_along(models)], models)
  }
  sign_colors <- if (is.null(args$colors)) c(tam_cols()[1L], tail(tam_cols(), 1L)) else args$colors
  if (is.function(sign_colors)) sign_colors <- sign_colors(2L)
  views <- lapply(types, function(view) {
    chosen <- .observation_residual_data(data, view)
    selection <- stats::setNames(chosen$status$residual_type, chosen$status$model)
    chosen <- chosen$data
    valid <- chosen[is.finite(chosen$.residual_value), , drop = FALSE]
    xrange <- .residual_axis_range(valid, if (kind == "bubbles") ~year else x)
    yrange <- .residual_axis_range(valid, if (kind == "bubbles") y else ~.residual_value)
    cells <- lapply(frame$levels, function(level) {
      d <- if (animated) data[frame$values == level, , drop = FALSE] else data
      cell <- .residual_plot_cell(d, view, kind, args, x, xlab, y, ylab, title, selection)
      if (!is.null(xrange)) cell$layout$xaxis <- utils::modifyList(cell$layout$xaxis, list(range = xrange, autorange = FALSE))
      if (!is.null(yrange)) cell$layout$yaxis <- utils::modifyList(cell$layout$yaxis, list(range = yrange, autorange = FALSE))
      for (j in seq_along(cell$data)) {
        trace <- cell$data[[j]]
        if (kind == "bubbles" && length(trace$name) == 1L &&
            trace$name %in% c("+", "-") && is.character(sign_colors)) {
          trace$marker$color <- sign_colors[match(trace$name, c("+", "-"))]
        } else if (length(model_colors) && length(trace$name) == 1L && trace$name %in% names(model_colors)) {
          trace$marker$color <- unname(model_colors[trace$name])
        }
        cell$data[[j]] <- trace
      }
      cell
    })
    list(type = view, cells = cells)
  })
  # Stable trace slots clear markers when a model lacks the selected residuals.
  templates <- list()
  for (view in views) for (cell in view$cells) for (trace in cell$data) {
    key <- .residual_trace_key(trace)
    if (is.null(templates[[key]])) templates[[key]] <- trace
  }
  for (i in seq_along(views)) for (j in seq_along(views[[i]]$cells)) {
    cell <- views[[i]]$cells[[j]]
    keys <- vapply(cell$data, .residual_trace_key, character(1))
    cell$data <- lapply(names(templates), function(key) {
      found <- match(key, keys)
      if (!is.na(found)) return(cell$data[[found]])
      trace <- templates[[key]]
      trace$x <- trace$y <- numeric()
      trace$text <- character()
      trace$showlegend <- FALSE
      if (!is.null(trace$marker$size)) trace$marker$size <- numeric()
      trace
    })
    views[[i]]$cells[[j]] <- cell
  }
  config <- list(views = views, frames = frame$levels, animated = animated,
    prefix = frame$prefix, view = match(type, types) - 1L, frame = 0L)
  initial <- views[[match(type, types)]]$cells[[1L]]
  p <- plotly::plotly_build(plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(), y = numeric()))
  p$x$data <- initial$data
  p$x$layout <- utils::modifyList(p$x$layout, initial$layout)
  if (animated) {
    cells <- views[[match(type, types)]]$cells
    p$x$frames <- lapply(seq_along(frame$levels), function(i) list(
      name = frame$levels[i], data = cells[[i]]$data, layout = cells[[i]]$layout,
      traces = seq_along(cells[[i]]$data) - 1L))
  }
  script <- paste(readLines(system.file("js", "residual-selection.js", package = "tinyAM"), warn = FALSE), collapse = "\n")
  htmlwidgets::onRender(p, script, data = config)
}
