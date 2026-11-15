#' @name plot_tam
#' @rdname plot_tam
#' @export
plot_bubbles <- function(data, res = NULL, y = ~age, ylab = "Age", title = NULL,
                         residual_type = c("auto", "standardized", "osa", "quantile"), ...) {
  if (is.null(res)) {
    return(.plot_observation_residuals(data, match.arg(residual_type), "bubbles",
      list(...), y = y, ylab = ylab, title = title))
  }
  if (!is.null(res)) {
    if (!inherits(res, "formula") || length(res) != 2L) cli::cli_abort("res must be a one-sided formula selecting a numeric column.")
    values <- stats::model.frame(res, data, na.action = stats::na.pass)
    if (ncol(values) != 1L) cli::cli_abort("res must select one signed numeric column.")
    data$res <- values[[1L]]
    if (!is.numeric(data$res)) cli::cli_abort("res must select a signed numeric column.")
    if (is.null(title)) title <- "Process deviations"
  }
  data$sign <- ifelse(data$res > 0, "+", "-")
  # Break size ties deterministically without advancing the caller's RNG.
  data$abs_res <- abs(data$res) + seq_len(nrow(data)) * .Machine$double.eps
  data <- data[!is.na(data$res), ]
  # Unused factor levels create empty Plotly animation frames and can leave
  # the initial bubble chart blank when some models lack this process.
  data <- droplevels(data)
  data$.bubble_text <- if ("kind" %in% names(data)) paste0(data$kind, ": ", round(data$res, 2)) else as.character(round(data$res, 2))

  args <- list(...)
  if (is.null(args$colors)) {
    args$colors <- c(utils::head(tam_cols(), 1), tail(tam_cols(), 1))
  }
  if (is.null(args$text)) args$text <- ~.bubble_text

  p <- do.call(
    plotly::plot_ly,
    c(
      list(
        data = data,
        x = ~year,
        y = y,
        marker = list(
          line = list(width = 0),
          size     = ~abs_res,
          sizemin  = 1,
          sizeref  = 0.01,
          sizemode = "area"
        ),
        color = ~sign
      ),
      args
    )
  )

  p |>
    plotly::add_markers() |>
    plotly::layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = ylab)
    )
}


#' @name plot_tam
#' @rdname plot_tam
#' @export
plot_resid <- function(
  data,
  x = ~year,
  xlab = "Year",
  residual_type = c("auto", "standardized", "osa", "quantile"),
  ...
) {
  .plot_observation_residuals(data, match.arg(residual_type), "scatter",
    list(...), x = x, xlab = xlab)
}

# Observation plots consume retained residuals; selecting a view never invokes
# residuals() or oneStepPredict(). Missing OSA values stay missing within a model.
.observation_residual_data <- function(data, type, selection = NULL) {
  values <- function(name) if (name %in% names(data)) data[[name]] else rep(NA_real_, nrow(data))
  standardized <- values("std_res")
  osa <- values("osa_res")
  quantile <- values("quantile_res")
  pld <- seq_len(nrow(data)) %in% which(.is_pld_family(data$family))
  quantile_family <- pld | seq_len(nrow(data)) %in% which(data$family %in% c("gamma", "gengamma"))
  quantile[!quantile_family] <- standardized[!quantile_family]
  eligible <- rep(TRUE, nrow(data))
  if ("obs" %in% names(data)) eligible <- eligible & is.finite(data$obs) & (data$obs > 0 | pld & data$obs == 0)
  if ("is_proj" %in% names(data)) eligible <- eligible & !is.na(data$is_proj) & !data$is_proj
  osa[!eligible] <- standardized[!eligible] <- NA_real_
  osa[quantile_family] <- NA_real_
  if ("obs" %in% names(data)) standardized[which(data$obs == 0)] <- NA_real_
  quantile[!eligible] <- NA_real_
  quantile[which(data$is_structural %in% TRUE)] <- NA_real_
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
    has_quantile_family <- any(quantile_family[rows])
    selected <- if (type == "auto" && !is.null(selection)) selection[[name]] else if (type == "auto") {
      if (has_quantile_family) "quantile" else if (any(is.finite(osa[rows]))) "osa" else "standardized"
    } else type
    data$.residual_value[rows] <<- switch(selected, osa = osa[rows],
      quantile = quantile[rows], standardized = standardized[rows])
    data$.residual_type[rows] <<- selected
    label <- switch(selected, osa = paste0("OSA residuals (", method[rows], ")"),
      quantile = "Quantile residuals",
      standardized = ifelse(pld[rows], "Standardized positive-component residuals", "Standardized residuals"))
    data$.residual_source[rows] <<- label
    finite <- rows[is.finite(data$.residual_value[rows])]
    data.frame(model = name, residual_type = selected,
      osa_method = if (selected == "osa" && length(finite)) paste(unique(method[finite]), collapse = ", ") else NA_character_,
      has_quantile_family = has_quantile_family,
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

.zero_frequency_data <- function(data) {
  pld <- .is_pld_family(data$family)
  if (!any(pld)) return(NULL)
  rows <- which(pld & data$is_observed &
    !data$is_structural & is.finite(data$p_zero))
  if (!length(rows)) return(NULL)
  d <- data[rows, , drop = FALSE]
  groups <- split(d, list(d$model, d$age), drop = TRUE)
  .bind_rows(lapply(groups, function(x) data.frame(
    model = x$model[1L], age = x$age[1L], obs = mean(x$obs == 0),
    pred = mean(x$p_zero), n = nrow(x)
  )))
}

.residual_frame_values <- function(data, frame) {
  if (is.null(frame)) return(list(values = rep("Model", nrow(data)), levels = "Model", label = "Model"))
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
  label <- if (inherits(frame, "formula")) all.vars(frame)[1L] else if (length(frame) == 1L &&
    is.character(frame) && frame %in% names(data)) frame else "Frame"
  if (!length(label) || is.na(label)) label <- "Frame"
  list(values = as.character(value), levels = levels,
    label = paste0(toupper(substr(label, 1L, 1L)), substring(label, 2L)))
}

.residual_plot_label <- function(status, type) {
  available <- status[status$n > 0L, , drop = FALSE]
  kinds <- unique(available$residual_type)
  if (identical(kinds, "quantile") || type == "quantile") {
    "Quantile residuals"
  } else if (identical(kinds, "osa")) "OSA residuals" else if (identical(kinds, "standardized")) {
    "Standardized residuals"
  } else if (type == "osa") "OSA residuals" else if (type == "standardized") {
    "Standardized residuals"
  } else "Residuals"
}

.residual_plot_annotations <- function(status, type) {
  available <- status$n > 0L
  if (!any(available)) {
    text <- if (type == "osa" || all(status$residual_type == "osa")) "OSA residuals unavailable for this model/panel" else
      "Residuals unavailable for this model/panel"
    return(list(list(text = text, x = .5, y = .5, xref = "paper", yref = "paper", showarrow = FALSE)))
  }
  if (all(available)) return(list())
  message <- if (type == "osa") "OSA residuals unavailable" else "Residuals unavailable"
  text <- paste0(htmltools::htmlEscape(status$model[!available]), ": ", message, collapse = "<br>")
  list(list(text = text, x = 0, y = 1,
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
      y = y, ylab = ylab, title = if (is.null(title)) "" else title), args))
  } else {
    d$.residual_hover <- paste0(htmltools::htmlEscape(d$.residual_source), ": ", signif(d$.residual_value, 4))
    if (is.null(args$text)) args$text <- ~.residual_hover
    if (is.null(args$hovertemplate)) args$hovertemplate <- "%{x}<br>%{text}<extra>%{fullData.name}</extra>"
    p <- do.call(plotly::plot_ly, c(list(data = d, x = x, y = ~.residual_value), args)) |>
      plotly::add_markers(alpha = .8)
  }
  p <- plotly::plotly_build(p)
  if (!nrow(d)) p$x$data <- list()
  p$x$layout$xaxis$title <- if (kind == "bubbles") "Year" else xlab
  p$x$layout$yaxis$title <- if (kind == "bubbles") ylab else label
  if (kind == "bubbles") p$x$layout$title <- if (is.null(title)) "" else title
  p$x$layout$annotations <- .residual_plot_annotations(selected$status, type)
  list(data = p$x$data, layout = p$x$layout, status = selected$status)
}

.plot_observation_residuals <- function(data, type, kind, args,
                                      x = ~year, xlab = "Year", y = ~age,
                                      ylab = "Age", title = NULL) {
  if (!is.data.frame(data)) cli::cli_abort("Residual plots require a data frame.")
  # Keep explicit Plotly row attributes aligned when frames or OSA gaps subset
  # the data. Formulas and scalar style arguments already have the right scope.
  for (field in c("text", "hovertext", "customdata", "ids", "key", "hovertemplate",
                  "color", "symbol", "size", "opacity")) {
    value <- args[[field]]
    rowwise <- !inherits(value, "formula") &&
      ((is.matrix(value) && nrow(value) == nrow(data)) ||
       ((is.atomic(value) || is.list(value)) && is.null(dim(value)) &&
        length(value) > 1L && length(value) == nrow(data)))
    if (rowwise) {
      column <- paste0(".tam_residual_", field)
      data[[column]] <- value
      args[[field]] <- stats::reformulate(column)
    }
  }
  frame <- .residual_frame_values(data, args$frame)
  has_frames <- !is.null(args$frame) && length(frame$levels) > 1L
  model_frame <- identical(args$frame, "model") ||
    (inherits(args$frame, "formula") && identical(all.vars(args$frame), "model"))
  args$frame <- NULL
  available <- .observation_residual_data(data, "auto")$status
  types <- c("standardized", "osa")
  if (any(available$has_quantile_family) || type == "quantile") types <- c(types, "quantile")
  # Resolve once for the complete panel, never from a partially observed frame.
  default_type <- if (any(available$has_quantile_family)) "quantile" else
    if (all(available$residual_type == "osa")) "osa" else "standardized"
  model_rows <- match(frame$levels, available$model)
  defaults <- if (model_frame) available$residual_type[model_rows] else
    rep(default_type, length(frame$levels))
  defaults[is.na(defaults)] <- "standardized"
  combine_quantile <- if (model_frame) !available$has_quantile_family[model_rows] else
    rep(!any(available$has_quantile_family), length(frame$levels))
  combine_quantile[is.na(combine_quantile)] <- TRUE
  default_views <- match(defaults, types) - 1L
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
      d <- if (has_frames) data[frame$values == level, , drop = FALSE] else data
      cell <- .residual_plot_cell(d, view, kind, args, x, xlab, y, ylab, title, selection)
      if (!is.null(xrange)) cell$layout$xaxis <- utils::modifyList(cell$layout$xaxis, list(range = xrange, autorange = FALSE))
      if (!is.null(yrange)) cell$layout$yaxis <- utils::modifyList(cell$layout$yaxis, list(range = yrange, autorange = FALSE))
      for (j in seq_along(cell$data)) {
        trace <- cell$data[[j]]
        if (kind == "bubbles") trace$marker$line <- list(width = 0)
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
  initial_view <- if (type == "auto") default_views[1L] else match(type, types) - 1L
  selectors <- which(types != "osa" | any(available$residual_type == "osa"))
  # Show OSA on explicit selection even when its retained values are unavailable.
  if (type == "osa") selectors <- seq_along(types)
  config <- list(views = views, frames = frame$levels, frameLabel = frame$label, view = initial_view,
    automatic = type == "auto", defaultViews = default_views,
    compactHeader = kind == "bubbles" && is.null(title),
    selectors = selectors - 1L, labels = c("Standardized", "OSA", "Quantile")[seq_along(types)],
    combineQuantile = combine_quantile,
    showSelector = (has_frames && model_frame) || any(available$residual_type == "osa") ||
      any(available$has_quantile_family))
  initial <- views[[initial_view + 1L]]$cells[[1L]]
  p <- plotly::plotly_build(plotly::plot_ly(type = "scatter", mode = "markers", x = numeric(), y = numeric()))
  p$x$data <- initial$data
  p$x$layout <- utils::modifyList(p$x$layout, initial$layout)
  script <- paste(readLines(system.file("js", "residual-selection.js", package = "tinyAM"), warn = FALSE), collapse = "\n")
  htmlwidgets::onRender(p, script, data = config)
}
