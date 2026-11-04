#' @rdname plot_tam
#' @export
plot_heatmap <- function(
  data,
  zlab = "Estimate",
  title = "",
  ...
) {
  unit_label <- .unit_label(data, "Estimate")
  if (missing(zlab)) zlab <- unit_label
  args <- list(...)
  if ("frame" %in% names(args) && !is.null(args$frame)) {
    cli::cli_abort(
      "Argument {.arg frame} is not supported for heatmaps.
       Supply a single frame at a time, or split/filter your data by the grouping
       variable and call {.fn plot_heatmap} per group."
    )
  }

  plot_ly(data,
    x = ~year, y = ~age, z = ~est,
    colorbar = list(title = zlab), ...
  ) |>
    add_heatmap() |>
    layout(
      title = title,
      xaxis = list(title = "Year"),
      yaxis = list(title = "Age")
    )
}


#' Explore age-year patterns
#'
#' Display a complete or sparse age-year table as a heatmap, surface or contour.
#' Buttons change the metric, display type and value scale; an optional model
#' selector compares models using the same ranges within each metric.
#'
#' @param data A data frame with numeric `year`, `age` and `est` columns.
#'   Optional `unit`, `definition`, `notes` and `is_proj` columns enrich hover text.
#' @param metric Optional one-sided formula selecting a metric column.
#' @param frame Optional one-sided formula selecting a model column. Factor
#'   levels determine the initial model and selector order.
#' @param type Initial display: `"heatmap"`, `"surface"` or `"contour"`.
#' @param signed Use symmetric orange-white-navy ranges centred on zero and
#'   disable the log control? Default `FALSE` preserves ordinary matrix plots.
#' @param y One-sided formula selecting the vertical coordinate, default `~age`.
#'   Character/factor coordinates use ordered positions with each model's own
#'   labels, suitable for unique M state blocks.
#' @param ylab Vertical-axis label, default `"Age"`. A named character vector
#'   can supply a different label for each metric.
#' @param y_labels Optional one-sided formula selecting display labels for the
#'   vertical coordinates, useful for model-specific state blocks.
#' @details Linear shows original estimates, including zero. Log displays
#'   log10 estimates; nonpositive values are gaps. Hover text always shows
#'   original values and units. Ranges are shared across models separately for
#'   each metric. All models of a metric must use the same units.
#'   Colours use the built-in [tam_cols()] ramp from low to high values.
#'   Colour-ramp titles place units on a second line; axis labels remain on one line.
#'   Missing cells remain missing. Surface and contour rendering interpolates
#'   only for display and does not create additional source estimates.
#' @return A Plotly object with independent controls for each widget.
#' @examples
#' d <- expand.grid(year = 2000:2005, age = 1:5)
#' d$est <- exp(.1 * (d$year - 2000) - .3 * d$age)
#' plot_matrix(d)
#' @seealso [plot_heatmap()], [vis_tam()]
#' @export
plot_matrix <- function(data, metric = NULL, frame = NULL,
                        type = c("heatmap", "surface", "contour"),
                        signed = FALSE, y = ~age, ylab = "Age", y_labels = NULL) {
  type <- match.arg(type)
  if (!is.logical(signed) || length(signed) != 1L || is.na(signed)) cli::cli_abort("signed must be TRUE or FALSE.")
  if (!inherits(y, "formula") || length(all.vars(y)) != 1L || !all.vars(y) %in% names(data)) {
    cli::cli_abort("y must select one existing vertical-axis column.")
  }
  vertical <- data[[all.vars(y)]]
  categorical <- !is.numeric(vertical)
  if (anyNA(vertical)) cli::cli_abort("Vertical coordinates cannot be missing.")
  data$.vertical_label <- as.character(vertical)
  if (!is.null(y_labels)) {
    if (!inherits(y_labels, "formula") || length(all.vars(y_labels)) != 1L ||
        !all.vars(y_labels) %in% names(data)) cli::cli_abort("y_labels must select one existing label column.")
    data$.vertical_label <- as.character(data[[all.vars(y_labels)]])
  }
  data$age <- if (categorical) match(vertical, unique(vertical)) else vertical
  if (!is.data.frame(data) || !all(c("year", "age", "est") %in% names(data)) ||
      !nrow(data)) cli::cli_abort("{.arg data} must contain year, age and est rows.")
  if (!all(vapply(data[c("year", "age", "est")], is.numeric, logical(1))) ||
      any(!is.finite(data$year)) || any(!is.finite(data$age))) {
    cli::cli_abort("year, age and est must be numeric; year and age must be finite.")
  }
  column <- function(form, default) {
    if (is.null(form)) return(rep(default, nrow(data)))
    if (!inherits(form, "formula") || length(form) != 2L ||
        length(all.vars(form)) != 1L || !is.symbol(form[[2L]]) ||
        !all.vars(form) %in% names(data)) {
      cli::cli_abort("metric and frame must each select one existing column.")
    }
    data[[all.vars(form)]]
  }
  metrics <- column(metric, "Estimate")
  models <- column(frame, "Model")
  if (anyNA(metrics) || anyNA(models)) cli::cli_abort("Metric and model labels cannot be missing.")
  labels <- function(x) if (is.factor(x)) levels(droplevels(x)) else unique(as.character(x))
  metric_names <- labels(metrics)
  model_names <- labels(models)
  keys <- data.frame(year = data$year, age = data$age, metric = metrics, model = models)
  if (anyDuplicated(keys)) cli::cli_abort("Duplicate year-age cells within a metric and model.")
  rows <- function(x) unname(lapply(seq_len(nrow(x)), function(i) unname(x[i, ])))
  span <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(0, 1))
    r <- range(x)
    if (diff(r) == 0) r <- r + c(-1, 1) * max(abs(r[1]) * .01, .01)
    unname(r)
  }
  spec <- lapply(metric_names, function(name) {
    d <- data[as.character(metrics) == name, , drop = FALSE]
    axis_label <- if (length(ylab) == 1L) unname(ylab) else unname(ylab[name])
    if (length(axis_label) != 1L || is.na(axis_label)) cli::cli_abort("ylab must supply a label for each metric.")
    units <- if ("unit" %in% names(d)) unique(as.character(d$unit[is.finite(d$est) &
      !is.na(d$unit) & nzchar(d$unit)])) else character()
    if (length(units) > 1L) cli::cli_abort("Metric {.val {name}} has different units; convert them before plotting.")
    unit <- if (length(units)) units[1L] else ""
    years <- sort(unique(d$year))
    ages <- sort(unique(d$age))
    grid <- expand.grid(age = ages, year = years)
    cells <- lapply(model_names, function(model) {
      z <- d[as.character(models[as.character(metrics) == name]) == model, , drop = FALSE]
      cell_ages <- ages
      ticktext <- as.character(ages)
      if (categorical) {
        ticks <- unique(z$.vertical_label)
        # Retain each model's own ordered state blocks rather than a union of ages.
        z$age <- match(z$.vertical_label, ticks)
        cell_ages <- seq_along(ticks)
        ticktext <- ticks
        if (!length(cell_ages)) { cell_ages <- 1; ticktext <- "Unavailable" }
      } else if (!is.null(y_labels)) {
        cell_ages <- sort(unique(z$age))
        ticktext <- z$.vertical_label[match(cell_ages, z$age)]
        if (!length(cell_ages)) { cell_ages <- 1; ticktext <- "" }
      }
      grid <- expand.grid(age = cell_ages, year = years)
      index <- match(paste(grid$year, grid$age), paste(z$year, z$age))
      values <- z$est[index]
      values[!is.finite(values)] <- NA_real_
      text <- paste0("Model: ", model, "<br>", name, ": ", format(values, trim = TRUE),
        if (nzchar(unit)) paste0(" ", unit) else "", "<br>Year: ", grid$year,
        if ("is_proj" %in% names(z)) ifelse(!is.na(z$is_proj[index]) & z$is_proj[index], " (projection)", "") else "",
        "<br>", axis_label, ": ",
        if (categorical || !is.null(y_labels)) ticktext[match(grid$age, cell_ages)] else grid$age,
        .forecast_hover(z)[index])
      for (field in intersect(c("definition", "notes"), names(z))) {
        value <- as.character(z[[field]][index])
        text <- paste0(text, ifelse(is.na(value) | !nzchar(value), "", paste0("<br>", value)))
      }
      natural <- matrix(values, nrow = length(cell_ages))
      logged <- natural
      logged[logged <= 0] <- NA_real_
      logged <- log10(logged)
      list(linear = rows(natural), log = rows(logged),
           ages = cell_ages, ticktext = ticktext,
           text = rows(matrix(text, nrow = length(cell_ages))))
    })
    list(name = name, unit = unit, years = years, ages = ages, cells = cells, ylab = axis_label,
      linear = if (signed) c(-1, 1) * max(c(abs(d$est[is.finite(d$est)]), .01)) else span(d$est),
      log = span(log10(d$est[is.finite(d$est) & d$est > 0])))
  })
  colours <- tam_cols()
  colorscale <- unname(Map(list, seq(0, 1, length.out = length(colours)), colours))
  if (signed) colorscale <- list(list(0, tail(tam_cols(), 1L)), list(.5, "#FFFFFF"), list(1, tam_cols()[1L]))
  config <- list(metrics = spec, models = model_names, type = type, colorscale = colorscale,
                 signed = signed, ylab = ylab, categorical = categorical || !is.null(y_labels))
  script <- paste(readLines(system.file("js", "matrix.js", package = "tinyAM"),
                           warn = FALSE), collapse = "\n")
  initial_metric <- if (signed) which(vapply(spec, function(metric) any(vapply(metric$cells,
    function(cell) any(is.finite(unlist(cell$linear))), logical(1))), logical(1)))[1L] else 1L
  if (is.na(initial_metric)) initial_metric <- 1L
  config$initial_metric <- initial_metric - 1L
  first <- spec[[initial_metric]]
  initial <- if (signed) which(vapply(first$cells, function(cell) any(is.finite(unlist(cell$linear))), logical(1)))[1L] else 1L
  if (is.na(initial)) initial <- 1L
  cell <- first$cells[[initial]]
  p <- plotly::plotly_build(plotly::plot_ly(x = first$years, y = cell$ages,
    z = matrix(0, nrow = length(cell$ages), ncol = length(first$years)),
    type = type))
  p$x$data[[1L]]$z <- cell$linear
  p$x$data[[1L]]$colorscale <- colorscale
  htmlwidgets::onRender(p, script, data = config)
}
