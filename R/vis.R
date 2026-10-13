
#' Make a flexdashboard for visualizing model fits
#'
#' @inheritParams tidy_tam
#' @param model_list   A **named list** of fitted TAM objects or precomputed
#'                     reference objects of class `tam_ref`. Names label models.
#' @param output_file  Name of file to export using [rmarkdown::render()].
#'                     If `NULL`, a temporary HTML file is rendered.
#'                     The file opens in your browser only when `open_file = TRUE`.
#' @param open_file    Logical. Open rendered html file?
#' @param background Optional Markdown text shown on a Background dashboard page.
#' @param render_args  Named list of additional arguments passed to
#'                     [rmarkdown::render()].
#' @param ...          One or more TAM fits or reporting references.
#'                     Supply these or `model_list`, not both. When supplying models
#'                     through `...`, their object names are used to label models
#'                     (even when a single model is supplied).
#' @details Supply models via `...` or `model_list`, but not both. The models must
#'          form a uniquely named list so they can be labeled in the dashboard.
#'          Additional arguments for [rmarkdown::render()] can be passed through
#'          `render_args`, which must itself be a (named) list.
#'          Precomputed assessment references may be mixed with tinyAM fits; they are
#'          not refittable TAM objects. Supply `background` to add a page
#'          describing assessment assumptions and translation choices.
#'          The Parameters menu separates Fixed and Random pages. Fixed plots
#'          group similar quantities, with initial abundance shown separately.
#'          Stock-recruit models add a relationship page with parent-year pairs,
#'          median curves at declared reference covariates and conditional
#'          recruitment predictions. Random parameters also show signed
#'          recruitment residuals. See [tidy_recruitment()].
#'          Random plots include latent states and formula-effect trends.
#'          Age-year plots offer Annual trends (year axis, age selector) and
#'          Age profiles (age axis, year selector). Each plot starts in its
#'          specified orientation and retains each view's controls when switching.
#'          F-at-age and M-at-age include available process means and their
#'          separate confidence ribbons; Average F remains the population summary.
#' @return Used for its side effects: writes an HTML dashboard and optionally
#'   opens it in the browser. Supply `output_file` to retain a known file path.
#'
#' @example inst/examples/example_fits.R
#' @examples
#' if (interactive()) {
#'   ## Build dashboard ----
#'   vis_tam(fits)
#' }
#'
#'
#' @importFrom utils browseURL
#'
#' @export
vis_tam <- function(..., model_list = NULL, interval = 0.95, output_file = NULL,
                    open_file = TRUE, background = NULL, render_args = list()) {
  pkg <- c("knitr", "rmarkdown", "flexdashboard")
  missing <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install_call <- sprintf(
      "install.packages(c(%s))",
      paste(sprintf("'%s'", missing), collapse = ", ")
    )

    cli::cli_abort(c(
      "Required package(s) not installed: {paste(missing, collapse = ', ')}.",
      "i" = "Install with: {.code {install_call}}"
    ))
  }

  if (!is.null(background)) {
    if (!is.character(background) || !length(background) || anyNA(background)) {
      cli::cli_abort("{.arg background} must be NULL or Markdown text.")
    }
    unmarked_utf8 <- Encoding(background) == "unknown" & validUTF8(background)
    Encoding(background)[unmarked_utf8] <- "UTF-8"
    background <- paste(enc2utf8(background), collapse = "\n")
    if (!nzchar(trimws(background))) background <- NULL
  }

  render_args <- .validate_named_list(render_args, arg = "render_args", allow_empty = TRUE,
                                      require_unique = FALSE)

  dots <- list(...)
  dot_expr <- as.list(substitute(list(...)))[-1]

  fits_info <- .dots_or_list(dots, dot_expr, model_list = model_list, list_arg_name = "model_list")
  fits <- fits_info$fits

  rmd_file <- system.file("rmd", "vis_tam.Rmd", package = "tinyAM")
  rmd_env <- new.env(parent = globalenv())
  rmd_env$fits <- fits
  rmd_env$interval <- interval
  rmd_env$background <- background

  if (is.null(output_file)) {
    output_file <- tempfile(pattern = "vis_tam_", fileext = ".html")
  }
  output_dir <- normalizePath(dirname(output_file))
  output_name <- basename(output_file)
  render_call <- c(
    list(
      input = rmd_file,
      output_file = output_name,
      output_dir = output_dir,
      envir = rmd_env
    ),
    render_args
  )
  do.call(rmarkdown::render, render_call)

  if (open_file) utils::browseURL(output_file)

}

# Two independent widgets retain their animation and control state when hidden.
.plot_age_year <- function(data, plot = plot_trend, x = ~year, xlab = "Year",
                           frame = NULL, ...) {
  args <- list(...)
  eligible <- all(c("year", "age") %in% names(data)) &&
    ((identical(all.vars(x), "year") && identical(all.vars(frame), "age")) ||
       (identical(all.vars(x), "age") && identical(all.vars(frame), "year")))
  if (!eligible || !nrow(data)) {
    return(do.call(plot, c(list(data = data, x = x, xlab = xlab, frame = frame), args)))
  }
  data <- data[order(data$year, data$age), , drop = FALSE]
  projection <- if ("is_proj" %in% names(data)) data$is_proj else rep(FALSE, nrow(data))
  years <- unique(data$year)
  suffix <- vapply(years, function(year) {
    flags <- projection[data$year == year]
    if (all(flags)) " (projection)" else if (any(flags)) " (includes projections)" else ""
  }, character(1))
  year_labels <- paste0(data$year, suffix[match(data$year, years)])
  data$.view_year <- factor(year_labels, levels = unique(year_labels))
  if (is.null(args$text)) {
    data$.view_hover <- paste0("Year: ", data$year,
      ifelse(projection, " (projection)", ""), "<br>Age: ", data$age)
    args$text <- ~.view_hover
  }
  make_view <- function(x, xlab, frame, prefix, view) {
    animated <- length(unique(data[[all.vars(frame)]])) > 1L
    if (!animated) frame <- NULL
    p <- do.call(plot, c(list(data = data, x = x, xlab = xlab, frame = frame), args))
    if (animated) p <- plotly::animation_slider(p,
      currentvalue = list(prefix = prefix), pad = list(t = 45))
    p <- .age_year_menu(p, view)
    p$width <- "100%"
    p$height <- "100%"
    p
  }
  annual <- make_view(~year, "Year", ~age, "Age: ", "annual")
  profile <- make_view(~age, "Age", ~.view_year, "Year: ", "profile")
  initial <- if (identical(all.vars(frame), "age")) "annual" else "profile"
  panel <- function(view, widget) htmltools::tags$div(
    class = "tam-age-year-panel", `data-view` = view,
    style = paste0("flex:1;min-height:0;", if (view != initial) "display:none;"), widget)
  script <- paste(readLines(system.file("js", "age-year.js", package = "tinyAM"),
                           warn = FALSE), collapse = "\n")
  htmltools::browsable(htmltools::tags$div(
    class = "tam-age-year", `data-default-view` = initial,
    style = "display:flex;flex-direction:column;width:100%;height:100%;min-height:400px;",
    panel("annual", annual), panel("profile", profile),
    htmltools::tags$script(htmltools::HTML(script))))
}

.age_year_menu <- function(p, view) {
  menus <- plotly::plotly_build(p)$x$layout$updatemenus
  controls <- Filter(function(menu) any(unlist(lapply(menu$buttons, `[[`, "label")) %in%
    c("Linear", "Log", "Show CI", "Hide CI")), menus)
  menu <- if (length(controls)) controls[[length(controls)]] else .log_linear_buttons[[1L]]
  menu$name <- "tam-age-year-view"
  menu$pad$t <- if (is.null(menu$pad$t)) 75 else menu$pad$t + 75
  menu$active <- match(view, c("annual", "profile")) - 1L
  menu$buttons <- list(
    list(method = "skip", execute = FALSE, args = list("annual"), label = "Annual trends"),
    list(method = "skip", execute = FALSE, args = list("profile"), label = "Age profiles"))
  if (is.null(p$x$layout)) p$x$layout <- list()
  p$x$layout$updatemenus <- c(menus, list(menu))
  p
}

.mortality_plot_data <- function(data, mean, fits, component) {
  for (nm in c("est", "lwr", "upr")) data[[paste0("mean_", nm)]] <- NA_real_
  if (is.null(mean) || !nrow(mean)) return(data)
  if (component == "F") {
    no_mean <- names(fits)[vapply(fits, function(fit) {
      inherits(fit, "tam_fit") && is.null(fit$dat$F_settings$mu_form)
    }, logical(1))]
    mean <- mean[!mean$model %in% no_mean, , drop = FALSE]
  }
  keys <- c("model", "year", "age", "is_proj")
  key <- function(d) do.call(paste, c(d[keys], sep = "\r"))
  rows <- match(key(data), key(mean))
  for (nm in intersect(c("est", "lwr", "upr"), names(mean))) {
    data[[paste0("mean_", nm)]] <- mean[[nm]][rows]
  }
  data
}

.plot_mortality <- function(data, ...) {
  for (nm in c("lwr", "upr")) {
    if (is.null(data[[nm]])) data[[nm]] <- NA_real_
  }
  p <- plot_trend(data, add_buttons = FALSE, ...)
  if (any(is.finite(data$mean_est))) {
    if (any(is.finite(data$mean_upr))) {
      p <- plotly::add_ribbons(p, ymin = ~mean_lwr, ymax = ~mean_upr,
        opacity = .15, line = list(width = 0), showlegend = FALSE,
        hoverinfo = "skip")
    }
    # Missing reference means must not drop or shift animation traces.
    data$mean_est[is.na(data$mean_est)] <- Inf
    p <- plotly::add_lines(p, data = data, y = ~mean_est, showlegend = FALSE,
      line = list(width = .5),
      hovertemplate = "Process mean: %{y}<br>%{text}<extra>%{fullData.name}</extra>")
    values <- unlist(data[c("est", "upr", "mean_est", "mean_upr")])
    p <- plotly::layout(p, yaxis = list(range = c(0, max(c(0, values[is.finite(values)])) * 1.05)))
  }
  .plot_trend_buttons(p)
}

.fixed_parameter_groups <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) return(list())
  par <- data$par
  category <- rep("Other", nrow(data))
  category[grepl("^sd_", par)] <- "Process SDs"
  category[grepl("^sd_beta_", par)] <- "Process SD coefficients"
  category[par %in% c("sd_catch", "sd_index")] <- "Observation SDs"
  category[par %in% c("q", "logit_q")] <- "Catchability"
  category[par %in% "dq" | grepl("^q_(a50|slope)_", par)] <- "Catchability curves"
  category[grepl("^mu_", par)] <- "Mortality means"
  category[par == "dmu_f" | grepl("^F_(a50|slope)_", par)] <- "Mortality means"
  category[grepl("^(phi_|rho_)", par)] <- "Correlations"
  category[par %in% c("rec_beta", "sr_alpha", "sr_beta")] <- "Recruitment"
  category[par %in% c("r0", "r_init", "n0")] <- "Initial abundance"
  category <- factor(category, levels = c("Process SDs", "Process SD coefficients", "Observation SDs",
    "Catchability", "Catchability curves", "Mortality means", "Correlations",
    "Recruitment", "Initial abundance", "Other"))
  split(data, category, drop = TRUE)
}

.formula_effect_label <- function(effect) {
  component <- unique(effect$component)
  if (component == "q") component <- "Catchability"
  if (component == "sd_catch") component <- "Catch SD"
  if (component == "sd_index") component <- "Survey SD"
  process <- switch(unique(effect$process), rw = "random walk", ar1 = "AR1", iid = "IID")
  variable <- unique(effect$variable)
  by <- unique(effect$by)
  if (!is.na(by)) variable <- paste(variable, "by", by)
  paste0(component, " ", process, " (", gsub("_", " ", variable, fixed = TRUE), ")")
}
