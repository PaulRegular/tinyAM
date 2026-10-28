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
#' @param model Exact model name in a `tam_retro` result. A single model is
#'   selected automatically; interactive sessions offer a chooser for several
#'   models. Scripts must specify a model when several are available.
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
#'          recruitment deviations. See [tidy_recruitment()].
#'          The Random page separates Recruitment, N, F and M states from
#'          Process deviations and Joint-sample residuals. Optional initial
#'          abundance plots are inside N, followed by optional formula-effect tabs.
#'          [plot_process_deviations()] and [plot_process_residuals()] reproduce
#'          individual diagnostic plots from a single tidy result, selecting
#'          `component = "Recruitment"` or `"NFM"`. The dashboard supplies dropdown
#'          tabs for these views; N/F/M plots have component buttons and a model slider.
#'          Shared M blocks are not treated as separate biological ages.
#'          `process_seed` reproduces the same Gaussian joint draw as [tidy_tam()].
#'          Random plots include latent states and formula-effect trends.
#'          Age-year plots offer Annual trends (year axis, age selector) and
#'          Age profiles (age axis, year selector). Each plot starts in its
#'          specified orientation and retains each view's controls when switching.
#'          F-at-age and M-at-age include available process means and their
#'          separate confidence ribbons; Average F remains the population summary.
#'          Reproduce ordinary dashboard plots with [plot_trend()] and
#'          [plot_obs_pred()], or use [plot_age_year()] for the two orientations.
#'          The small F/M mean-overlay helpers and their line styling live in
#'          `inst/rmd/vis_tam.Rmd`.
#'          A `tam_retro` dashboard compares the selected model's successful
#'          peels with its unpeeled reference. Failed or missing attempts are
#'          listed on Background. Hindcast forecasts remain visible, with origin
#'          and horizon in hover text. A common projection separator is omitted
#'          when origins differ. The result is never refitted or modified.
#'          Population patterns and Mortality patterns use [plot_matrix()] to
#'          explore age-year heatmaps, surfaces and contours, with shared ranges
#'          across models for each metric. Log changes values, not year or age.
#'          Hover text retains natural units, spawning definitions and projection
#'          labels. Unknown reference units and unavailable cells remain missing.
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
#' @importFrom utils browseURL
#'
#' @export
vis_tam <- function(..., model_list = NULL, interval = 0.95, output_file = NULL,
                    open_file = TRUE, background = NULL, render_args = list(), model = NULL,
                    process_seed = 123456L) {
  first <- if (length(list(...))) list(...)[[1L]] else NULL
  UseMethod("vis_tam", first)
}

#' @rdname vis_tam
#' @export
vis_tam.default <- function(..., model_list = NULL, interval = 0.95, output_file = NULL,
                            open_file = TRUE, background = NULL, render_args = list(), model = NULL,
                            process_seed = 123456L) {
  dots <- list(...)
  dot_expr <- as.list(substitute(list(...)))[-1]
  if (any(vapply(c(dots, model_list), inherits, logical(1), "tam_retro"))) {
    cli::cli_abort("Supply one tam_retro object alone, without additional fits or model_list.")
  }
  if (!is.null(model)) cli::cli_abort("model is only used with a tam_retro result.")
  fits <- .dots_or_list(dots, dot_expr, model_list = model_list, list_arg_name = "model_list")$fits
  .render_vis_tam(fits, interval, output_file, open_file, background, render_args,
    process_seed = process_seed)
}

.retro_dashboard <- function(retro, model = NULL, is_interactive = interactive(),
                             choose = utils::select.list) {
  models <- names(retro$fits)
  if (is.null(model)) {
    if (length(models) == 1L) {
      model <- models
    } else {
      if (!is_interactive) cli::cli_abort(c(
        "Specify model: retrospective dashboards display one model at a time.",
        "i" = "Available models: {paste(models, collapse = ', ')}."
      ))
      cli::cli_inform("The dashboard displays one model at a time; select a model.")
      model <- choose(models, title = "Retrospective model")
      if (!length(model) || !nzchar(model)) return(NULL)
    }
  }
  if (!is.character(model) || length(model) != 1L || is.na(model) || !model %in% models) {
    cli::cli_abort("model must be one exact name: {paste(models, collapse = ', ')}.")
  }
  diagnostics <- retro$fold_diagnostics[retro$fold_diagnostics$model == model, , drop = FALSE]
  diagnostics <- diagnostics[order(!diagnostics$is_reference, -diagnostics$fold), , drop = FALSE]
  candidates <- retro$fits[[model]][as.character(diagnostics$fold)]
  keep <- !is.na(diagnostics$is_converged) & diagnostics$is_converged &
    vapply(candidates, inherits, logical(1), "tam_fit")
  if (!any(keep)) cli::cli_abort(c(
    "No successful fits are available for {model}.",
    "i" = "Inspect fold_diagnostics for errors and numerical failures."
  ))
  labels <- as.character(diagnostics$fold[keep])
  fits <- stats::setNames(candidates[keep], labels)
  origins <- stats::setNames(diagnostics$fold[keep], labels)
  omitted <- diagnostics[!keep, , drop = FALSE]
  background <- c("### Retrospective comparison", "",
    paste0("Model: **", model, "**. Traces show the full-data reference and successful peels."),
    paste0("Omitted attempts: **", nrow(omitted), "**."))
  if (!any(diagnostics$is_reference & keep)) background <- c(background, "The full-data reference is unavailable.")
  if (any(diagnostics$is_reference & keep)) background <- c(background,
    paste0("Full-data reference year: **", diagnostics$fold[diagnostics$is_reference & keep], "**."))
  if (nrow(omitted)) {
    reason <- omitted$error
    reason[is.na(reason) | !nzchar(reason)] <- "Missing fit or failed numerical checks; see fold_diagnostics."
    background <- c(background, "", paste0("- ", omitted$fold, ": ", gsub("[\r\n]+", " ", reason)))
  }
  list(fits = fits, origins = origins, background = background)
}

#' @rdname vis_tam
#' @export
vis_tam.tam_retro <- function(..., model_list = NULL, interval = 0.95, output_file = NULL,
                              open_file = TRUE, background = NULL, render_args = list(), model = NULL,
                              process_seed = 123456L) {
  dots <- list(...)
  if (length(dots) != 1L || !is.null(model_list)) {
    cli::cli_abort("Supply one tam_retro object alone, without additional fits or model_list.")
  }
  view <- .retro_dashboard(dots[[1L]], model)
  if (is.null(view)) return(invisible(NULL))
  .render_vis_tam(view$fits, interval, output_file, open_file,
    c(view$background, "", background), render_args, retro_origins = view$origins,
    process_seed = process_seed)
}

.retro_plot_metadata <- function(x, origins) {
  if (is.data.frame(x)) {
    if (all(c("model", "year") %in% names(x))) {
      x$fold <- unname(origins[as.character(x$model)])
      x$horizon <- pmax(0, x$year - x$fold)
    }
  } else if (is.list(x)) x[] <- lapply(x, .retro_plot_metadata, origins = origins)
  x
}

.render_vis_tam <- function(fits, interval, output_file, open_file, background,
                            render_args, retro_origins = NULL, process_seed = 123456L) {
  process_seed <- .validate_process_seed(process_seed)
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

  render_args <- .validate_named_list(render_args,
    arg = "render_args", allow_empty = TRUE,
    require_unique = FALSE
  )

  rmd_file <- system.file("rmd", "vis_tam.Rmd", package = "tinyAM")
  rmd_env <- new.env(parent = globalenv())
  rmd_env$fits <- fits
  rmd_env$interval <- interval
  rmd_env$process_seed <- process_seed
  rmd_env$background <- background
  rmd_env$retro_origins <- retro_origins

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

.fixed_parameter_groups <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) {
    return(list())
  }
  par <- data$par
  category <- rep("Other", nrow(data))
  category[grepl("^sd_", par)] <- "Process SDs"
  category[grepl("^sd_beta_", par)] <- "Process SD coefficients"
  category[par %in% c("sd_catch", "sd_index")] <- "Observation SDs"
  category[par %in% c("q", "logit_q", "q_power")] <- "Catchability"
  category[par %in% "dq" | grepl("^q_(a50|slope)_", par)] <- "Catchability curves"
  category[grepl("^mu_", par)] <- "Mortality means"
  category[par == "mu_n"] <- "Cohort means"
  category[par == "dmu_f" | grepl("^F_(a50|slope)_", par)] <- "Mortality means"
  category[grepl("^(phi_|rho_)", par)] <- "Correlations"
  category[par %in% c("rec_beta", "sr_alpha", "sr_beta")] <- "Recruitment"
  category[par %in% c("r0", "r_init", "n0")] <- "Initial abundance"
  category[par == "d_n0"] <- "Initial steps"
  category <- factor(category, levels = c(
    "Process SDs", "Process SD coefficients", "Observation SDs",
    "Catchability", "Catchability curves", "Mortality means", "Cohort means", "Correlations",
    "Recruitment", "Initial abundance", "Initial steps", "Other"
  ))
  split(data, category, drop = TRUE)
}

.formula_effect_label <- function(effect) {
  component <- unique(effect$component)
  if (component == "q") component <- "Catchability"
  if (component == "sd_catch") component <- "Catch SD"
  if (component == "sd_index") component <- "Survey SD"
  process <- switch(unique(effect$process),
    rw = "random walk",
    ar1 = "AR1",
    iid = "IID"
  )
  variable <- unique(effect$variable)
  by <- unique(effect$by)
  if (!is.na(by)) variable <- paste(variable, "by", by)
  paste0(component, " ", process, " (", gsub("_", " ", variable, fixed = TRUE), ")")
}
