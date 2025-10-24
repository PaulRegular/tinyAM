
#' Make a flexdashboard for visualizing model fits
#'
#' @inheritParams tidy_tam
#' @param model_list   A **named list** of fitted TAM objects (e.g., `fits` list returned by
#'                     [fit_retro()]). Names are used to label models.
#' @param output_file  Name of file to export using [rmarkdown::render()].
#'                     If `NULL`, a temporary HTML file will be rendered and opened
#'                     automatically in your default browser.
#' @param open_file    Logical. Open rendered html file?
#' @param render_args  Named list of additional arguments passed to
#'                     [rmarkdown::render()].
#' @param ...          One or more fitted TAM objects (as returned by [fit_tam()]).
#'                     Ignored if `model_list` is provided. When supplying models
#'                     through `...`, their object names are used to label models
#'                     (even when a single model is supplied).
#' @details Supply models via `...` or `model_list`, but not both. The models must
#'          form a uniquely named list so they can be labeled in the dashboard.
#'          Additional arguments for [rmarkdown::render()] can be passed through
#'          `render_args`, which must itself be a (named) list.
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
                    open_file = TRUE, render_args = list()) {

  pkg <- c("knitr", "rmarkdown", "flexdashboard")
  missing <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install_call <- sprintf(
      "install.packages(c(%s))",
      paste(sprintf("'%s'", missing), collapse = ", ")
    )
    cli::cli_abort(c(
      "Package{?s} {.pkg {missing}} {cli::qty(length(missing))}is/are required for {.fn vis_tam}.",
      "i" = "Install with: {.code {install_call}}"
    ))
  }

  if (!is.list(render_args)) {
    cli::cli_abort("`render_args` must be a list.")
  }

  if (length(render_args) && (is.null(names(render_args)) || any(!nzchar(names(render_args))))) {
    cli::cli_abort("`render_args` must be a named list.")
  }

  dots <- list(...)
  dot_expr <- as.list(substitute(list(...)))[-1]

  if (!is.null(model_list) && length(dots)) {
    cli::cli_abort("Supply models either through `...` or `model_list`, not both.")
  }

  if (!is.null(model_list)) {
    fits <- model_list
  } else if (length(dots)) {
    is_tam_fit <- function(x) {
      is.list(x) && !is.null(names(x)) && all(c("dat", "opt", "pop") %in% names(x))
    }

    if (length(dots) == 1L && is.list(dots[[1]]) && !is_tam_fit(dots[[1]])) {
      fits <- dots[[1]]
    } else {
      current_names <- names(dots)
      if (is.null(current_names)) {
        current_names <- rep("", length(dots))
      }
      dot_names <- vapply(dot_expr, deparse1, character(1))
      missing <- !nzchar(current_names)
      current_names[missing] <- dot_names[missing]
      names(dots) <- current_names
      fits <- dots
    }
  } else {
    cli::cli_abort("No models supplied.")
  }

  if (!is.list(fits) || !length(fits)) {
    cli::cli_abort("Supplied models must be provided as a non-empty list of fits.")
  }

  if (is.null(names(fits)) || any(!nzchar(names(fits)))) {
    cli::cli_abort("Supplied models must form a named list.")
  }

  if (anyDuplicated(names(fits))) {
    cli::cli_abort("Model names must be unique.")
  }

  rmd_file <- system.file("rmd", "vis_tam.Rmd", package = "tinyAM")
  rmd_env <- new.env(parent = globalenv())
  rmd_env$fits <- fits
  rmd_env$interval <- interval

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
