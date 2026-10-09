
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
#'          Random plots include latent states and formula effects; RW changes
#'          and effects multiplied by numeric covariates are explained alongside
#'          their plots.
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

.fixed_parameter_groups <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) return(list())
  par <- data$par
  category <- rep("Other", nrow(data))
  category[grepl("^sd_", par)] <- "Process SDs"
  category[par %in% c("sd_catch", "sd_index")] <- "Observation SDs"
  category[par %in% c("q", "logit_q")] <- "Catchability"
  category[par %in% "dq" | grepl("^q_(a50|slope)_", par)] <- "Catchability curves"
  category[grepl("^mu_", par)] <- "Mortality means"
  category[grepl("^phi_", par)] <- "Correlations"
  category[par %in% c("r0", "n0")] <- "Initial abundance"
  category <- factor(category, levels = c("Process SDs", "Observation SDs",
    "Catchability", "Catchability curves", "Mortality means", "Correlations",
    "Initial abundance", "Other"))
  split(data, category, drop = TRUE)
}
