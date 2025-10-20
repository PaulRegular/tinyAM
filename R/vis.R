
#' Make a flexdashboard for visualizing model fits
#'
#' @param fits         A **named list** of fitted TAM objects (e.g., `fits` list returned by
#'                     [fit_retro()]). Names are used to label models.
#' @param output_file  Name of file to export using [rmarkdown::render()].
#'                     If `NULL`, a temporary HTML file will be rendered and opened
#'                     automatically in your default browser.
#' @param open_file    Logical. Open rendered html file?
#' @param ...          Additional arguments to send to [rmarkdown::render()].
#'
#' @example inst/examples/example_fits.R
#' @examples
#' if (interactive) {
#'   ## Build dashboard ----
#'   vis_tam(fits)
#' }
#'
#'
#' @importFrom utils browseURL
#'
#' @export
vis_tam <- function(fits = NULL, output_file = NULL, open_file = TRUE, ...) {

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

  rmd_file <- system.file("rmd", "vis_tam.Rmd", package = "tinyAM")
  rmd_env <- new.env(parent = globalenv())
  rmd_env$fits <- fits

  if (is.null(output_file)) {
    output_file <- tempfile(pattern = "vis_tam_", fileext = ".html")
  }
  output_dir <- normalizePath(dirname(output_file))
  output_name <- basename(output_file)
  rmarkdown::render(input = rmd_file,
                    output_file = output_name,
                    output_dir = output_dir,
                    envir = rmd_env,
                    ...)

  if (open_file) utils::browseURL(output_file)

}
