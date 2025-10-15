
#' Make a flexdashboard for visualizing model fits
#'
#' @param fits         A **named list** of fitted TAM objects. Required to be named;
#'                     the names are used as label values (even for length-1 lists).
#' @param output_file  Name of file to export using [rmarkdown::render()].
#'                     If `NULL`, flexdashboard will be rendered using [rmarkdown::run()]
#' @param ...          Additional arguments to send to [rmarkdown::run()] or
#'                     [rmarkdown::render()].
#'
#' @examples
#' if (interactive()) {
#'   N_dev <- fit_tam(
#'     cod_obs, years = 1983:2024, ages = 2:14,
#'     N_settings = list(process = "iid", init_N0 = FALSE),
#'     M_settings = list(process = "off", assumption = ~ I(0.3)),
#'     proj_settings = list(n_proj = 5, n_mean = 5, F_mult = 1),
#'     silent = TRUE
#'   )
#'   M_dev <- update(
#'     N_dev,
#'     N_settings = list(process = "off", init_N0 = TRUE),
#'     M_settings = list(process = "ar1", assumption = ~ I(0.3),
#'                       age_breaks = c(2, 14))
#'   )
#'
#'   fits <- list("N_dev" = N_dev, "M_dev" = M_dev)
#'   vis_tam(fits)
#' }
#'
#' @export
vis_tam <- function(fits = NULL, output_file = NULL, ...) {

  pkg <- c("rmarkdown", "flexdashboard")
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

  rmd_env <- new.env()
  rmd_env$fits <- fits

  if (is.null(output_file)) {
    rmarkdown::run(file = rmd_file,
                   render_args = list(envir = rmd_env), ...)
  } else {
    output_dir <- normalizePath(dirname(output_file))
    output_file <- basename(output_file)
    rmarkdown::render(input = rmd_file,
                      output_file = output_file,
                      output_dir = output_dir, ...)
  }

}
