#' Normalize TAM model inputs supplied via `...` or `model_list`
#'
#' @param dots A list captured from `...` inside a calling function.
#' @param dot_expr The unevaluated expressions supplied through `...`.
#' @param model_list Optional explicit list argument passed to the calling
#'   function.
#' @param list_arg_name Character name of the explicit list argument, used in
#'   error messages.
#'
#' @return A list with elements `fits` (the validated, named list of models) and
#'   `using_dots` (logical indicating whether the models originated from
#'   `...`).
#'
#' @noRd
.dots_or_list <- function(dots, dot_expr, model_list = NULL, list_arg_name = "model_list") {
  using_dots <- is.null(model_list)

  if (!using_dots && length(dots)) {
    cli::cli_abort("Supply models either through `...` or `{list_arg_name}`, not both.")
  }

  if (using_dots) {
    if (!length(dots)) {
      cli::cli_abort("No models supplied.")
    }

    if (length(dots) == 1L && is.list(dots[[1]]) && !.is_tam_fit(dots[[1]])) {
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
    fits <- model_list
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

  list(fits = fits, using_dots = using_dots)
}

.is_tam_fit <- function(x) {
  is.list(x) && !is.null(names(x)) && all(c("dat", "opt", "pop") %in% names(x))
}

