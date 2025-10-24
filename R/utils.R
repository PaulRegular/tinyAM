#' Validate a (possibly empty) named list input
#'
#' @param x Object to validate.
#' @param arg Character label used in error messages.
#' @param allow_empty Logical; should zero-length lists be allowed?
#' @param require_unique Logical; should the names be required to be unique?
#'
#' @return The validated list.
#'
#' @noRd
.validate_named_list <- function(x, arg, allow_empty = FALSE, require_unique = TRUE) {
  if (!is.list(x)) {
    cli::cli_abort("{.arg {arg}} must be a list.")
  }

  if (!length(x)) {
    if (allow_empty) {
      return(x)
    }
    cli::cli_abort("{.arg {arg}} must be a non-empty list.")
  }

  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    cli::cli_abort("{.arg {arg}} must be a named list.")
  }

  if (require_unique && anyDuplicated(nm)) {
    cli::cli_abort("{.arg {arg}} names must be unique.")
  }

  x
}

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
#' @details Ensures that the resulting list is non-empty (unless the explicit
#'   list argument is `NULL`), uniquely named, and that every element resembles a
#'   TAM fit (i.e. contains components `dat`, `opt`, and `pop`).
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

  arg_label <- if (using_dots) "models" else list_arg_name
  fits <- .validate_named_list(fits, arg = arg_label)

  is_fit <- vapply(fits, .is_tam_fit, logical(1))
  if (any(!is_fit)) {
    bad <- names(fits)[!is_fit]
    cli::cli_abort(c(
      "All supplied models must be {.cls tam_fit} objects produced by {.fn fit_tam}.",
      "x" = "Problematic element{?s}: {cli::format_inline('{.val {bad}}')}"
    ))
  }

  list(fits = fits, using_dots = using_dots)
}
