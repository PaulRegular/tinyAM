
#' Check the structure of a TAM obs list
#'
#' @description
#' Validates the minimal structure required by [make_dat()] and [fit_tam()]:
#' the `obs` object must be a list containing data.frames `catch`, `index`,
#' `weight`, and `maturity`, each with columns `year`, `age`, and `obs`.
#' Types must be numeric/integerish for `year` and `age`, and numeric for `obs`.
#' If present, `index$samp_time` must be numeric, between 0 to 1. Duplicate
#' `(year, age)` rows are warned.
#'
#' On failure, the function aborts with a `cli` error. On success, it returns
#' `TRUE` (invisibly).
#'
#' @param obs A named list expected to contain `catch`, `index`, `weight`,
#'   and `maturity` data.frames with the columns described above.
#'
#' @return `TRUE` (invisibly) if validation passes; otherwise aborts.
#' @examples
#' data(northern_cod_data)
#' check_obs(northern_cod_data)  # returns TRUE (invisibly) if valid
#'
#' @importFrom cli cli_abort cli_warn
#' @export
check_obs <- function(obs) {
  # high-level shape
  if (!is.list(obs)) {
    cli::cli_abort(c("x" = "`obs` must be a list."))
  }

  required_tabs <- c("catch", "index", "weight", "maturity")
  missing_tabs  <- setdiff(required_tabs, names(obs))
  if (length(missing_tabs)) {
    cli::cli_abort(c(
      "{.strong Invalid obs list}",
      "x" = "Missing required tables: {paste(missing_tabs, collapse = \", \")}",
      "i" = "Expected names: {paste(required_tabs, collapse = \", \")}"
    ))
  }

  # per-table checks
  check_one <- function(x, nm) {
    if (!is.data.frame(x)) {
      cli::cli_abort(c(
        "{.strong Invalid obs ${nm}}",
        "x" = "`{nm}` must be a data.frame (got {class(x)[1]})."
      ))
    }
    need <- c("year", "age", "obs")
    miss <- setdiff(need, names(x))
    if (length(miss)) {
      cli::cli_abort(c(
        "{.strong Invalid obs ${nm}}",
        "x" = "`{nm}` is missing required column(s): {paste(miss, collapse = \", \")}"
      ))
    }
    # type checks (allow integerish for year/age; obs must be numeric)
    if (!is.numeric(x$year)) {
      cli::cli_abort(c("{.strong Invalid obs ${nm}}", "x" = "`{nm}$year` must be numeric/integer."))
    }
    if (!is.numeric(x$age)) {
      cli::cli_abort(c("{.strong Invalid obs ${nm}}", "x" = "`{nm}$age` must be numeric/integer."))
    }
    if (!is.numeric(x$obs)) {
      cli::cli_abort(c("{.strong Invalid obs ${nm}}", "x" = "`{nm}$obs` must be numeric (NA allowed)."))
    }
    # duplicate (year, age) rows: warn
    if (any(duplicated(x[c("year", "age")]))) {
      cli::cli_warn(c(
        "!" = "`{nm}` has duplicate (year, age) rows; downstream joins may be ambiguous."
      ))
    }
    # index-specific soft check for samp_time
    if (nm == "index" && "samp_time" %in% names(x)) {
      if (!is.numeric(x$samp_time) || any(x$samp_time < 0 | x$samp_time > 1, na.rm = TRUE)) {
        cli::cli_abort(c(
          "{.strong Invalid obs index}",
          "x" = "`index$samp_time` must be numeric in [0, 1] (NA allowed)."
        ))
      }
    } else if (nm == "index" && !("samp_time" %in% names(x))) {
      cli::cli_warn(c("!" = "`index$samp_time` is missing; model will treat it as NA unless supplied."))
    }
    invisible(TRUE)
  }

  # run checks
  for (nm in required_tabs) check_one(obs[[nm]], nm)

  invisible(TRUE)
}
