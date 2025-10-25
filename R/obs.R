
#' Northern cod example data for TAM
#'
#' @description
#' A compact, tidy-ish list of inputs prepared for the Tiny Assessment Model (TAM),
#' derived from data used in the assessment of Northern cod, NAFO Divisions 2J3KL.
#' This is intended for vignettes, examples, and tests.
#'
#' @format
#' A named list of four data frames:
#'
#' - **catch**: commercial catch-at-age observations and blocking variables.
#'   Columns:
#'   - `year`: numeric calendar year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric observed catch in numbers (thousands; may include `NA`).
#'   - `F_y_block`: factor; year blocks for fishing-mortality mean structure
#'     (e.g., `"1954-1991"`, `"1992-1997"`, ...).
#'   - `F_a_block`: factor; age blocks for fishing-mortality mean structure
#'     (e.g., `"0-1"`, `"2-3"`, ...).
#'
#' - **index**: survey (index) observations and blocking variables.
#'   Columns:
#'   - `year`: integer year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric survey index-at-age (may include `NA`).
#'   - `survey`: character or factor; survey name for catchability blocks and labeling plots.
#'   - `samp_time`: numeric sampling time within the year (fraction of a year; e.g., `0.8`).
#'   - `q_block`: factor; age blocks for catchability.
#'
#' - **weight**: weight-at-age and mortality assumption used for \eqn{M} mean structure.
#'   Columns:
#'   - `year`: integer year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric stock weight-at-age (kg).
#'   - `M_assumption`: numeric assumed natural mortality level by age-year
#'     (used when constructing `M` mean via formula).
#'
#' - **maturity**: maturity-at-age (proportion mature).
#'   Columns:
#'   - `year`: numeric year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric proportion mature (0-1).
#'
#' @details
#' The object was assembled by:
#'
#' 1. Converting NCAM arrays to tidy data frames.
#' 2. Creating blocking variables for `F` by year (`F_y_block`) and age
#'    (`F_a_block`) using [cut_years()] and [cut_ages()].
#' 3. Joining weight-at-age with a natural mortality assumption column
#'    (`M_assumption`) derived from `NCAM::inputs$nm`.
#'
#' Exact dimensions may change if the source data evolve. At the time of inclusion,
#' the elements contained approximately:
#'
#' - `catch`: 1,065 rows x 6 columns
#' - `index`: 1,065 rows x 6 columns
#' - `weight`: 1,065 rows x 4 columns
#' - `maturity`: 1,215 rows x 3 columns
#'
#' @source
#' Derived from the \pkg{NCAM} example inputs
#' (`NCAM::inputs`). See the NCAM package for provenance and licensing.
#'
#' @examples
#' data(cod_obs)
#' names(cod_obs)
#'
#' # Peek at catch blocks
#' head(cod_obs$catch)
#' with(cod_obs$catch, table(F_y_block, F_a_block))
#'
#' # Build a minimal TAM dataset
#' dat <- make_dat(
#'   obs = cod_obs,
#'   years = 1983:2024,
#'   ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", mu_form = NULL,
#'                     assumption = ~ I(0.2), age_breaks = NULL),
#'   obs_settings = list(sd_catch_form = ~ 1, sd_index_form = ~ 1, q_form = ~ q_block)
#' )
#'
#' @docType data
#' @keywords datasets
#' @name cod_obs
#' @usage data(cod_obs)
"cod_obs"



#' Check the structure of a TAM obs list
#'
#' @description
#' Validates the minimal structure required by [make_dat()] and [fit_tam()].
#' The `obs` object must be a list containing data.frames `catch`, `index`,
#' `weight`, and `maturity`. Each table must include columns `year`, `age`,
#' and `obs`. Types must be numeric (integerish allowed) for `year` and `age`,
#' and numeric for `obs`.
#'
#' Additional requirements:
#' - `index$samp_time` is **required**, numeric in **\[0, 1\]**, and **no NA**.
#' - `weight$obs` and `maturity$obs` must not contain `NA`.
#' - `catch`, `weight`, and `maturity` must contain a row for **every**
#'   `(year, age)` combination over the global modeled range
#'   (from min to max year and age across all tables).
#'
#' On failure the function aborts with a `cli` error. On success, it returns
#' `TRUE` (invisibly).
#'
#' @param obs A named list expected to contain data.frames `catch`, `index`,
#'   `weight`, and `maturity` with the columns described above.
#'
#' @return `TRUE` (invisibly) if validation passes; otherwise aborts.
#' @examples
#' data(cod_obs)
#' check_obs(cod_obs)  # returns TRUE (invisibly) if valid
#'
#' @importFrom cli cli_abort cli_warn
#' @export
check_obs <- function(obs) {
  if (!is.list(obs)) {
    cli::cli_abort(c(
      "{.strong Invalid input}",
      "x" = "`obs` must be a list."
    ))
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

  # helper: integerish (allow numeric with no fractional part)
  .is_integerish <- function(x) is.numeric(x) && all(is.finite(x)) && all(x == as.integer(x))

  # per-table checks
  check_one <- function(x, nm) {
    if (!is.data.frame(x)) {
      cli::cli_abort(c(
        "{.strong Invalid obs {nm}}",
        "x" = "`{nm}` must be a data.frame (got {class(x)[1]})."
      ))
    }

    need <- c("year", "age", "obs")
    miss <- setdiff(need, names(x))
    if (length(miss)) {
      cli::cli_abort(c(
        "{.strong Invalid obs {nm}}",
        "x" = "`{nm}` is missing required column(s): {paste(miss, collapse = \", \")}"
      ))
    }

    # type and NA checks
    if (!.is_integerish(x$year)) {
      cli::cli_abort(c("{.strong Invalid obs {nm}}",
                       "x" = "`{nm}$year` must be numeric/integer with no NA and no fractional values."))
    }
    if (!.is_integerish(x$age)) {
      cli::cli_abort(c("{.strong Invalid obs {nm}}",
                       "x" = "`{nm}$age` must be numeric/integer with no NA and no fractional values."))
    }
    if (!is.numeric(x$obs)) {
      cli::cli_abort(c("{.strong Invalid obs {nm}}",
                       "x" = "`{nm}$obs` must be numeric (NA allowed in `catch`/`index`, not in `weight`/`maturity`)."))
    }

    # table-specific NA policy for obs
    if (nm %in% c("weight", "maturity") && any(is.na(x$obs))) {
      cli::cli_abort(c(
        "{.strong Missing values in {nm}$obs}",
        "x" = "`{nm}$obs` must be available for all years and ages (no NA)."
      ))
    }

    # index survey and samp_time: required, numeric in [0,1], no NA
    if (nm == "index") {
      if (!("survey" %in% names(x))) {
        cli::cli_abort(c(
          "{.strong Missing `index$survey`}",
          "x" = "Please supply survey name(s)."
        ))
      }
      if (!("samp_time" %in% names(x))) {
        cli::cli_abort(c(
          "{.strong Missing `index$samp_time`}",
          "x" = "Please supply survey timing as a fractional year in [0, 1] (no NA)."
        ))
      }
      ok_num <- is.numeric(x$samp_time)
      ok_rng <- all(x$samp_time >= 0 & x$samp_time <= 1, na.rm = TRUE)
      ok_na  <- !any(is.na(x$samp_time))
      if (!(ok_num && ok_rng && ok_na)) {
        cli::cli_abort(c(
          "{.strong Invalid `index$samp_time`}",
          "x" = "`index$samp_time` must be numeric in [0, 1] with no NA."
        ))
      }
    }

    invisible(TRUE)
  }

  for (nm in required_tabs) check_one(obs[[nm]], nm)

  # global year/age range (min..max), then require full grid in chosen tables
  obs_years <- unlist(lapply(obs, function(d) d$year))
  obs_ages  <- unlist(lapply(obs, function(d) d$age))
  years <- seq.int(min(obs_years, na.rm = TRUE), max(obs_years, na.rm = TRUE))
  ages  <- seq.int(min(obs_ages,  na.rm = TRUE), max(obs_ages,  na.rm = TRUE))
  grid_size <- length(years) * length(ages)

  must_be_full <- c("catch", "weight", "maturity")
  for (nm in must_be_full) {
    xy <- unique(obs[[nm]][, c("year", "age")])
    if (nrow(xy) != grid_size) {
      cli::cli_abort(c(
        "{.strong `{nm}` coverage is incomplete}",
        "x" = "`{nm}` must contain a row for every (year, age) in the modeled range: {years[1]}-{years[length(years)]} x {ages[1]}-{ages[length(ages)]}.",
        "i" = "Current unique (year, age) rows in `{nm}`: {nrow(xy)} (expected {grid_size})."
      ))
    }
  }

  invisible(TRUE)
}

