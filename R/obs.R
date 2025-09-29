
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
#'   - `obs`: numeric observed catch in numbers (may include `NA`).
#'   - `sd_obs_block`: character; observation-variance block label (`"catch"`).
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
#'   - `samp_time`: numeric sampling time within the year (fraction of a year; e.g., `0.8`).
#'   - `q_block`: factor; age blocks for catchability.
#'   - `sd_obs_block`: character; observation-variance block label (`"index"`).
#'
#' - **weight**: weight-at-age and mortality assumption used for \eqn{M} mean structure.
#'   Columns:
#'   - `year`: integer year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric weight-at-age (same units as source; e.g., kg).
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
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
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

    # index samp_time: required, numeric in [0,1], no NA
    if (nm == "index") {
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



#' Append projection rows across TAM obs tables
#'
#' @description
#' Appends `n_proj` projection years to each present table among `catch`,
#' `index`, `weight`, `maturity`.
#' For each table:
#' - averages `weight$obs` and `maturity$obs` over the last `n_mean` years by `age`;
#' - sets `catch$obs` and `index$obs` to `NA` across projection years;
#' - copies all other columns from the terminal year;
#' - adds `is_proj` column.
#'
#' @param obs A named list containing data.frames `catch`, `index`, `weight`,
#'   and `maturity` with columns `year`, `age`, `obs` (e.g., [cod_obs]).
#' @param n_proj Integer; number of projection years to append (default `3`).
#' @param n_mean Integer; number of terminal years to average for `obs`
#'   (default `3`). Mean `obs` values are replicated across `is_proj` years
#'   for `weight`, and `maturity` tables. Values other than `obs`
#'   are not averaged and `catch$obs` and `index$obs` values are `NA` across `is_proj`
#'   years.
#' @param quiet Logical; if `FALSE` (default) prints brief messages about what is done.
#'
#' @return The same list structure as `obs`, with projection rows (and `is_proj`)
#'   appended to each core table.
#'
#' @examples
#' data(cod_obs)
#' max_year <- max(cod_obs$catch$year)
#' obs2 <- add_proj_rows(cod_obs, n_proj = 3, n_mean = 3)
#' lapply(obs2, subset, year %in% c(max_year, max_year + 1))
#'
#' @seealso [check_obs()]
#' @export
add_proj_rows <- function(obs, n_proj = 3, n_mean = 3, quiet = FALSE) {
  check_obs(obs)

  has_proj <- sapply(names(obs), function(nm) "is_proj" %in% names(obs[[nm]]))
  if (any(has_proj)) {
    tabs <- paste(names(obs)[has_proj], collapse = ", ")
    cli::cli_abort(c(
      "{.strong Projections already present}",
      "x" = "`is_proj` column already exists in: {tabs}."
    ))
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "Adding {n_proj} projection year(s) using {n_mean}-year mean of {.field weight$obs}, and {.field maturity$obs} by age.",
      "i" = "Projected {.field catch$obs}, and {.field index$obs} set to NA.",
      "i" = "Non-core columns are copied from the terminal year",
      "i" = "For advanced assumptions, manually add {.field is_proj} rows."
    ))
  }

  .add_one <- function(x) {
    max_year   <- max(x$year, na.rm = TRUE)
    proj_years <- seq.int(max_year + 1L, max_year + n_proj)
    mean_years <- seq.int(max_year - n_mean + 1L, max_year)
    aux <- x[x$year == max_year, setdiff(names(x), c("year", "obs")), drop = FALSE]
    mean_obs <- aggregate(obs ~ age, mean, dat = x, subset = year %in% mean_years) |>
      merge(aux, by = "age")
    proj_grid <- expand.grid(year = proj_years, age = mean_obs$age)
    proj_rows <- merge(proj_grid, mean_obs, by = "age")[, names(x)]
    proj_rows$is_proj <- TRUE
    x$is_proj <- FALSE
    rbind(x, proj_rows)
  }
  obs_with_proj <- lapply(obs, .add_one)
  obs_with_proj$catch$obs[obs_with_proj$catch$is_proj] <- NA
  obs_with_proj$index$obs[obs_with_proj$index$is_proj] <- NA
  obs_with_proj
}



