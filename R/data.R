#' Northern cod assessment example data
#'
#' @description
#' Age-structured fishery, survey, biological, and covariate data for Northern
#' cod (*Gadus morhua*) in Northwest Atlantic Fisheries Organization (NAFO)
#' Divisions 2J3KL. The data are adapted from Northern cod assessment inputs
#' and arranged for use with [fit_tam()]. They provide a reproducible example
#' for fitting and comparing alternative observation and population-process
#' models; they are not a reproduction of the full official assessment.
#'
#' @format
#' A `tam_data` object containing five data frames:
#'
#' - **`catch`**: Commercial catch numbers at age (`year`, `age`, `obs`), in
#'   thousands of fish; observations can be missing.
#' - **`index`**: DFO autumn research-vessel (RV) abundance index at age
#'   (`year`, `age`, `obs`), with `survey`, survey timing `samp_time = 0.8`,
#'   and catchability age group `q_block`. Survey observations have been expanded using
#'   12,048,556 trawl units and expressed in thousands of fish.
#' - **`weight`**: Annual stock weight at age (`obs`, kg).
#' - **`maturity`**: Annual proportion mature at age (`obs`, from 0 to 1).
#' - **`pop_covariates`**: A complete year-age table with `F_y_block` and
#'   `F_a_block` (candidate fishing-mortality mean groups), `M_a_block`
#'   (candidate natural-mortality age groups), `allometric_M_assumption`
#'   (a weight-based baseline M, per year).
#'
#' Catch and RV indices have `unit = "thousand fish"`; individual weights have
#' `unit = "kg"`. Multiplying thousands of fish by kg per fish yields tonnes.
#' The additional columns are candidate explanatory variables, not fitted
#' effects unless referenced in the corresponding model formulas.
#'
#' @details
#' **Assessment provenance.** The data originate from inputs to the Northern
#' cod assessment model (NCAM), used in the extended framework (xteNCAM).
#' The stock's 2026 assessment, based on observations available through 2025,
#' is documented by Regular et al. (in press) and summarized in DFO (2026).
#' The framework methods and additional data sources are documented by
#' DFO (2024) and Regular et al. (2025a, 2025b). Only a subset of the assessment
#' inputs is included here: for example, tagging observations, juvenile and
#' Sentinel surveys, and the full Smith Sound likelihood are not supplied as
#' separate observations. Differences between results from tinyAM and the
#' official assessment should therefore not be interpreted as discrepancies in
#' the published stock assessment.
#'
#' **Natural mortality assumptions.** `allometric_M_assumption` uses
#' `M = M0 * weight^(-0.305)`, following the weight-scaling approach of
#' Lorenzen (1996). The multiplier `M0` is calibrated so that the arithmetic
#' mean of these M values over source weight records for ages 5-14 is 0.30 per
#' year, before restricting the example to the catch years. The mean over the
#' packaged subset can therefore differ from 0.30.
#' This is a constructed example baseline, not an independently estimated
#' natural-mortality series.
#'
#' **Mortality blocks.** `M_a_block` defines broad age groups that can be used
#' to allow natural mortality responses to differ among younger and older fish.
#' The fishing-mortality variables `F_y_block` and `F_a_block` provide analogous
#' options for modelling fishing mortality. `F_y_block` divides the time series
#' into periods broadly reflecting changes in the Northern cod fishery, including
#' the 1992 moratorium, subsequent limited fishing and closures, the
#' establishment of the stewardship fishery in 2006, and the return to
#' commercial fishing in 2024. `F_a_block` distinguishes individual ages,
#' allowing age-specific fishing mortality or selectivity to be modelled.
#' These groupings are illustrative model-design variables, not observations
#' or assumptions imposed by the accepted assessment. They have no effect
#' unless explicitly included in the corresponding model settings, and
#' alternative groupings may be defined to explore competing hypotheses.
#'
#' @source
#' Derived from `NCAM::inputs`, with covariate preparation in
#' `data-raw/cod_obs.R`. See the cited assessment and framework publications
#' for the underlying surveys, fishery data, model structure, and provenance.
#'
#' @references
#' DFO. 2024. Northern (2J3KL) Atlantic Cod Assessment Framework.
#' DFO Can. Sci. Advis. Sec. Sci. Advis. Rep. 2024/046.
#' \url{https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-6/Fs70-6-2024-046-eng.pdf}
#'
#' DFO. 2026. NAFO Divisions 2J3KL Northern Cod
#' (*Gadus morhua*) Stock Assessment to 2026.
#' DFO Can. Sci. Advis. Sec. Sci. Advis. Rep. 2026/030.
#' \url{https://csas-scas.dfo-mpo.gc.ca/publications-publications/f7f156a8-b0ef-4840-b2e0-4573fb325faa}
#'
#' Regular, P.M., Robertson, G.J., Kumar, R., Varkey, D.A., Gregory, R.S.,
#' Lewis, R.S., Skanes, K., Gullage, N., Koen-Alonso, M., and Dwyer, K.S.
#' 2025a. Extending the Northern Cod Assessment Model - Part I: Bridging
#' Gaps with Additional Data and Model Variations. DFO Can. Sci. Advis.
#' Sec. Res. Doc. 2025/034.
#' \url{https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-5/Fs70-5-2025-034-eng.pdf}
#'
#' Regular, P.M., Kumar, R., Varkey, D.A., Koen-Alonso, M., and Stenson,
#' G.B. 2025b. Extending the Northern Cod Assessment Model - Part II:
#' Quantifying the Impact of Capelin and Seals. DFO Can. Sci. Advis.
#' Sec. Res. Doc. 2025/035.
#' \url{https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-5/Fs70-5-2025-035-eng.pdf}
#'
#' Regular, P.M., Skanes, K., Novaczek, E., Robertson, G.J., Rideout, R.M.,
#' Trueman, S., Kumar, R., Gregory, R.S., Koen-Alonso, M., and Munro, H.
#' In press. Assessment of the Northern Cod (*Gadus morhua*) Stock in NAFO
#' Divisions 2J3KL in 2026. DFO Can. Sci. Advis. Sec. Res. Doc.
#'
#' Lorenzen, K. 1996. The relationship between body weight and natural
#' mortality in juvenile and adult fish: a comparison of natural
#' ecosystems and aquaculture. J. Fish Biol. 49: 627-642.
#' \doi{10.1111/j.1095-8649.1996.tb00060.x}
#'
#' @examples
#' data(cod_obs)
#' names(cod_obs)
#' print(cod_obs)
#'
#' @seealso [as_tam_data()], [fit_tam()], [fit_hindcast()]
#' @docType data
#' @keywords datasets
#' @name cod_obs
#' @usage data(cod_obs)
"cod_obs"



#' Organize tinyAM input data
#'
#' @description
#' Add the optional `tam_data` class to a list of tinyAM input tables. Printing
#' shows a compact preview of each table; table values and types are unchanged.
#' Plain lists can also be passed directly to [prepare_tam()] and [fit_tam()].
#'
#' @details
#' The required tables are `catch`, `index`, `weight` and `maturity`.
#' An optional `pop_covariates` table supplies columns used by N, F and M
#' formulas, including recruitment covariates and supplied M assumptions.
#' It has one row per input year and age, with `year` and `age` keys and no
#' required `obs` column. Values may be missing where a formula does not use
#' them. Recruitment covariates are read at the youngest modeled age.
#'
#' Population formulas use these columns and intrinsic `year` and `age`, never
#' search other tables. Process-SD designs must remain constant across years
#' and within shared M age blocks. For a model plus group, the supplied row at
#' its terminal age is used; covariates are not averaged or summed. Projections
#' carry terminal covariates forward by age. Modeled SSB is calculated internally.
#'
#' Catch/index observation-SD covariates and survey catchability covariates stay
#' in their respective observation tables. Weight, maturity and spawning timing
#' stay in the biological tables. Validation does not add the class, and inputs
#' are checked again during preparation after any edits.
#'
#' The optional `weight$type` key distinguishes `"stock"` and `"catch"`
#' average individual weights. Stock weights are required and determine biomass
#' and SSB. Catch weights determine yield and catch targets; if absent, stock
#' weights serve both purposes. Each supplied type needs the complete same
#' year-age grid, including retained biological ages beyond a modeled plus age.
#' Missing catch-weight cells never fall back individually to stock weights.
#' For a prepared wide table with `stock` and `catch` columns, the tidy form is:
#' ```r
#' tidyr::pivot_longer(wide, c(stock, catch),
#'   names_to = "type", values_to = "obs")
#' ```
#'
#' Optional `unit` columns label catch numbers, average individual weight and each
#' survey index. Use one unit throughout catch and weight, and one per survey.
#' Values are never converted during fitting. For example, `"thousand fish"`
#' catch and `"kg"` weight give abundance in thousands and biomass/SSB in
#' tonnes. Unknown units are left unlabeled; descriptive index units such as
#' `"fish per tow"` are retained. Reported tables and default plot labels carry
#' these units. Explicit plot labels take precedence. Matrix views convert
#' recognized abundance/mass units to fish/kg for comparisons, independently
#' of aggregate comparison scales; missing reference units remain unavailable.
#'
#' @param data A named list of input data.frames; see [check_tam_data()].
#' @return The same list with class `tam_data` added.
#' @examples
#' data <- as_tam_data(cod_obs)
#' data
#' data$pop_covariates$age_group <- factor(data$pop_covariates$age)
#' check_tam_data(data)
#' @export
as_tam_data <- function(data) {
  check_tam_data(data)
  class(data) <- unique(c("tam_data", class(data)))
  data
}

#' Construct tinyAM input data from tables
#'
#' Assemble prepared input tables and validate them. The returned `tam_data`
#' list prints compact table previews. Plain lists remain valid model inputs;
#' use [check_tam_data()] to validate them or [as_tam_data()] to add the class.
#'
#' @param catch,index,weight,maturity Required input tables; see [check_tam_data()].
#' @param pop_covariates Optional year-by-age population covariates.
#' @return A validated `tam_data` list, with table values and attributes unchanged.
#' @examples
#' dat <- make_tam_data(cod_obs$catch, cod_obs$index,
#'                      cod_obs$weight, cod_obs$maturity,
#'                      cod_obs$pop_covariates)
#' check_tam_data(dat)
#' @export
make_tam_data <- function(catch, index, weight, maturity, pop_covariates = NULL) {
  data <- list(catch = catch, index = index, weight = weight, maturity = maturity)
  if (!is.null(pop_covariates)) data$pop_covariates <- pop_covariates
  as_tam_data(data)
}

#' @rdname as_tam_data
#' @param x A `tam_data` object.
#' @param n Number of first and last rows to show per table (default five).
#' @param ... Arguments passed to table printing.
#' @return The print method returns `x` invisibly.
#' @export
print.tam_data <- function(x, n = 5, ...) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n)) {
    cli::cli_abort("n must be one positive integer.")
  }
  cli::cli_h1("tinyAM data")
  for (nm in names(x)) {
    d <- x[[nm]]
    cli::cli_h2("{nm}: {nrow(d)} rows, {ncol(d)} columns")
    if (nrow(d)) {
      cli::cli_text("Years {min(d$year)}\u2013{max(d$year)}; ages {min(d$age)}\u2013{max(d$age)}")
      rows <- if (nrow(d) <= 2 * n) {
        seq_len(nrow(d))
      } else {
        c(seq_len(n), seq.int(nrow(d) - n + 1L, nrow(d)))
      }
      preview <- d[rows, , drop = FALSE]
      if (nrow(d) > 2 * n) {
        preview <- format(as.data.frame(preview), ...)
        separator <- preview[1L, , drop = FALSE]
        separator[] <- "\u22ee"
        rownames(separator) <- "\u22ee"
        preview <- rbind(
          preview[seq_len(n), , drop = FALSE], separator,
          preview[seq.int(n + 1L, nrow(preview)), , drop = FALSE]
        )
      }
      print(preview, ...)
    }
  }
  invisible(x)
}

# One model-surface frame, ordered year within age, for every population formula.
.population_data <- function(dat, formula = NULL, label = "Population formula") {
  d <- dat$pop_data
  if (!is.null(formula)) {
    if (.contains_power(formula)) cli::cli_abort("power() is only supported in index_settings$q_form.")
    variables <- function(x) {
      if (!is.call(x)) {
        return(x)
      }
      parts <- as.list(x)
      if (.formula_call_name(x) %in% c("iid", "rw", "ar1")) {
        parts <- as.list(match.call(get(.formula_call_name(x), mode = "function"), x))
        parts <- parts[!names(parts) %in% c("sd", "phi")]
      }
      as.call(lapply(parts, variables))
    }
    missing <- setdiff(all.vars(variables(formula[[2L]])), names(d))
    if (length(missing)) {
      cli::cli_abort("{label} needs {paste(missing, collapse = ', ')} in data$pop_covariates.")
    }
  }
  d
}

#' Check tinyAM input data
#'
#' @description
#' Validates the minimal structure required by [prepare_tam()] and [fit_tam()].
#' The `obs` object must be a list containing data.frames `catch`, `index`,
#' `weight`, and `maturity`. Each table must include columns `year`, `age`,
#' and `obs`. Types must be numeric (integerish allowed) for `year` and `age`,
#' and numeric for `obs`.
#'
#' Additional requirements:
#' - `index$samp_time` is **required**, numeric in **\[0, 1\]**, and **no NA**.
#' - `weight$obs` and `maturity$obs` must not contain `NA`.
#' - Optional `maturity$spawn_time` must contain finite numeric fractions in
#'   `[0, 1]` for every year and age. Absence means zero; see [prepare_tam()].
#' - Non-missing observations must be finite and non-negative; maturity must
#'   be a proportion in `[0, 1]`. Zero catch/index values are accepted here but
#'   treated as missing by [prepare_tam()], because the observation model is on logs.
#' - `catch` and `maturity` must contain exactly one row for **every**
#'   `(year, age)` combination over the global modeled range
#'   (from min to max year and age across all tables).
#'   Weight requires the same grid within each supplied type.
#' - Survey tables may omit year-age combinations.
#' - Optional `weight$type` contains `"stock"` or `"catch"`. Stock weights
#'   are required. With types, the unique key is `(year, age, type)` and every
#'   supplied type must cover the complete same grid. Without types, or with
#'   stock alone, stock weights also supply yield and catch targets.
#' - Optional character `unit` columns must be consistent throughout catch and
#'   weight, and within each survey. Missing columns are valid; a supplied column
#'   requires nonempty labels on every row. Units never change input values.
#'   See [as_tam_data()] for units of derived quantities.
#'
#' On failure the function aborts with a `cli` error. On success, it returns
#' `TRUE` (invisibly).
#'
#' @param data A named list expected to contain data.frames `catch`, `index`,
#'   `weight`, and `maturity` with the columns described above.
#'
#' @return `TRUE` (invisibly) if validation passes; otherwise aborts.
#' @examples
#' data(cod_obs)
#' check_obs(cod_obs) # returns TRUE (invisibly) if valid
#'
#' @importFrom cli cli_abort cli_warn
#' @export
check_tam_data <- function(data) {
  obs <- data
  if (!is.list(obs)) {
    cli::cli_abort(c(
      "{.strong Invalid input}",
      "x" = "`obs` must be a list."
    ))
  }

  if (is.null(names(obs)) || anyNA(names(obs)) || any(!nzchar(names(obs))) || anyDuplicated(names(obs))) {
    cli::cli_abort("data must have unique, nonempty table names.")
  }

  required_tabs <- c("catch", "index", "weight", "maturity")
  missing_tabs <- setdiff(required_tabs, names(obs))
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
                       "x" = "`{nm}$year` must be numeric/integer with no NA and no fractional values."
      ))
    }
    if (!.is_integerish(x$age)) {
      cli::cli_abort(c("{.strong Invalid obs {nm}}",
                       "x" = "`{nm}$age` must be numeric/integer with no NA and no fractional values."
      ))
    }
    if (!is.numeric(x$obs)) {
      cli::cli_abort(c("{.strong Invalid obs {nm}}",
                       "x" = "`{nm}$obs` must be numeric (NA allowed in `catch`/`index`, not in `weight`/`maturity`)."
      ))
    }

    # table-specific NA policy for obs
    if (nm %in% c("weight", "maturity") && any(is.na(x$obs))) {
      cli::cli_abort(c(
        "{.strong Missing values in {nm}$obs}",
        "x" = "`{nm}$obs` must be available for all years and ages (no NA)."
      ))
    }
    if (any(!is.finite(x$obs[!is.na(x$obs)])) || any(x$obs < 0, na.rm = TRUE)) {
      cli::cli_abort("{nm}$obs must contain finite, non-negative values (NA allowed for catch and index).")
    }
    if (nm == "maturity" && any(x$obs > 1)) {
      cli::cli_abort("maturity$obs must contain proportions between 0 and 1.")
    }

    if (nm == "maturity" && "spawn_time" %in% names(x)) {
      time <- x$spawn_time
      if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0 | time > 1)) {
        cli::cli_abort("maturity$spawn_time must contain finite fractions between 0 and 1 for every year and age.")
      }
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
      ok_na <- !any(is.na(x$samp_time))
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
  if ("type" %in% names(obs$weight)) {
    type <- as.character(obs$weight$type)
    if (anyNA(type) || any(!type %in% c("stock", "catch"))) {
      cli::cli_abort("weight$type must contain only 'stock' or 'catch', without missing types.")
    }
    if (!"stock" %in% type) cli::cli_abort("Stock weights are required; catch-only weights are not supported.")
  }

  # Covariates cannot change the range supplied by observations and biology.
  obs_years <- unlist(lapply(obs[required_tabs], function(d) d$year))
  obs_ages <- unlist(lapply(obs[required_tabs], function(d) d$age))
  years <- seq.int(min(obs_years, na.rm = TRUE), max(obs_years, na.rm = TRUE))
  ages <- seq.int(min(obs_ages, na.rm = TRUE), max(obs_ages, na.rm = TRUE))
  grid_size <- length(years) * length(ages)

  if (!is.null(obs$pop_covariates)) {
    d <- obs$pop_covariates
    if (!is.data.frame(d) || !all(c("year", "age") %in% names(d))) {
      cli::cli_abort("pop_covariates must be a data.frame with year and age columns.")
    }
    if (!.is_integerish(d$year) || !.is_integerish(d$age)) {
      cli::cli_abort("pop_covariates year and age must be finite numeric integer values.")
    }
    if (any(!d$year %in% years) || any(!d$age %in% ages)) {
      cli::cli_abort("pop_covariates must use the year-age range supplied by catch, index, weight and maturity.")
    }
  }
  must_be_full <- c(
    "catch", "weight", "maturity",
    if (!is.null(obs$pop_covariates)) "pop_covariates"
  )
  for (nm in must_be_full) {
    keys <- c("year", "age", if (nm == "weight" && "type" %in% names(obs$weight)) "type")
    if (anyDuplicated(obs[[nm]][, keys])) {
      cli::cli_abort("{.arg {nm}} must have exactly one row per {paste(keys, collapse = ' and ')}; duplicate rows would misalign the model surfaces.")
    }
    xy <- unique(obs[[nm]][, c("year", "age")])
    sizes <- if (nm == "weight" && "type" %in% names(obs$weight)) {
      vapply(split(obs$weight, obs$weight$type, drop = TRUE), nrow, integer(1))
    } else nrow(xy)
    if (any(sizes != grid_size)) {
      cli::cli_abort(c(
        "{.strong `{nm}` coverage is incomplete}",
        "x" = "`{nm}` must contain a row for every (year, age) in the modeled range: {years[1]}-{years[length(years)]} x {ages[1]}-{ages[length(ages)]}.",
        "i" = if (length(keys) == 3L) "Every weight type must cover the complete stock-weight grid, including retained plus ages." else
          "Current unique (year, age) rows in `{nm}`: {nrow(xy)} (expected {grid_size})."
      ))
    }
  }

  .data_units(obs)
  invisible(TRUE)
}

# Select stock or catch weights; stock is the fallback only when catch is absent.
.select_weights <- function(weight, type = "stock") {
  if (!"type" %in% names(weight)) return(weight)
  if (type == "catch" && !any(weight$type == "catch")) type <- "stock"
  weight[weight$type == type, , drop = FALSE]
}

#' Check the structure of a TAM obs list
#'
#' A convenience wrapper for [check_tam_data()].
#' @param obs A named list of tinyAM input tables.
#' @inherit check_tam_data details
#' @return `TRUE` (invisibly) if validation passes; otherwise aborts.
#' @examples
#' check_tam_data(cod_obs)
#' @export
check_obs <- function(obs) check_tam_data(obs)

