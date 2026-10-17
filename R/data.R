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
#' Optional `unit` columns label catch numbers, stock weight per fish and each
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
