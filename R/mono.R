#' Monotonic survey catchability in a formula
#'
#' Use `mono()` as an additive term in `index_settings$q_form`. It is a
#' formula marker, not a numeric transformation. Numeric values are ordered
#' increasingly; factors (including ordered factors) use their declared levels.
#' Supply factor levels in the scientifically intended order.
#'
#' The first represented level has no increment. Subsequent log-q levels add
#' cumulative positive magnitudes
#' `dq = exp(log_dq)`. These are increments on the log-q scale, not absolute q.
#' Steps start at 0.05, giving a nearly flat initial curve. Separate levels have
#' strictly positive step magnitudes, which can approach zero; exact plateaus
#' are represented by pooling observations into the same level.
#'
#' `by` gives each group independent steps, using only levels represented in
#' that group, in their declared order. Each group needs at least two levels.
#' Baselines come from ordinary formula terms: `~ mono(x, by = survey)` shares
#' an intercept, whereas `~ survey + mono(x, by = survey)` gives separate
#' survey baselines. Other ordinary covariates can modify q, so monotonicity
#' holds with those covariates held constant. Interactions involving `mono()`,
#' transformed arguments, and ordinary effects of the same `x` are unsupported.
#'
#' @param x Name of a numeric or factor column in the index observations.
#' @param by Optional name of a categorical grouping column.
#' @return A formula marker; calling this function directly raises an error.
#' @seealso [make_dat()], [make_par()], [tidy_obs_pred()], [tidy_par()]
#' @examples
#' ~ q_block # ordinary unconstrained q
#' ~ mono(q_block)
#' ~ survey + mono(q_block, by = survey)
#' @export
mono <- function(x, by = NULL) {
  cli::cli_abort("mono() is only supported as an additive term in {.arg index_settings$q_form}.")
}

.parse_q_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg q_form} must be a formula.")
  }
  contains_mono <- function(x) {
    if (!is.call(x)) return(FALSE)
    if (identical(x[[1L]], as.name("mono"))) return(TRUE)
    if (as.character(x[[1L]])[1L] %in% c("::", ":::") &&
        identical(x[[3L]], as.name("mono"))) return(TRUE)
    any(vapply(as.list(x), contains_mono, logical(1)))
  }
  # Preserve the original model.matrix path, including contrasts and attributes.
  if (!contains_mono(formula)) {
    return(list(q_modmat = stats::model.matrix(formula, data = data)))
  }
  if (length(formula) != 2L) cli::cli_abort("Monotonic {.arg q_form} must be one-sided.")
  specs <- list()
  strip <- function(x) {
    if (!contains_mono(x)) return(x)
    if (is.call(x) && identical(x[[1L]], as.name("mono"))) {
      spec <- tryCatch(match.call(mono, x), error = function(e) {
        cli::cli_abort("Invalid mono() arguments: {conditionMessage(e)}")
      })
      spec <- as.list(spec)[-1L]
      column <- function(arg, label) {
        if (!is.symbol(arg) || !as.character(arg) %in% names(data)) {
          cli::cli_abort("mono() {.arg {label}} must name an existing index column.")
        }
        as.character(arg)
      }
      variable <- column(spec$x, "x")
      by <- if (is.null(spec$by)) NULL else column(spec$by, "by")
      specs[[length(specs) + 1L]] <<- list(variable = variable, by = by)
      return(NULL)
    }
    if (is.call(x) && identical(x[[1L]], as.name("("))) return(strip(x[[2L]]))
    if (is.call(x) && identical(x[[1L]], as.name("-")) && length(x) == 3L &&
        !contains_mono(x[[3L]])) {
      left <- strip(x[[2L]])
      return(call("-", if (is.null(left)) 1 else left, x[[3L]]))
    }
    if (is.call(x) && identical(x[[1L]], as.name("+"))) {
      args <- lapply(as.list(x)[-1L], strip)
      args <- Filter(Negate(is.null), args)
      if (!length(args)) return(NULL)
      if (length(args) == 1L) return(args[[1L]])
      return(as.call(c(list(as.name("+")), args)))
    }
    cli::cli_abort("mono() must be an additive term, without interactions or transformations.")
  }
  ordinary <- formula
  rhs <- strip(formula[[2L]])
  ordinary[[2L]] <- if (is.null(rhs)) 1 else rhs
  variables <- vapply(specs, `[[`, character(1), "variable")
  ordinary_vars <- all.vars(stats::terms(ordinary, data = data))
  if (anyDuplicated(variables) || any(variables %in% ordinary_vars) ||
      any(variables %in% unlist(lapply(specs, `[[`, "by")))) {
    cli::cli_abort("A mono() variable cannot also have ordinary, duplicate, or grouping effects in {.arg q_form}.")
  }
  design <- lapply(specs, .make_mono_q_design, data = data)
  q_modmat <- stats::model.matrix(ordinary, data = data)
  if (nrow(q_modmat) != nrow(data)) {
    cli::cli_abort("Missing ordinary q covariates would misalign mono() observation rows.")
  }
  list(q_modmat = q_modmat,
       q_mono_modmat = do.call(cbind, lapply(design, `[[`, "matrix")),
       q_mono_steps = do.call(rbind, lapply(design, `[[`, "steps")))
}

.make_mono_q_design <- function(spec, data) {
  x <- data[[spec$variable]]
  if (!(is.numeric(x) || is.factor(x)) || !is.null(dim(x)) || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x)))) {
    cli::cli_abort("mono() {.arg x} must be numeric or a factor, with no missing or non-finite values.")
  }
  ordered_levels <- if (is.factor(x)) levels(x) else sort(unique(x))
  group <- if (is.null(spec$by)) rep("all", length(x)) else data[[spec$by]]
  if (!is.atomic(group) || !is.null(dim(group)) || anyNA(group)) {
    cli::cli_abort("mono() {.arg by} must be a categorical column with no missing values.")
  }
  groups <- if (is.factor(group)) levels(droplevels(group)) else unique(group)
  matrices <- steps <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    rows <- which(group == groups[i])
    lev <- ordered_levels[ordered_levels %in% x[rows]]
    if (length(lev) < 2L) {
      cli::cli_abort("mono() needs at least two observed levels in each group ({groups[i]}).")
    }
    k <- length(lev) - 1L
    mat <- matrix(0, nrow(data), k)
    mat[rows, ] <- outer(match(x[rows], lev), seq_len(k), `>`)
    label <- paste0(spec$variable,
                    if (!is.null(spec$by)) paste0("[", spec$by, "=", groups[i], "]"),
                    ":", utils::head(lev, -1L), "->", utils::tail(lev, -1L))
    colnames(mat) <- label
    matrices[[i]] <- mat
    steps[[i]] <- data.frame(coef = label, variable = spec$variable,
      from_level = as.character(utils::head(lev, -1L)), to_level = as.character(utils::tail(lev, -1L)),
      by = if (is.null(spec$by)) NA_character_ else spec$by,
      by_level = if (is.null(spec$by)) NA_character_ else as.character(groups[i]))
  }
  list(matrix = do.call(cbind, matrices), steps = do.call(rbind, steps))
}
