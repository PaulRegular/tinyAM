#' Abundance-dependent catchability
#'
#' Use `power(N)` in `index_settings$q_form` to allow catchability to change
#' with abundance. Ordinary formula terms specify baseline catchability.
#' `N` inside this marker is the model's abundance at survey time,
#' \eqn{A=N\exp(-sZ)}, not an observation column or the observed index.
#'
#' @details With baseline \eqn{\log q_0=X\beta} (including other q effects),
#' \deqn{\log q_{effective}=\log q_0+(p-1)\log A,\qquad
#' \log I\sim N(\log q_0+p\log A,\sigma^2).}
#' The exponent is positive, \eqn{p=\exp(\gamma)}; it starts at one.
#' At one, the index is proportional to abundance and q is constant;
#' above one q increases with abundance,
#' and below one q increases as abundance declines. Error is multiplicative
#' on the observation scale. Abundance is not automatically centered or rescaled;
#' baseline q therefore depends on its units and on the exponent.
#'
#' One additive power term is supported, with the log q link only. It can be
#' combined with ordinary q terms and existing catchability curves and Gaussian
#' effects. Its coefficients are fixed effects, with no additional process SD.
#' Do not place it in interactions or numeric transformations.
#'
#' Without `by`, one exponent is shared across all index observations. With
#' `by`, labels define separate exponents, even when labels are numeric. Missing
#' labels mean no adjustment (p = 1); these observations remain in the likelihood.
#' Use a separate grouping column to apply powers to selected survey ages.
#' Groups need fitted positive observations. Limited abundance variation and
#' flexible q effects can make exponents weakly identified; inspect uncertainty.
#'
#' Projections and simulations use their own current abundance, not terminal
#' effective q. [tidy_obs_pred()] reports effective `q`, `q_base` and `q_power`;
#' [tidy_par()] reports the exponent with log-scale SEs. For fixed-effect
#' coefficient tables, baseline q units depend on p; effective q retains the
#' index-unit/population-unit ratio.
#' Compare `q_base` with a source baseline-q table; baseline q alone is not an
#' effective-q reference when p differs from one.
#'
#' @param x Exactly `N`, the model's abundance adjusted to survey time.
#' @param by Optional index-table grouping column. Missing labels fix p to one.
#' @return A formula marker; direct calls raise an error.
#' @examples
#' ~ q_key + power(N, by = q_key)
#' ~ survey + power(N)
#' # Set q_power_group to NA for survey ages without an estimated power:
#' ~ q_key + power(N, by = q_power_group)
#' @seealso [prepare_tam()], [formula_effects], [catchability_curves]
#' @export
power <- function(x, by = NULL) {
  cli::cli_abort("power() is a formula marker; use power(N) in index_settings$q_form with q_link = 'log'.")
}

.contains_power <- function(x) {
  is.call(x) && (.formula_call_name(x) == "power" ||
    any(vapply(as.list(x)[-1L], .contains_power, logical(1))))
}

.make_q_power <- function(spec, data) {
  group <- if (is.null(spec$by)) rep("all", nrow(data)) else data[[spec$by]]
  if (!is.atomic(group) || !is.null(dim(group)) ||
      (is.numeric(group) && any(!is.finite(group[!is.na(group)])))) {
    cli::cli_abort("power() by must contain categorical labels or missing values.")
  }
  group <- as.character(group)
  if (any(!is.na(group) & !nzchar(trimws(group)))) cli::cli_abort("power() group labels must be nonempty or NA.")
  levels <- unique(group[!is.na(group)])
  index <- match(group, levels, nomatch = 0L)
  active <- is.finite(data$obs) & data$obs > 0
  if ("is_proj" %in% names(data)) active <- active & !data$is_proj
  absent <- setdiff(levels, group[active])
  if (length(absent)) cli::cli_abort("power() groups need fitted positive observations: {paste(absent, collapse = ', ')}.")
  list(by = spec$by, levels = levels, index = index)
}

.q_power_values <- function(par, dat) {
  "[<-" <- RTMB::ADoverload("[<-")
  index <- dat$q_power$index
  values <- rep(1, length(index))
  active <- index > 0L
  if (any(active)) values[active] <- exp(par$log_q_power[index[active]])
  values
}

.q_power_starts <- function(values, before, after) {
  if (is.null(after$q_power)) return(NULL)
  old <- before$obs$index; new <- after$obs$index
  key <- function(d) paste(d$survey, d$year, d$age, sep = ":")
  rows <- match(key(new), key(old))
  shared <- !is.na(rows)
  # Zero indices need an explicit sentinel: indexing a vector with zero drops rows.
  previous <- c(NA_character_, before$q_power$levels)[before$q_power$index[rows[shared]] + 1L]
  current <- c(NA_character_, after$q_power$levels)[after$q_power$index[shared] + 1L]
  retained <- names(values)[vapply(names(values), function(g) {
    any(current == g, na.rm = TRUE) &&
      identical(!is.na(previous) & previous == g, !is.na(current) & current == g)
  }, logical(1))]
  values[retained]
}

.q_power_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  spec <- fit$dat$q_power
  if (is.null(spec)) return(out)
  d <- fit$dat$obs$index
  active <- is.finite(d$obs) & d$obs > 0 & !d$is_proj
  ii <- fit$dat$obs_map$type == "index"
  abundance <- fit$rep$log_pred[ii] - fit$rep$log_q_obs
  for (i in seq_along(spec$levels)) {
    x <- abundance[active & spec$index == i]
    if (length(x) && all(is.finite(x)) && diff(range(x)) < log(1.1)) {
      out <- rbind(out, data.frame(issue = "q_power_support", detail = paste0(
        "Sampled abundance spans less than 10% in power group ", spec$levels[i],
        "; its exponent may be difficult to distinguish from baseline catchability.")))
    }
  }
  out
}
