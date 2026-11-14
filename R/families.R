#' Observation families for catch and survey data
#'
#' The default lognormal family fits positive observations. The Poisson-link
#' delta lognormal family also uses sampled zeros to inform catch or survey
#' predictions.
#'
#' @param link Link for the lognormal observation prediction. Only `"log"` is supported.
#' @param y50 Optional positive detection scale in observation units. At this
#'   predicted level the probability of observing zero is one half. The default
#'   estimates a separate value for each selected fleet or survey.
#'
#' @details
#' Supply a family in `catch_settings$family` or `index_settings$family`.
#' A single family applies to the whole observation type. A named list selects
#' families by fleet or survey, with the default lognormal for unlisted series.
#' Unlabelled catches use the name `"Fishery"`.
#'
#' With `lognormal()`, the prediction is the conditional median and sampled
#' zeros are excluded. With `pld_lognormal()`, the prediction \eqn{\widetilde Y}
#' is the arithmetic mean over zeros and positives conditional on population
#' states. The probability of zero is
#' \deqn{P(Y=0)=\exp[-\log(2)\widetilde Y/\widetilde Y_{50}].}
#' Writing \eqn{p = P(Y>0)}, the mean of positive observations is
#' \eqn{\widetilde Y/p}. Their lognormal log location is
#' \eqn{\log(\widetilde Y/p)-\sigma^2/2}. The observation SD settings specify
#' \eqn{\sigma}, the log-SD of positive observations. The `q_link` setting
#' continues to determine how the catchability formula maps to catchability.
#'
#' Use zero for a sampled nondetection and `NA` for a missing observation.
#' Estimating the detection scale requires both sampled zeros and positives.
#' A fixed scale can be used with all-positive data. Each fitted series needs
#' at least one positive observation. Catch at structurally unfished ages is
#' deterministic. Supply pooled terminal-age observations for PLD series at
#' the modeled plus age. Biological inputs can extend to older ages.
#' If only pooled catches are available with older biological ages, use
#' `proj_settings = list(F_mult = 1)` in [fit_hindcast()]. Observed-catch
#' hindcasts require catches by biological age to compute held-out yield.
#'
#' PLD assumes independent observation errors conditional on the fitted states
#' and a stable detection scale across years and ages within each series.
#' Sampling effort, the number of fish aged and common sampling units can
#' alter detection or induce dependence. Check how fitted zero probabilities
#' and positive magnitudes agree with the observations. At very low predictions
#' the positive mean approaches \eqn{\widetilde Y_{50}/\log(2)}.
#'
#' @return A family specification for [fit_tam()] or [prepare_tam()].
#' @references
#' Thorson, J. T. (2018). Three problems with the conventional delta-model for
#' biomass sampling data, and a computationally efficient alternative.
#' Canadian Journal of Fisheries and Aquatic Sciences.
#' \doi{10.1139/cjfas-2017-0266}
#' @seealso [tidy_obs_pred()], [residuals.tam_fit()]
#' @examples
#' lognormal()
#' pld_lognormal()
#' list(RV = pld_lognormal(y50 = 100), Sentinel = lognormal())
#' @name observation_families
#' @export
lognormal <- function(link = "log") {
  if (!identical(link, "log")) cli::cli_abort("The lognormal family supports {.code link = 'log'}.")
  structure(list(name = "lognormal", link = link), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
pld_lognormal <- function(y50 = NULL) {
  if (!is.null(y50) && (!is.numeric(y50) || length(y50) != 1L || !is.finite(y50) || y50 <= 0)) {
    cli::cli_abort("{.arg y50} must be NULL or one positive finite number.")
  }
  structure(list(name = "pld_lognormal", y50 = y50), class = "tam_observation_family")
}

.obs_families <- function(x, series, setting) {
  if (is.null(x)) x <- lognormal()
  out <- stats::setNames(rep(list(lognormal()), length(series)), series)
  if (inherits(x, "tam_observation_family")) {
    return(stats::setNames(rep(list(x), length(series)), series))
  }
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))) ||
      !all(vapply(x, inherits, logical(1), "tam_observation_family"))) {
    cli::cli_abort("{setting} must be a family constructor or a list of families named by series.")
  }
  unknown <- setdiff(names(x), series)
  if (length(unknown)) cli::cli_abort("{setting} has unknown series {.val {unknown}}.")
  out[names(x)] <- x
  out
}

.observation_series <- function(d, type) {
  if (type == "index") as.character(d$survey) else
    if (is.null(d$fleet)) rep("Fishery", nrow(d)) else as.character(d$fleet)
}

.obs_has_pld <- function(dat) length(dat$pld$series) > 0L

.prepare_observation_rows <- function(dat, d, families) {
  series <- ifelse(d$type == "catch", .observation_series(d, "catch"),
    .observation_series(d, "index"))
  d$family <- "lognormal"
  d$is_structural <- FALSE
  catch <- d$type == "catch"
  column <- match(paste(series[catch], d$age[catch]), paste(dat$F_columns$fleet, dat$F_columns$age))
  d$is_structural[catch] <- dat$F_settings$state_index[column] == 0L
  d$is_structural[catch] <- d$is_structural[catch] |
    d$year[catch] %in% dat$proj_years[dat$proj_settings$catch == 0]
  d$is_missing <- is.na(d$obs)
  d$is_observed <- !d$is_proj & !d$is_structural & !d$is_missing & d$obs > 0
  d$log_obs <- rep(NA_real_, nrow(d))
  d$log_obs[d$is_observed] <- log(d$obs[d$is_observed])

  pld_series <- list()
  pld_group <- integer(nrow(d))
  for (type in names(families)) {
    for (name in names(families[[type]])) {
      spec <- families[[type]][[name]]
      if (spec$name != "pld_lognormal") next
      rows <- which(d$type == type & series == name)
      sample <- rows[!d$is_missing[rows] & !d$is_proj[rows] & !d$is_structural[rows]]
      positive <- sample[d$obs[sample] > 0]
      zeros <- sample[d$obs[sample] == 0]
      label <- paste(type, name)
      if (!length(positive)) cli::cli_abort("PLD series {.val {label}} needs at least one positive observation in the fitting years.")
      if (is.null(spec$y50) && !length(zeros)) {
        cli::cli_abort(c(
          "PLD series {.val {label}} has no eligible sampled zeros.",
          "i" = "Use {.fn lognormal} or supply a fixed {.arg y50} with {.fn pld_lognormal}."
        ))
      }
      d$family[rows] <- spec$name
      d$is_observed[sample] <- TRUE
      d$log_obs[rows[!d$is_missing[rows] & d$obs[rows] == 0]] <- -Inf
      unit <- if (is.null(d$unit)) NA_character_ else as.character(d$unit[positive[1L]])
      pld_series[[length(pld_series) + 1L]] <- list(type = type, name = name, unit = unit,
        fixed = !is.null(spec$y50), log_start = if (is.null(spec$y50)) {
          log(log(2)) + stats::median(d$log_obs[positive])
        } else log(spec$y50), positive = length(positive), zero = length(zeros))
      pld_group[rows] <- length(pld_series)
    }
  }

  dat$obs_values <- d$obs
  dat$obs_map <- d[, setdiff(names(d), c("obs", "log_obs"))]
  dat$log_obs <- d$log_obs
  dat$is_missing <- is.na(d$log_obs)
  dat$is_observed <- d$is_observed
  # Only positive lognormal observations enter the Gaussian/OSA tape.
  lognormal_rows <- d$is_observed & pld_group == 0L
  dat$observed <- d$log_obs[lognormal_rows]
  if (length(pld_series)) {
    dat$pld <- list(series = pld_series, group = pld_group,
      lognormal = lognormal_rows,
      positive = which(d$is_observed & pld_group > 0L & d$obs > 0),
      zero = which(d$is_observed & pld_group > 0L & d$obs == 0))
  }
  dat
}

.y50_parameters <- function(dat, type) {
  specs <- Filter(function(x) x$type == type && !x$fixed, dat$pld$series)
  stats::setNames(vapply(specs, `[[`, numeric(1), "log_start"),
    vapply(specs, `[[`, character(1), "name"))
}

.y50_starts <- function(value, before, after, type) {
  old <- Filter(function(x) x$type == type && !x$fixed, before$pld$series)
  new <- Filter(function(x) x$type == type && !x$fixed, after$pld$series)
  retained <- vapply(new, function(series) {
    any(vapply(old, function(x) identical(x$name, series$name) &&
      identical(x$unit, series$unit), logical(1)))
  }, logical(1))
  value[intersect(names(value), vapply(new[retained], `[[`, character(1), "name"))]
}

.tidy_y50 <- function(fit, interval) {
  if (!.obs_has_pld(fit$dat)) return(NULL)
  estimate <- .tam_parameter_summary(fit, "Estimate")
  se <- .tam_parameter_summary(fit, "Std. Error")
  z <- stats::qnorm((1 + interval) / 2)
  stack_list(lapply(fit$dat$pld$series, function(s) {
    parameter <- paste0("log_y50_", s$type)
    e <- if (s$fixed) s$log_start else unname(estimate[[parameter]][s$name])
    err <- if (s$fixed) NA_real_ else unname(se[[parameter]][s$name])
    data.frame(par = paste0("y50_", s$type), coef = s$name,
      fleet = if (s$type == "catch") s$name else NA_character_,
      survey = if (s$type == "index") s$name else NA_character_,
      unit = s$unit, is_fixed = s$fixed, est = exp(e), se = err, se_scale = "log",
      lwr = exp(e - z * err), upr = exp(e + z * err))
  }), label = NULL)
}

.y50_advisories <- function(fit) {
  if (!.obs_has_pld(fit$dat)) return(NULL)
  table <- .tidy_y50(fit, 0.95)
  weak <- vapply(fit$dat$pld$series, function(s) min(s$positive, s$zero) < 5L, logical(1)) |
    !is.finite(table$se) | table$se > 1
  table <- table[weak & !table$is_fixed, , drop = FALSE]
  data.frame(issue = rep("y50_support", nrow(table)), detail = if (nrow(table)) {
    paste("Detection scale for", table$coef,
      "has few zeros or positives, or imprecise uncertainty. Inspect y50 and zero calibration.")
  } else character())
}

# Smooth evaluation of the analytic low-intensity limit. Adding exp(-35) keeps
# the positive-mean ratio finite even when the original intensity underflows.
# The same adjustment is removed from log p. The upper transform prevents
# overflow in the positive tail, where log p has already rounded to zero.
.pld_positive <- function(log_pred, log_y50, sd) {
  log_scale <- log_y50 - log(log(2))
  log_intensity <- log_pred - log_scale
  zero <- rep(0, length(log_intensity))
  shifted <- RTMB::logspace_add(log_intensity, zero - 35)
  safe <- exp(-RTMB::logspace_add(-shifted, zero - 36))
  log_p_shifted <- RTMB::logspace_sub(zero, -safe)
  list(
    log_intensity = log_intensity,
    log_p_positive = log_p_shifted - RTMB::logspace_add(zero, -35 - log_intensity),
    log_median = log_scale + (shifted - log_p_shifted) - sd^2 / 2
  )
}

.pld_predictions <- function(dat, par, log_pred, sd, zero_catch) {
  "[<-" <- RTMB::ADoverload("[<-")
  log_p_zero <- rep(-Inf, length(log_pred))
  log_p_positive <- numeric(length(log_pred))
  log_median <- log_pred
  enabled <- dat$pld$group > 0L
  log_p_zero[enabled & zero_catch] <- 0
  log_p_positive[enabled & zero_catch] <- -Inf
  for (j in seq_along(dat$pld$series)) {
    s <- dat$pld$series[[j]]
    log_y50 <- if (s$fixed) s$log_start else par[[paste0("log_y50_", s$type)]][[s$name]]
    rows <- which(dat$pld$group == j & !zero_catch)
    if (!length(rows)) next
    p <- .pld_positive(log_pred[rows], log_y50, sd[rows])
    log_p_zero[rows] <- -exp(p$log_intensity)
    log_p_positive[rows] <- p$log_p_positive
    log_median[rows] <- p$log_median
  }
  list(log_p_zero = log_p_zero, log_p_positive = log_p_positive, log_pred_positive = log_median)
}
