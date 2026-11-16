#' Observation families for catch and survey data
#'
#' Lognormal, Gamma and generalized gamma families fit positive observations.
#' Poisson-link delta (PLD) families also use sampled zeros.
#'
#' @param link Link for the assessment prediction. Only `"log"` is supported.
#'   For PLD families, the link applies to the overall assessment prediction,
#'   the arithmetic mean across zeros and positives. The prediction on the
#'   link scale is \eqn{\log\widetilde{Y}_i}. The complementary log-log link
#'   applies to the encounter probability described in Details.
#' @param y50 Optional positive soft detection scale in the observation's units
#'   (the predicted level at which half of observations are expected to be zero).
#'   By default, the model estimates a separate scale for each selected fleet
#'   or survey.
#' @param Q Optional fixed Prentice shape parameter for generalized gamma.
#'   The default estimates a separate signed shape for each selected series.
#'
#' @details
#' Supply a family in `catch_settings$family` or `index_settings$family`.
#' A single family applies to the whole observation type. A named list selects
#' families by fleet or survey, with the default lognormal for unlisted series.
#' Unlabelled catches use the name `"Fishery"`.
#'
#' Use [stats::Gamma()] with `link = "log"` for Gamma observations. For that
#' family and `gengamma()`, the prediction is the arithmetic mean.
#' Positive-only families exclude sampled zeros.
#'
#' ### Poisson-link delta
#'
#' Let \eqn{Y_i} be a recorded catch or survey index for an age and year, and
#' let \eqn{\widetilde{Y}_i} be its prediction from the population model. The
#' catch prediction comes from the Baranov catch equation. The survey prediction
#' is abundance at sampling time scaled by catchability. Poisson-link delta
#' (PLD) uses this prediction to explain both the chance of recording a positive
#' value and the size of positive observations (Thorson 2018, see also the
#' [sdmTMB Poisson-link delta guide](https://sdmtmb.github.io/sdmTMB/articles/poisson-link.html)).
#'
#' tinyAM adapts Thorson's parameterization by using the assessment prediction
#' as the overall mean. Each fleet or survey \eqn{j} has a soft detection scale
#' \eqn{\widetilde{Y}_{50,j}} (the predicted catch or index at which half of
#' observations are expected to be zero). The prediction and the scale use the
#' observation's units.
#'
#' **Encounter.** Conditional on the population states and parameters, recording
#' a positive value is a binomial outcome with one trial and probability
#' \eqn{p_{+,i}}. The probability of recording zero is the chance of zero
#' events under a Poisson distribution with mean \eqn{\lambda_i}.
#' \deqn{\lambda_i = \log(2)\frac{\widetilde{Y}_i}{\widetilde{Y}_{50,j}},
#' \qquad p_{0,i} = \Pr(Y_i = 0) = e^{-\lambda_i},
#' \qquad p_{+,i} = \Pr(Y_i > 0) = 1-p_{0,i}.}{
#' lambda_i = log(2) * Ytilde_i / Ytilde_50,j,
#' p_0,i = Pr(Y_i = 0) = exp(-lambda_i),
#' p_+,i = Pr(Y_i > 0) = 1 - p_0,i.}
#' Doubling the prediction from the soft detection scale to twice that level
#' lowers the chance of recording zero from 50% to 25%. The encounter component
#' uses a complementary log-log link for the probability of a positive.
#' \deqn{\operatorname{cloglog}(p_{+,i}) =
#' \log\!\big[-\log(1-p_{+,i})\big] = \log\lambda_i.}{
#' cloglog(p_+,i) = log(-log(1 - p_+,i)) = log(lambda_i).}
#'
#' **Positive observations.** The positive component uses a log link for its
#' arithmetic mean \eqn{m_{+,i}}, with
#' \eqn{\log m_{+,i} = \log\widetilde{Y}_i - \log p_{+,i}}. Dividing by the
#' probability of a positive keeps the average across zeros and positives equal
#' to the assessment prediction.
#' \deqn{m_{+,i} = \frac{\widetilde{Y}_i}{p_{+,i}},
#' \qquad p_{+,i}\,m_{+,i} = \widetilde{Y}_i.}{
#' m_+,i = Ytilde_i / p_+,i, and p_+,i * m_+,i = Ytilde_i.}
#'
#' The positive distribution can be lognormal, Gamma, or generalized gamma.
#' For example, with a lognormal positive component,
#' \deqn{\log Y_i \mid Y_i > 0 \sim
#' N\!\left(\log m_{+,i}-\frac{\sigma_i^2}{2},\ \sigma_i^2\right).}{
#' log(Y_i) given Y_i > 0 follows Normal(log(m_+,i) - sigma_i^2 / 2, sigma_i^2).}
#' Here \eqn{\sigma_i} controls the spread on the log scale. Subtracting
#' \eqn{\sigma_i^2/2} gives the positive observations an arithmetic mean of
#' \eqn{m_{+,i}}. The default, positive-only `lognormal()` family uses the
#' assessment prediction as its conditional median.
#'
#' **Likelihood.** A recorded zero contributes the log probability of recording
#' zero. A positive contributes the log probability of recording a positive
#' plus the log density of its size. Let \eqn{g_{+,i}} be the chosen positive
#' density on the log-observation scale.
#' \deqn{\ell_i = \begin{cases}
#' \log p_{0,i}, & Y_i = 0, \\[4pt]
#' \log p_{+,i} + \log g_{+,i}(\log Y_i), & Y_i > 0.
#' \end{cases}}{
#' ell_i = log(p_0,i) if Y_i = 0, or
#' log(p_+,i) + log(g_+,i(log(Y_i))) if Y_i > 0.}
#'
#' The `q_link` setting determines how the catchability formula maps to
#' catchability in the survey prediction.
#'
#' ### Dispersion, shape and observation support
#'
#' For lognormal families, the observation SD settings specify the log-SD of
#' positives. For Gamma families, they specify the coefficient of variation of
#' positives. Gamma shape is the inverse squared coefficient of variation.
#' For generalized gamma, these settings specify the Prentice scale
#' \eqn{\sigma}. Use one constant scale per generalized gamma series. If a
#' fixed \eqn{Q} is negative, supply its scale, with no fitted multiplier,
#' satisfying \eqn{1+\sigma Q>0}. Estimated shapes satisfy this finite-mean
#' condition by construction.
#'
#' Generalized gamma becomes Gamma at \eqn{Q=\sigma} and mean-parameterized
#' lognormal at \eqn{Q=0}. The latter is the same distribution family as
#' `lognormal()`, with the population prediction representing the mean rather
#' than the median. `pld_gengamma(Q = 0)` gives `pld_lognormal()`.
#' Gamma, generalized gamma and PLD observations are conditionally independent.
#' Use separate SD coefficients for series whose dispersion scales differ.
#' Within-year age correlation (`cor_struct = "ar1"`) can be used for ordinary
#' lognormal series, including lognormal series in a mixed-family fit.
#' Standardized and quantile residuals use marginal observation distributions
#' and retain that correlation. Use a Gaussian one-step-ahead algorithm for
#' correlated observations in an all-lognormal fit. See [residuals.tam_fit()].
#'
#' Use zero for a sampled nondetection and `NA` for a missing observation.
#' Estimating the soft detection scale requires both sampled zeros and positives.
#' A fixed scale can be used with all-positive data. Each fitted series needs
#' at least one positive observation. Catch at structurally unfished ages is
#' deterministic.
#'
#' For a native-age series, [prepare_tam()] sums the observations from the
#' modeled plus age through the oldest biological age, separately by year and
#' series. Each age needs one sampled value, including any sampled zeros. An
#' absent age or `NA` makes that terminal observation missing. The chosen
#' positive-component family applies to the pooled value.
#'
#' A series with no recorded values above the modeled plus age in the supplied
#' data is treated as already pooled at that age. Older catch rows may be
#' present with missing observations to complete the input grid. The input
#' resolution is determined across all supplied years and is used in each
#' retrospective fold. Biological inputs can extend to older ages.
#'
#' Observed-catch hindcasts use the original counts and weights by biological
#' age to compute held-out yield. If only pooled counts are available with older
#' biological ages, use `proj_settings = list(F_mult = 1)` in [fit_hindcast()].
#'
#' PLD assumes independent observation errors conditional on the fitted states
#' and a stable soft detection scale across years and ages within each series.
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
#'
#' Prentice, R. L. (1974). A log gamma model and its maximum likelihood
#' estimation. Biometrika. \doi{10.1093/biomet/61.3.539}
#'
#' Dunic, J. C., Conner, J., Anderson, S. C. and Thorson, J. T. (2025).
#' The generalized gamma is a flexible distribution that outperforms
#' alternatives when modelling catch rate data. ICES Journal of Marine Science.
#' \doi{10.1093/icesjms/fsaf040}
#'
#' Albertsen, C. M., Nielsen, A. and Thygesen, U. H. (2017).
#' Choosing the observational likelihood in state-space stock assessment
#' models. Canadian Journal of Fisheries and Aquatic Sciences.
#' \doi{10.1139/cjfas-2015-0532}
#'
#' Monnahan et al. (2025). Testing the validity of lognormal likelihoods for
#' abundance indices in stock assessment models using the generalized gamma
#' distribution. ICES Journal of Marine Science. \doi{10.1093/icesjms/fsaf071}
#' @seealso [tidy_obs_pred()], [residuals.tam_fit()]
#' @examples
#' lognormal()
#' stats::Gamma(link = "log")
#' gengamma()
#' pld_lognormal()
#' pld_gamma(y50 = 100)
#' pld_gengamma(Q = 0)
#' list(RV = pld_lognormal(y50 = 100), Sentinel = lognormal())
#' @name observation_families
#' @export
lognormal <- function(link = "log") {
  if (!identical(link, "log")) cli::cli_abort("The lognormal family supports {.code link = 'log'}.")
  structure(list(name = "lognormal", link = link), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
pld_lognormal <- function(y50 = NULL, link = "log") {
  .validate_pld_link(link)
  .validate_y50(y50)
  structure(list(name = "pld_lognormal", y50 = y50, link = link), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
pld_gamma <- function(y50 = NULL, link = "log") {
  .validate_pld_link(link)
  .validate_y50(y50)
  structure(list(name = "pld_gamma", y50 = y50, link = link), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
gengamma <- function(link = "log", Q = NULL) {
  if (!identical(link, "log")) cli::cli_abort("The generalized gamma family supports {.code link = 'log'}.")
  .validate_gengamma_Q(Q)
  structure(list(name = "gengamma", link = link, Q = Q), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
pld_gengamma <- function(y50 = NULL, Q = NULL, link = "log") {
  .validate_pld_link(link)
  .validate_y50(y50)
  .validate_gengamma_Q(Q)
  structure(list(name = "pld_gengamma", y50 = y50, Q = Q, link = link), class = "tam_observation_family")
}

.validate_pld_link <- function(link) {
  if (!identical(link, "log")) {
    cli::cli_abort("Poisson-link delta families support {.code link = 'log'} for the overall assessment prediction.")
  }
}

.validate_y50 <- function(y50) {
  if (!is.null(y50) && (!is.numeric(y50) || length(y50) != 1L || !is.finite(y50) || y50 <= 0)) {
    cli::cli_abort("{.arg y50} must be NULL or one positive finite number.")
  }
}

.validate_gengamma_Q <- function(Q) {
  if (!is.null(Q) && (!is.numeric(Q) || length(Q) != 1L || !is.finite(Q))) {
    cli::cli_abort("{.arg Q} must be NULL or one finite number.")
  }
}

.is_pld_family <- function(family) family %in% c("pld_lognormal", "pld_gamma", "pld_gengamma")

.is_gengamma_family <- function(family) family %in% c("gengamma", "pld_gengamma")

.positive_family <- function(family) sub("^pld_", "", family)

.is_family_spec <- function(x) {
  inherits(x, "tam_observation_family") || inherits(x, "family")
}

.normalize_family <- function(x, setting) {
  if (inherits(x, "tam_observation_family")) return(x)
  if (!identical(x$family, "Gamma") || !identical(x$link, "log")) {
    cli::cli_abort("{setting} supports {.code stats::Gamma(link = 'log')} and tinyAM family constructors.")
  }
  structure(list(name = "gamma", link = "log"), class = "tam_observation_family")
}

.obs_families <- function(x, series, setting) {
  if (is.null(x)) x <- lognormal()
  out <- stats::setNames(rep(list(lognormal()), length(series)), series)
  if (.is_family_spec(x)) {
    return(stats::setNames(rep(list(.normalize_family(x, setting)), length(series)), series))
  }
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))) ||
      !all(vapply(x, .is_family_spec, logical(1)))) {
    cli::cli_abort("{setting} must be a family or a list of families named by series.")
  }
  unknown <- setdiff(names(x), series)
  if (length(unknown)) cli::cli_abort("{setting} has unknown series {.val {unknown}}.")
  out[names(x)] <- lapply(x, .normalize_family, setting = setting)
  out
}

.observation_series <- function(d, type) {
  if (type == "index") as.character(d$survey) else
    if (is.null(d$fleet)) rep("Fishery", nrow(d)) else as.character(d$fleet)
}

.obs_has_pld <- function(dat) length(dat$pld$series) > 0L

.obs_has_non_lognormal <- function(dat) {
  .obs_has_pld(dat) || any(dat$obs_map$family != "lognormal", na.rm = TRUE)
}

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

  pld_series <- gg_series <- list()
  pld_group <- gg_group <- integer(nrow(d))
  for (type in names(families)) {
    for (name in names(families[[type]])) {
      spec <- families[[type]][[name]]
      if (spec$name == "lognormal") next
      rows <- which(d$type == type & series == name)
      sample <- rows[!d$is_missing[rows] & !d$is_proj[rows] & !d$is_structural[rows]]
      positive <- sample[d$obs[sample] > 0]
      label <- paste(type, name)
      pld <- .is_pld_family(spec$name)
      if (!length(positive)) {
        family <- if (pld) "PLD" else spec$name
        cli::cli_abort("{family} series {.val {label}} needs at least one positive observation in the fitting years.")
      }
      d$family[rows] <- spec$name
      structural_zeros <- rows[d$is_structural[rows] & !d$is_missing[rows] & d$obs[rows] == 0]
      d$log_obs[structural_zeros] <- -Inf
      if (.is_gengamma_family(spec$name)) {
        gg_series[[length(gg_series) + 1L]] <- list(type = type, name = name,
          fixed = !is.null(spec$Q), Q = spec$Q)
        gg_group[rows] <- length(gg_series)
      }
      if (!pld) next
      zeros <- sample[d$obs[sample] == 0]
      if (is.null(spec$y50) && !length(zeros)) {
        cli::cli_abort(c(
          "PLD series {.val {label}} has no eligible sampled zeros.",
          "i" = "Use a positive-only family or supply a fixed {.arg y50} with {.fn {spec$name}}."
        ))
      }
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
  lognormal_rows <- d$is_observed & d$family == "lognormal"
  dat$observed <- d$log_obs[lognormal_rows]
  if (length(gg_series)) dat$gengamma <- list(series = gg_series, group = gg_group)
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
    paste("Soft detection scale for", table$coef,
      "has few zeros or positives, or imprecise uncertainty. Inspect y50 and zero calibration.")
  } else character())
}


.gengamma_parameters <- function(dat, type) {
  specs <- Filter(function(x) x$type == type && !x$fixed, dat$gengamma$series)
  stats::setNames(numeric(length(specs)), vapply(specs, `[[`, character(1), "name"))
}

.gengamma_starts <- function(value, before, after, type) {
  old <- names(.gengamma_parameters(before, type))
  new <- names(.gengamma_parameters(after, type))
  value[intersect(names(value), intersect(old, new))]
}

.gengamma_shape <- function(dat, par, sd) {
  "[<-" <- RTMB::ADoverload("[<-")
  Q <- numeric(length(sd))
  for (j in seq_along(dat$gengamma$series)) {
    s <- dat$gengamma$series[[j]]
    rows <- which(dat$gengamma$group == j)
    if (!length(rows)) next
    if (s$fixed) {
      Q[rows] <- s$Q
    } else {
      eta <- par[[paste0("eta_gengamma_Q_", s$type)]][[s$name]]
      # Softplus grows linearly for large positive eta. The margin keeps the
      # finite-mean boundary strict when softplus rounds to zero.
      Q[rows] <- (1 - 16 * .Machine$double.eps) *
        (RTMB::logspace_add(0, eta) / log(2) - 1) / sd[rows[1L]]
    }
  }
  Q
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
  log_mean <- log_scale + (shifted - log_p_shifted)
  list(
    log_intensity = log_intensity,
    log_p_positive = log_p_shifted - RTMB::logspace_add(zero, -35 - log_intensity),
    log_mean = log_mean,
    log_median = log_mean - sd^2 / 2
  )
}

.pld_predictions <- function(dat, par, log_pred, sd, zero_catch) {
  "[<-" <- RTMB::ADoverload("[<-")
  log_p_zero <- rep(-Inf, length(log_pred))
  log_p_positive <- numeric(length(log_pred))
  log_median <- log_mean <- log_pred
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
    log_mean[rows] <- p$log_mean
    log_median[rows] <- p$log_median
  }
  list(log_p_zero = log_p_zero, log_p_positive = log_p_positive, log_pred_positive = log_median,
    log_mean_positive = log_mean)
}

.gamma_log_density <- function(log_y, log_mean, sd) {
  RTMB::dlgamma(log_y - log_mean, shape = 1 / sd^2, scale = sd^2, log = TRUE)
}

# The log-Gamma normalizing constant after its leading Stirling terms are
# removed. Use inverse shape to include the infinite-shape limit at Q = 0.
.gengamma_stirling <- function(inv_shape) {
  "if" <- RTMB::ADoverload("if")
  inv_shape_safe <- if (inv_shape < .1) inv_shape * 0 + .1 else inv_shape
  shape <- 1 / inv_shape_safe
  if (inv_shape < .1) {
    x <- inv_shape^2
    inv_shape * (1/12 + x * (-1/360 + x * (1/1260 + x * (-1/1680 +
      x * (1/1188 - 691/360360 * x)))))
  } else {
    lgamma(shape) - ((shape - .5) * log(shape) - shape + .5 * log(2 * pi))
  }
}

# ((1+x) log(1+x) - x) / x^2, with its signed series at zero.
.gengamma_mean_ratio <- function(x) {
  "if" <- RTMB::ADoverload("if")
  safe_x <- if (abs(x) < .05) x * 0 + .05 else x
  if (abs(x) < .05) {
    value <- 1 / (13 * 14)
    for (i in 12:1) value <- 1 / (i * (i + 1)) - x * value
    value
  } else {
    ((1 + safe_x) * log1p(safe_x) - safe_x) / safe_x^2
  }
}

# log E[exp(sigma * log(G/k)/Q)] with G ~ Gamma(k, 1), k = Q^-2.
# Taking the difference after Stirling's leading terms cancel avoids
# subtracting two large log Gamma values when Q is close to zero.
.gengamma_mean_adjustment <- function(sd, Q) {
  t <- sd * Q
  inv_shape <- Q^2
  sd^2 * .gengamma_mean_ratio(t) - .5 * log1p(t) +
    .gengamma_stirling(inv_shape / (1 + t)) - .gengamma_stirling(inv_shape)
}

# Evaluate (exp(z) - 1 - z) / z^2 without subtracting nearly equal
# quantities near zero. The Taylor series and regular expression agree to
# numerical precision in the range where both receive appreciable weight.
.gengamma_exponential_ratio <- function(z) {
  value <- 1 / factorial(18)
  for (i in 17:2) value <- 1 / factorial(i) + z * value
  z2 <- z^2
  z18 <- z2^9
  denominator <- 1 + z2 * z18
  value / denominator + (z18 / denominator) * (expm1(z) - z)
}

.gengamma_log_density <- function(log_y, log_mean, sd, Q) {
  if (!inherits(Q, "advector") && Q == 0) {
    return(RTMB::dnorm(log_y, log_mean - sd^2 / 2, sd = sd, log = TRUE))
  }
  location <- log_mean - .gengamma_mean_adjustment(sd, Q)
  w <- (log_y - location) / sd
  -log(sd) - .5 * log(2 * pi) - .gengamma_stirling(Q^2) -
    w^2 * .gengamma_exponential_ratio(Q * w)
}

# Numeric helpers for observation summaries, residuals and simulation.
# They take the log arithmetic mean of the positive component.
.gengamma_log_location <- function(log_mean, sd, Q) {
  if (length(log_mean) && isTRUE(all(sd == sd[1L] & Q == Q[1L]))) {
    return(log_mean - .gengamma_mean_adjustment(sd[1L], Q[1L]))
  }
  log_mean - vapply(seq_along(log_mean), function(i) {
    if (is.na(sd[i]) || is.na(Q[i])) return(NA_real_)
    .gengamma_mean_adjustment(sd[i], Q[i])
  }, numeric(1))
}

.log_gamma_relative_median <- function(shape) {
  median <- stats::qgamma(.5, shape = shape)
  out <- log(median / shape)
  tiny <- which(median == 0 & is.finite(shape) & shape > 0)
  out[tiny] <- (log(.5) + lgamma(shape[tiny] + 1)) / shape[tiny] - log(shape[tiny])
  out
}

.digamma_minus_log_shape <- function(shape) {
  out <- digamma(shape) - log(shape)
  large <- which(shape > 1e4)
  inv <- 1 / shape[large]
  out[large] <- inv * (-.5 + inv * (-1/12 + inv^2 / 120))
  out
}

.positive_log_summary <- function(family, log_mean, sd, Q = 0) {
  n <- length(log_mean)
  sd <- rep_len(sd, n)
  family <- .positive_family(family)
  if (family == "lognormal") {
    location <- log_mean - sd^2 / 2
    return(list(mean_log = location, sd_log = sd, median_log = location))
  }
  if (family == "gamma") {
    shape <- 1 / sd^2
    return(list(mean_log = log_mean + .digamma_minus_log_shape(shape),
      sd_log = sqrt(trigamma(shape)), median_log = log_mean + .log_gamma_relative_median(shape)))
  }
  Q <- rep_len(Q, n)
  location <- .gengamma_log_location(log_mean, sd, Q)
  near <- abs(Q) < 1e-7
  mean_log <- location - sd * Q / 2
  sd_log <- sd * sqrt(1 + Q^2 / 2)
  median_log <- location - sd * Q / 3
  i <- which(!near)
  shape <- Q[i]^-2
  mean_log[i] <- location[i] + sd[i] * .digamma_minus_log_shape(shape) / Q[i]
  sd_log[i] <- sd[i] * sqrt(shape * trigamma(shape))
  median_log[i] <- location[i] + sd[i] * .log_gamma_relative_median(shape) / Q[i]
  list(mean_log = mean_log, sd_log = sd_log, median_log = median_log)
}

.log_gamma_cdf <- function(log_x, shape, lower.tail) {
  out <- stats::pgamma(exp(log_x), shape = shape, lower.tail = lower.tail, log.p = TRUE)
  tiny <- which(is.finite(log_x) & log_x < log(.Machine$double.xmin))
  lower <- shape[tiny] * log_x[tiny] - lgamma(shape[tiny] + 1)
  out[tiny] <- if (lower.tail) lower else log(-expm1(lower))
  out
}

.positive_log_cdf <- function(family, log_y, log_mean, sd, Q = 0, lower.tail = TRUE) {
  n <- length(log_y)
  log_mean <- rep_len(log_mean, n)
  sd <- rep_len(sd, n)
  family <- .positive_family(family)
  if (family == "lognormal") {
    return(stats::pnorm(log_y, log_mean - sd^2 / 2, sd = sd, lower.tail = lower.tail, log.p = TRUE))
  }
  if (family == "gamma") {
    shape <- 1 / sd^2
    return(.log_gamma_cdf((log_y - log_mean) + log(shape), shape, lower.tail))
  }
  Q <- rep_len(Q, n)
  w <- (log_y - .gengamma_log_location(log_mean, sd, Q)) / sd
  out <- rep(NA_real_, n)
  normal_limit <- Q == 0 | (abs(Q) < 1e-7 & abs(Q * w) < .01)
  near <- which(normal_limit)
  z <- w[near]
  finite <- is.finite(z) & Q[near] != 0
  # First-order Prentice approximation retains the skew on either side
  # of zero without relying on differences of huge gamma quantiles.
  z[finite] <- z[finite] + Q[near][finite] * (z[finite]^2 + 2) / 6
  out[near] <- stats::pnorm(z, lower.tail = lower.tail, log.p = TRUE)
  for (negative in c(FALSE, TRUE)) {
    i <- which(!normal_limit & (Q < 0) == negative)
    shape <- Q[i]^-2
    out[i] <- .log_gamma_cdf(log(shape) + Q[i] * w[i], shape, xor(lower.tail, negative))
  }
  out
}

.log_gamma_relative_draw <- function(shape) {
  # Boost shapes below one before drawing so very small positive draws
  # can be retained on the log scale.
  small <- which(shape < 1)
  boosted_shape <- shape
  boosted_shape[small] <- boosted_shape[small] + 1
  draw <- log(stats::rgamma(length(shape), shape = boosted_shape) / shape)
  draw[small] <- draw[small] + log(stats::runif(length(small))) / shape[small]
  draw
}

.draw_positive_log_obs <- function(family, log_mean, sd, Q = 0) {
  n <- length(log_mean)
  sd <- rep_len(sd, n)
  family <- .positive_family(family)
  if (family == "lognormal") {
    return(stats::rnorm(n, log_mean - sd^2 / 2, sd = sd))
  }
  if (family == "gamma") {
    shape <- 1 / sd^2
    return(log_mean + .log_gamma_relative_draw(shape))
  }
  Q <- rep_len(Q, n)
  location <- .gengamma_log_location(log_mean, sd, Q)
  out <- numeric(n)
  near <- which(abs(Q) < 1e-7)
  z <- stats::rnorm(length(near))
  out[near] <- location[near] + sd[near] * (z - Q[near] * (z^2 + 2) / 6)
  i <- which(abs(Q) >= 1e-7)
  shape <- Q[i]^-2
  out[i] <- location[i] + (sd[i] / Q[i]) * .log_gamma_relative_draw(shape)
  out
}
