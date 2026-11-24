#' Observation families for catch and survey data
#'
#' The default lognormal family fits positive observations. Gamma and
#' generalized gamma are also available for positive observations. Delta and
#' Tweedie families also use sampled zeros. The observation family describes
#' variation around the predicted catch or survey index.
#'
#' @param type Link for the chance of recording a positive observation.
#'   Choose `"logit-link"` (default) or `"poisson-link"`. Both use the assessment
#'   prediction as the arithmetic mean across zeros and positives.
#' @param y50 Optional positive soft detection scale in the observation's units:
#'   the predicted level at which half of observations are expected to be zero.
#'   The default estimates a separate scale for each selected fleet or survey.
#' @param Q Optional fixed Prentice shape parameter for generalized gamma.
#'   The default estimates a separate signed shape for each selected series.
#' @param power Optional fixed Tweedie variance power strictly between 1 and 2.
#'   The default estimates a separate power for each selected fleet or survey.
#'
#' @details
#' Supply a family in `catch_settings$family` or `index_settings$family`.
#' A single family applies to the whole observation type. A named list selects
#' families by fleet or survey; unlisted series use `lognormal()`.
#' Unlabelled catches use the name `"Fishery"`.
#'
#' Use [stats::Gamma()] with `link = "log"` for positive Gamma observations.
#' The prediction has a different interpretation across families:
#'
#' | Family | Population prediction | Sampled zeros |
#' |:-------|:----------------------|:--------------|
#' | `lognormal()` | Median of positives | Excluded |
#' | Gamma, generalized gamma | Arithmetic mean of positives | Excluded |
#' | Delta, Tweedie | Arithmetic mean across zeros and positives | Included |
#'
#' Let \eqn{Y_i} be an observed catch or survey index for an age and year,
#' and let \eqn{\widetilde{Y}_i} be the population-model prediction:
#' Baranov catch for a fishery, or abundance at sampling time scaled by
#' catchability for a survey. The index \eqn{j} identifies its fleet or survey.
#' All probabilities, densities, means and medians below are conditional on
#' the model states and parameters.
#'
#' ### Delta observations
#'
#' Delta models describe whether a sample records any fish and how much it
#' records when fish are observed (Stefánsson, 1996). In tinyAM, both parts use
#' the same predicted catch or survey index, \eqn{\widetilde{Y}_i}. This
#' prediction is the arithmetic mean across zeros and positives.
#'
#' Each fleet or survey has a soft detection scale, \eqn{\widetilde{Y}_{50,j}},
#' in the same units as the prediction. It is the predicted level at which
#' half of sampled observations are expected to be zero. Let
#' \eqn{r_i=\widetilde{Y}_i/\widetilde{Y}_{50,j}} be the prediction relative to
#' that scale, \eqn{p_{+,i}=\Pr(Y_i>0)} the probability of a positive, and
#' \eqn{p_{0,i}=\Pr(Y_i=0)=1-p_{+,i}} the probability of zero.
#'
#' With `type = "logit-link"`, the odds of recording a positive equal
#' \eqn{r_i}.
#' \deqn{\operatorname{logit}(p_{+,i})=\log r_i,\qquad
#' p_{+,i}=\frac{r_i}{1+r_i}.}
#'
#' With `type = "poisson-link"`, a positive represents at least one event
#' under a Poisson model with rate \eqn{\lambda_i}. This adapts the encounter
#' model of Thorson (2018).
#' \deqn{\lambda_i=\log(2)\,r_i,\qquad p_{+,i}=1-e^{-\lambda_i}.}
#' This is equivalent to a complementary log-log (cloglog) link.
#'
#' The assessment prediction averages zeros and positives. Dividing by the
#' probability of a positive gives the mean among positive observations.
#' \deqn{m_{+,i}=\widetilde{Y}_i/p_{+,i},\qquad
#' p_{+,i}m_{+,i}=\widetilde{Y}_i.}
#' For example, if the prediction is 10 and 25% of samples are expected to
#' be positive, the positives average \eqn{10/0.25=40}. Averaging them with
#' the remaining 75% zeros gives the overall mean of 10.
#'
#' As the prediction approaches zero, this mean approaches
#' \eqn{\widetilde{Y}_{50,j}} for the logit link and
#' \eqn{\widetilde{Y}_{50,j}/\log(2)} for the Poisson link.
#' Positives can follow a lognormal, Gamma or generalized gamma distribution.
#' For delta lognormal,
#' \deqn{\log Y_i\mid Y_i>0\sim
#' N\!\left(\log m_{+,i}-\sigma_i^2/2,\sigma_i^2\right).}
#' Let \eqn{g_{+,i}(z)} be the conditional density of \eqn{\log Y_i}
#' given a positive observation. The observation log-likelihood contribution is
#' \deqn{\ell_i=\begin{cases}
#' \log p_{0,i}, & Y_i=0,\\[4pt]
#' \log p_{+,i}+\log g_{+,i}(\log Y_i), & Y_i>0.
#' \end{cases}}
#' The `q_link` setting separately controls the catchability formula.
#'
#' ### Tweedie observations
#'
#' Fitting a Tweedie model requires the suggested `tweedie` package (version
#' 3.2.0 or later) for quantile residuals and positive-observation medians.
#' Install it with `install.packages("tweedie")`.
#'
#' A compound Poisson–Gamma model represents an observation as the sum of a
#' Poisson number of positive Gamma contributions (Foster and Bravington, 2013).
#' Zero contributions give an observed zero. The Tweedie family uses the overall arithmetic mean
#' \eqn{\widetilde{Y}_i}, dispersion \eqn{\phi_i>0}, and power
#' \eqn{1<p_j<2}:
#' \deqn{\operatorname{Var}(Y_i)=\phi_i\widetilde{Y}_i^{\,p_j},\qquad
#' \Pr(Y_i=0)=\exp\!\left[-\frac{\widetilde{Y}_i^{\,2-p_j}}{\phi_i(2-p_j)}\right].}
#' Its spread and zero probability are coupled; it has no separate `y50`.
#' The positive part is a mixture of Gamma distributions, rather than one
#' Gamma distribution. Power is constant within a fleet or survey.
#' Endpoints are excluded: this interface does not fit Poisson or Gamma limits.
#'
#' Dunn and Smyth (2005) describe numerical evaluation of the Tweedie density.
#'
#' Write \eqn{f_{\mathrm{Tw}}(Y_i)} for the full Tweedie density on the
#' response scale at a positive observation, with parameters
#' \eqn{\widetilde{Y}_i,\phi_i,p_j}. Then
#' \deqn{\ell_i=\begin{cases}
#' -\frac{\widetilde{Y}_i^{\,2-p_j}}{\phi_i(2-p_j)}, & Y_i=0,\\[6pt]
#' \log f_{\mathrm{Tw}}(Y_i)+\log Y_i, & Y_i>0.
#' \end{cases}}
#' The positive Tweedie density already includes the probability of a
#' positive, so no separate \eqn{\log p_{+,i}} is added. The \eqn{\log Y_i}
#' term expresses the positive density on the log-observation scale, matching
#' the other families; zero observations retain their discrete mass.
#'
#' ### Dispersion and shape
#'
#' Observation SD settings describe log-SD for lognormal positives,
#' coefficient of variation for Gamma positives, and Prentice scale
#' \eqn{\sigma} for generalized gamma positives. For Tweedie, they specify
#' \eqn{\sqrt{\phi_i}}, so the response SD is
#' \eqn{\sqrt{\phi_i}\widetilde{Y}_i^{\,p_j/2}}. These are different scales;
#' give series with different distribution families separate dispersion coefficients.
#'
#' Tweedie dispersion depends on the observation units. Multiplying observations
#' and predictions by \eqn{c} changes \eqn{\phi_i} to \eqn{c^{2-p_j}\phi_i};
#' it is neither a coefficient of variation nor a log-SD.
#'
#' Generalized gamma uses one constant Prentice scale per series. For fixed
#' negative `Q`, supply that scale without a fitted multiplier and require
#' \eqn{1+\sigma Q>0} for a finite mean. Estimated shapes satisfy this
#' condition by construction. It becomes Gamma at \eqn{Q=\sigma}, and
#' mean-parameterized lognormal at \eqn{Q=0}. `delta_gengamma(Q = 0)` gives
#' the corresponding `delta_lognormal()`.
#'
#' ### Observation support and diagnostics
#'
#' Use zero for a sampled nondetection and `NA` for a missing observation.
#' Each fitted non-lognormal series needs a positive observation. Estimating
#' delta `y50` also needs sampled zeros; a fixed `y50` permits all-positive
#' data. Tweedie power can be weakly estimated near either endpoint.
#' Catch at structurally unfished ages is deterministic.
#'
#' For native-age delta or Tweedie data, [prepare_tam()] sums the modeled plus
#' age through the oldest biological age separately by year and series. Every
#' age needs a sampled value, including zeros; an absent age or `NA` makes the
#' terminal observation missing. The specified family applies to that sum.
#' A series without recorded values above the modeled plus age is treated as
#' already pooled. This resolution is determined across all supplied years
#' and retained within retrospective folds.
#'
#' Observed-catch hindcasts use original counts and biological-age weights.
#' With already pooled catches and older biological ages, use
#' `proj_settings = list(F_mult = 1)` in [fit_hindcast()].
#'
#' Gamma, generalized gamma, delta and Tweedie errors are conditionally
#' independent. `cor_struct = "ar1"` applies only to ordinary lognormal
#' series, including those in mixed-family fits. Standardized and quantile
#' residuals use marginal distributions and retain any across-age dependence.
#' OSA is supported only for all-lognormal fits; see [residuals.tam_fit()].
#'
#' The log forecast scores from [fit_hindcast()] use positive observations
#' and conditional-positive medians for delta and Tweedie predictions.
#' These differ from the full-distribution median, which can be zero.
#' Information criteria require matching observation records and likelihood
#' measures; see [information_criteria] for the comparisons tinyAM supports.
#'
#' Delta assumes a stable detection scale within each series. Sampling effort,
#' ageing intensity and shared sampling units can change detection or induce
#' dependence that these independent observation models do not represent.
#' Inspect zero calibration alongside positive predictions and uncertainty.
#'
#' @return A family specification for [fit_tam()] or [prepare_tam()].
#' @references
#' Stefánsson, G. (1996). Analysis of groundfish survey abundance data:
#' combining the GLM and delta approaches. ICES Journal of Marine Science.
#' \doi{10.1006/jmsc.1996.0079}
#'
#' Thorson, J. T. (2018). Three problems with the conventional delta-model for
#' biomass sampling data, and a computationally efficient alternative.
#' Canadian Journal of Fisheries and Aquatic Sciences.
#' \doi{10.1139/cjfas-2017-0266}
#'
#' Foster, S. D. and Bravington, M. V. (2013). A Poisson–Gamma model for
#' analysis of ecological non-negative continuous data. Environmental and
#' Ecological Statistics. \doi{10.1007/s10651-012-0233-0}
#'
#' Prentice, R. L. (1974). A log gamma model and its maximum likelihood
#' estimation. Biometrika. \doi{10.1093/biomet/61.3.539}
#'
#' Dunn, P. K. and Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. Statistics and Computing.
#' \doi{10.1007/s11222-005-4070-y}
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
#' delta_lognormal(type = "logit-link")
#' delta_gamma(type = "poisson-link", y50 = 100)
#' delta_gengamma(Q = 0)
#' tweedie(power = 1.5)
#' list(RV = delta_lognormal(y50 = 100), Sentinel = lognormal())
#' @name observation_families
#' @export
lognormal <- function() {
  structure(list(name = "lognormal"), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
delta_lognormal <- function(type = "logit-link", y50 = NULL) {
  .delta_family("delta_lognormal", type, y50)
}

#' @rdname observation_families
#' @export
delta_gamma <- function(type = "logit-link", y50 = NULL) {
  .delta_family("delta_gamma", type, y50)
}

#' @rdname observation_families
#' @export
gengamma <- function(Q = NULL) {
  .validate_gengamma_Q(Q)
  structure(list(name = "gengamma", Q = Q), class = "tam_observation_family")
}

#' @rdname observation_families
#' @export
delta_gengamma <- function(type = "logit-link", y50 = NULL, Q = NULL) {
  .validate_gengamma_Q(Q)
  family <- .delta_family("delta_gengamma", type, y50)
  family$Q <- Q
  family
}

#' @rdname observation_families
#' @export
tweedie <- function(power = NULL) {
  if (!is.null(power) && (!is.numeric(power) || length(power) != 1L ||
      !is.finite(power) || power <= 1 || power >= 2)) {
    cli::cli_abort("{.arg power} must be NULL or one finite number strictly between 1 and 2.")
  }
  structure(list(name = "tweedie", power = power), class = "tam_observation_family")
}

.require_tweedie <- function() {
  if (!requireNamespace("tweedie", quietly = TRUE) ||
      package_version(getNamespaceVersion("tweedie")) < "3.2.0") {
    cli::cli_abort(c(
      "Tweedie observation models require the {.pkg tweedie} package (version 3.2.0 or later).",
      "i" = "Install or update it with {.code install.packages('tweedie')}."
    ))
  }
}

.delta_family <- function(name, type, y50) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("logit-link", "poisson-link")) {
    cli::cli_abort("{.arg type} must be 'logit-link' or 'poisson-link'.")
  }
  .validate_y50(y50)
  structure(list(name = name, delta_type = type, y50 = y50), class = "tam_observation_family")
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

.is_delta_family <- function(family) family %in% c("delta_lognormal", "delta_gamma", "delta_gengamma")

.is_zero_family <- function(family) .is_delta_family(family) | family %in% "tweedie"

.is_gengamma_family <- function(family) family %in% c("gengamma", "delta_gengamma")

.positive_family <- function(family) sub("^delta_", "", family)

.is_family_spec <- function(x) {
  inherits(x, "tam_observation_family") || inherits(x, "family")
}

.normalize_family <- function(x, setting) {
  if (inherits(x, "tam_observation_family")) {
    if (!is.character(x$name) || length(x$name) != 1L || is.na(x$name)) {
      cli::cli_abort("{setting} has an invalid observation family name.")
    }
    if (identical(x$name, "lognormal")) return(lognormal())
    if (identical(x$name, "gamma")) return(x)
    if (identical(x$name, "gengamma")) return(gengamma(Q = x$Q))
    if (identical(x$name, "tweedie")) {
      spec <- tweedie(power = x$power)
      .require_tweedie()
      return(spec)
    }
    if (.is_delta_family(x$name)) {
      family <- .delta_family(x$name, x$delta_type, x$y50)
      if (identical(x$name, "delta_gengamma")) {
        .validate_gengamma_Q(x$Q)
        family$Q <- x$Q
      }
      return(family)
    }
    cli::cli_abort("{setting} has an unsupported observation family.")
  }
  if (!identical(x$family, "Gamma") || !identical(x$link, "log")) {
    cli::cli_abort("{setting} supports {.code stats::Gamma(link = 'log')} and tinyAM family constructors.")
  }
  structure(list(name = "gamma"), class = "tam_observation_family")
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

.obs_has_delta <- function(dat) length(dat$delta$series) > 0L

.obs_has_tweedie <- function(dat) length(dat$tweedie$series) > 0L

.obs_has_non_lognormal <- function(dat) {
  .obs_has_delta(dat) || any(dat$obs_map$family != "lognormal", na.rm = TRUE)
}

.prepare_observation_rows <- function(dat, d, families) {
  series <- ifelse(d$type == "catch", .observation_series(d, "catch"),
    .observation_series(d, "index"))
  d$family <- "lognormal"
  d$delta_type <- NA_character_
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

  delta_series <- gg_series <- tweedie_series <- list()
  delta_group <- gg_group <- tweedie_group <- integer(nrow(d))
  for (type in names(families)) {
    for (name in names(families[[type]])) {
      spec <- families[[type]][[name]]
      if (spec$name == "lognormal") next
      rows <- which(d$type == type & series == name)
      sample <- rows[!d$is_missing[rows] & !d$is_proj[rows] & !d$is_structural[rows]]
      positive <- sample[d$obs[sample] > 0]
      label <- paste(type, name)
      delta <- .is_delta_family(spec$name)
      if (!length(positive)) {
        family <- if (delta) "Delta" else spec$name
        cli::cli_abort("{family} series {.val {label}} needs at least one positive observation in the fitting years.")
      }
      d$family[rows] <- spec$name
      if (delta) d$delta_type[rows] <- spec$delta_type
      structural_zeros <- rows[d$is_structural[rows] & !d$is_missing[rows] & d$obs[rows] == 0]
      d$log_obs[structural_zeros] <- -Inf
      if (.is_gengamma_family(spec$name)) {
        gg_series[[length(gg_series) + 1L]] <- list(type = type, name = name,
          fixed = !is.null(spec$Q), Q = spec$Q)
        gg_group[rows] <- length(gg_series)
      }
      if (!.is_zero_family(spec$name)) next
      zeros <- sample[d$obs[sample] == 0]
      if (delta && is.null(spec$y50) && !length(zeros)) {
        cli::cli_abort(c(
          "Delta series {.val {label}} has no eligible sampled zeros.",
          "i" = "Use a positive-only family or supply a fixed {.arg y50} with {.fn {spec$name}}."
        ))
      }
      d$is_observed[sample] <- TRUE
      d$log_obs[rows[!d$is_missing[rows] & d$obs[rows] == 0]] <- -Inf
      unit <- if (is.null(d$unit)) NA_character_ else as.character(d$unit[positive[1L]])
      if (identical(spec$name, "tweedie")) {
        tweedie_series[[length(tweedie_series) + 1L]] <- list(type = type, name = name,
          unit = unit, fixed = !is.null(spec$power), power = spec$power,
          positive = length(positive), zero = length(zeros))
        tweedie_group[rows] <- length(tweedie_series)
        next
      }
      delta_series[[length(delta_series) + 1L]] <- list(type = type, name = name, unit = unit,
        family = spec$name, delta_type = spec$delta_type,
        fixed = !is.null(spec$y50), log_start = if (is.null(spec$y50)) {
          stats::median(d$log_obs[positive]) +
            if (spec$delta_type == "poisson-link") log(log(2)) else 0
        } else log(spec$y50), positive = length(positive), zero = length(zeros))
      delta_group[rows] <- length(delta_series)
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
  if (length(tweedie_series)) dat$tweedie <- list(series = tweedie_series, group = tweedie_group)
  if (length(delta_series)) {
    dat$delta <- list(series = delta_series, group = delta_group,
      lognormal = lognormal_rows,
      positive = which(d$is_observed & delta_group > 0L & d$obs > 0),
      zero = which(d$is_observed & delta_group > 0L & d$obs == 0))
  }
  dat
}

.y50_parameters <- function(dat, type) {
  specs <- Filter(function(x) x$type == type && !x$fixed, dat$delta$series)
  stats::setNames(vapply(specs, `[[`, numeric(1), "log_start"),
    vapply(specs, `[[`, character(1), "name"))
}

.y50_starts <- function(value, before, after, type) {
  old <- Filter(function(x) x$type == type && !x$fixed, before$delta$series)
  new <- Filter(function(x) x$type == type && !x$fixed, after$delta$series)
  retained <- vapply(new, function(series) {
    any(vapply(old, function(x) identical(x$name, series$name) &&
      identical(x$unit, series$unit) && identical(x$delta_type, series$delta_type), logical(1)))
  }, logical(1))
  value[intersect(names(value), vapply(new[retained], `[[`, character(1), "name"))]
}

.tidy_y50 <- function(fit, interval) {
  if (!.obs_has_delta(fit$dat)) return(NULL)
  estimate <- .tam_parameter_summary(fit, "Estimate")
  se <- .tam_parameter_summary(fit, "Std. Error")
  z <- stats::qnorm((1 + interval) / 2)
  stack_list(lapply(fit$dat$delta$series, function(s) {
    parameter <- paste0("log_y50_", s$type)
    e <- if (s$fixed) s$log_start else unname(estimate[[parameter]][s$name])
    err <- if (s$fixed) NA_real_ else unname(se[[parameter]][s$name])
    data.frame(par = paste0("y50_", s$type), coef = s$name,
      family = s$family, delta_type = s$delta_type,
      fleet = if (s$type == "catch") s$name else NA_character_,
      survey = if (s$type == "index") s$name else NA_character_,
      unit = s$unit, is_fixed = s$fixed, est = exp(e), se = err, se_scale = "log",
      lwr = exp(e - z * err), upr = exp(e + z * err))
  }), label = NULL)
}

.y50_advisories <- function(fit) {
  if (!.obs_has_delta(fit$dat)) return(NULL)
  table <- .tidy_y50(fit, 0.95)
  weak <- vapply(fit$dat$delta$series, function(s) min(s$positive, s$zero) < 5L, logical(1)) |
    !is.finite(table$se) | table$se > 1
  table <- table[weak & !table$is_fixed, , drop = FALSE]
  data.frame(issue = rep("y50_support", nrow(table)), detail = if (nrow(table)) {
    paste("Soft detection scale for", table$coef,
      "has few zeros or positives, or imprecise uncertainty. Inspect y50 and zero calibration.")
  } else character())
}

.tweedie_parameters <- function(dat, type) {
  specs <- Filter(function(x) x$type == type && !x$fixed, dat$tweedie$series)
  stats::setNames(numeric(length(specs)), vapply(specs, `[[`, character(1), "name"))
}

.tweedie_starts <- function(value, before, after, type) {
  old <- Filter(function(x) x$type == type && !x$fixed, before$tweedie$series)
  new <- Filter(function(x) x$type == type && !x$fixed, after$tweedie$series)
  retained <- vapply(new, function(series) {
    any(vapply(old, function(x) identical(x$name, series$name) &&
      identical(x$unit, series$unit), logical(1)))
  }, logical(1))
  value[intersect(names(value), vapply(new[retained], `[[`, character(1), "name"))]
}

.tweedie_from_eta <- function(eta) {
  # Keep finite parameters strictly inside the compound-Poisson interval,
  # including when the inverse logit rounds to an endpoint.
  margin <- 16 * .Machine$double.eps
  1 + margin + (1 - 2 * margin) * plogis(eta)
}

.tweedie_power <- function(dat, par, sd) {
  "[<-" <- RTMB::ADoverload("[<-")
  power <- rep(NA_real_, length(sd))
  for (j in seq_along(dat$tweedie$series)) {
    s <- dat$tweedie$series[[j]]
    rows <- which(dat$tweedie$group == j)
    power[rows] <- if (s$fixed) s$power else
      .tweedie_from_eta(par[[paste0("logit_tweedie_power_", s$type)]][[s$name]])
  }
  power
}

.tidy_tweedie_power <- function(fit, interval) {
  if (!.obs_has_tweedie(fit$dat)) return(NULL)
  estimate <- .tam_parameter_summary(fit, "Estimate")
  se <- .tam_parameter_summary(fit, "Std. Error")
  z <- stats::qnorm((1 + interval) / 2)
  stack_list(lapply(fit$dat$tweedie$series, function(s) {
    parameter <- paste0("logit_tweedie_power_", s$type)
    e <- if (s$fixed) NA_real_ else unname(estimate[[parameter]][s$name])
    err <- if (s$fixed) NA_real_ else unname(se[[parameter]][s$name])
    data.frame(par = paste0("tweedie_power_", s$type), coef = s$name,
      family = "tweedie", delta_type = NA_character_,
      fleet = if (s$type == "catch") s$name else NA_character_,
      survey = if (s$type == "index") s$name else NA_character_,
      unit = "unitless", is_fixed = s$fixed,
      est = if (s$fixed) s$power else .tweedie_from_eta(e),
      se = err, se_scale = "logit(power - 1)",
      lwr = .tweedie_from_eta(e - z * err), upr = .tweedie_from_eta(e + z * err))
  }), label = NULL)
}

.tweedie_advisories <- function(fit) {
  table <- .tidy_tweedie_power(fit, 0.95)
  if (is.null(table)) return(NULL)
  weak <- !is.finite(table$se) | table$se > 1 | table$est < 1.05 | table$est > 1.95
  table <- table[weak & !table$is_fixed, , drop = FALSE]
  data.frame(issue = rep("tweedie_power_support", nrow(table)), detail = if (nrow(table)) {
    paste("Tweedie power for", table$coef,
      "is near a boundary or has imprecise uncertainty. Inspect power, dispersion and zero calibration.")
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
.poisson_link_positive <- function(log_pred, log_y50, sd) {
  log_scale <- log_y50 - log(log(2))
  log_intensity <- log_pred - log_scale
  zero <- rep(0, length(log_intensity))
  shifted <- RTMB::logspace_add(log_intensity, zero - 35)
  safe <- exp(-RTMB::logspace_add(-shifted, zero - 36))
  log_p_shifted <- RTMB::logspace_sub(zero, -safe)
  log_mean <- log_scale + (shifted - log_p_shifted)
  list(
    log_intensity = log_intensity,
    log_p_zero = -exp(log_intensity),
    log_p_positive = log_p_shifted - RTMB::logspace_add(zero, -35 - log_intensity),
    log_mean = log_mean,
    log_median = log_mean - sd^2 / 2
  )
}

.delta_positive <- function(log_pred, log_y50, sd, type = "logit-link") {
  if (identical(type, "poisson-link")) {
    return(.poisson_link_positive(log_pred, log_y50, sd))
  }
  contrast <- log_pred - log_y50
  log_mean <- RTMB::logspace_add(log_pred, log_y50)
  list(log_p_zero = -RTMB::logspace_add(0, contrast),
    log_p_positive = -RTMB::logspace_add(0, -contrast),
    log_mean = log_mean, log_median = log_mean - sd^2 / 2)
}

# A Tweedie is a compound-Poisson mixture, with Gamma sums conditional on
# encounter count. The low-intensity expression avoids subtracting large logs.
.tweedie_components <- function(log_pred, sd, power) {
  log_scale <- 2 * log(sd) + log(2 - power) + (power - 1) * log_pred
  p <- .poisson_link_positive(log_pred, log_scale + log(log(2)), sd)
  list(log_p_zero = p$log_p_zero, log_p_positive = p$log_p_positive,
    log_mean = p$log_mean)
}

.delta_predictions <- function(dat, par, log_pred, sd, zero_catch) {
  "[<-" <- RTMB::ADoverload("[<-")
  log_p_zero <- rep(-Inf, length(log_pred))
  log_p_positive <- numeric(length(log_pred))
  log_median <- log_mean <- log_pred
  enabled <- dat$delta$group > 0L
  log_p_zero[enabled & zero_catch] <- 0
  log_p_positive[enabled & zero_catch] <- -Inf
  for (j in seq_along(dat$delta$series)) {
    s <- dat$delta$series[[j]]
    log_y50 <- if (s$fixed) s$log_start else par[[paste0("log_y50_", s$type)]][[s$name]]
    rows <- which(dat$delta$group == j & !zero_catch)
    if (!length(rows)) next
    p <- .delta_positive(log_pred[rows], log_y50, sd[rows], s$delta_type)
    log_p_zero[rows] <- p$log_p_zero
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
