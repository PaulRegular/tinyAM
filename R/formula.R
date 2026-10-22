#' Random effects in model formulas
#'
#' @description
#' Allow catchability to vary among groups or through time without fitting an
#' unrelated fixed coefficient for every value. Use these markers as additive
#' terms in `index_settings$q_form`, `F_settings$mu_form`, or
#' `M_settings$mu_form`, `catch_settings$sd_form`, or `index_settings$sd_form`;
#' ordinary formula terms supply the baseline.
#' Recruitment uses a separate named residual process in `rec_settings`;
#' see [recruitment_formulas] for its initial-state convention.
#'
#' @details
#' `iid(x)` gives each level an independent, mean-zero Normal effect.
#' `(1 | group)` is shorthand for `iid(group)`. `rw(x)` allows changes to persist:
#' its first effect is zero and successive increments are independent Normal
#' variables. `ar1(x)` describes fluctuations that tend to return to zero; its
#' first state has the stationary Normal distribution. Estimated AR1 correlation
#' is between zero and one. Integer numeric coordinates preserve calendar gaps;
#' ordered factors use their declared order. RWs and AR1 need ordered coordinates.
#'
#' Gaussian effects are added on the selected log/logit scale, not the q scale.
#' Numeric `by` multiplies one shared trajectory by that covariate. Factor or
#' character `by` gives independent group trajectories sharing one SD and, for
#' AR1, one correlation. IID/AR1 effects are not forcibly centered. SD is estimated
#' unless supplied: marginal SD for IID, increment/innovation SD for RW/AR1.
#'
#' Terms must be additive, with bare column names. Unsupported aliases and
#' observation-specific IID effects indistinguishable from estimated observation
#' error are rejected. For sparse data, simpler formulas or supplied process SDs
#' may be necessary. Passing design checks does not establish biological
#' identifiability; inspect [check_tam()].
#' With a logit link, fitting warns when observation-specific IID variation and
#' observation SD are both estimated: the nonlinear link does not guarantee
#' useful separation. Groups with fewer than five observed effect levels are
#' also flagged when process SD or correlation is estimated. These are
#' descriptive cautions, not additional model restrictions. Projection years,
#' missing/zero observations and zero numeric multipliers do not supply support.
#'
#' Forecast states follow the same normalized process density and are integrated
#' out when there are no observations. At the fitted mode, new IID levels have
#' zero effect, RW holds its last level, and AR1 decays toward zero. Simulation
#' includes future process variation. These are conditional predictions on the
#' link scale, not marginal response-scale means.
#'
#' For F and M means, effects are added to log mean mortality. The absolute
#' latent states remain `log_f` and `log_m`; their residuals are calculated
#' relative to the resulting mean surfaces. A shared `rw(year)` or `ar1(year)`
#' can describe a population-wide change while `process = "iid"` represents
#' independent residual fluctuations by age (F) or age block (M).
#' M also supports structured means with `process = "off"`.
#' Existing supplied M is an offset: `~ 0 + rw(year)` modifies it without
#' estimating another baseline. Effects must be constant within each fitted
#' M age block. Use a categorical age column for separate trajectories;
#' numeric `by = age` instead multiplies one shared trajectory by age.
#'
#' Initially, one temporal mean term is permitted per F/M surface, with IID
#' residuals or M residuals off. Structured means with RW/AR1 residuals are
#' rejected pending evidence that both temporal components can be separated,
#' except for the age-group M random slopes described below.
#' Two estimated IID variances on the same states are also rejected.
#' Known design redundancies are checked, but sparse observations and the
#' trade-off between F, M, recruitment and catchability can still make a
#' numerically converged model weakly identified. Use recovery simulations
#' before interpreting separate process components.
#' Fitting warns when M mean and IID residual SDs are both estimated, or when
#' both F and M mean SDs are estimated. Recovery experiments found weak
#' separation of small M components and imprecise AR1 persistence; begin with
#' mean-only M and add residual variation only when supported by the data.
#' [check_tam()] retains these cautions and flags wide SD/correlation intervals.
#' Warnings do not change the fitted model or establish non-identifiability.
#'
#' ## M covariate slopes with AR1 residuals
#'
#' `M_settings$process = "ar1"` also permits one
#' `iid(age_group, by = covariate)` mean term. Each age-defined group receives
#' one slope, reused across years; the numeric covariate multiplies it.
#' Store both columns in `data$pop_covariates`. Group membership must be constant
#' through time at each age. Supplied M, ordinary terms and random slopes must
#' all agree within each M state block: `age_breaks` still pools actual M,
#' not just its unexplained deviations. Use finer breaks for distinct age slopes.
#'
#' For example, `~ 0 + capelin + iid(M_age_group, by = capelin)` gives
#' \deqn{\log M_{t,b} = \log M^{supplied}_{t,b}
#'       + (\beta + u_{g(b)})x_{t,b} + \eta_{t,b},
#'       \qquad u_g \sim N(0,\tau^2).}
#' The residual matrix \eqn{\eta} retains the existing stationary age/year AR1
#' distribution and innovation-scale SD. Random slopes are independent of it.
#' The ordinary covariate term estimates the common slope \eqn{\beta}; omit
#' that term to centre slopes on zero. No common slope, covariate centring,
#' sign restriction or covariate transformation is added automatically.
#' Supplied covariates remain external inputs, not functions recalculated from
#' the fitted population. Simulation draws one slope per group before drawing
#' the M states. Forecast covariates retain the existing terminal-value rule.
#'
#' This combination is experimental. Persistent covariates may be difficult to
#' distinguish from persistent M deviations. Fitting and [check_tam()] flag it,
#' and flag fewer than five informative slope groups when estimating their SD.
#' Several ages sharing one M state do not increase the number of independent
#' groups. Inspect uncertainty and compare a common fixed-slope model.
#' With a positive-definite uncertainty calculation, [check_tam()] also applies
#' the existing wide-SD and AR1-interval cautions to these M residual parameters.
#' Other Gaussian means with AR1 residuals, multiple random-slope terms and the
#' corresponding F extension remain unsupported.
#'
#' In observation SD formulas, Gaussian effects modify log SD:
#' \eqn{\log \sigma_i = \log \sigma_{supplied,i} + X_i\beta + z_i}.
#' A positive effect increases residual spread; it does not change the expected
#' log observation. For example, `~ iid(age)` partially pools SDs across ages,
#' while `~ survey + rw(age, by = survey)` gives separate smooth age patterns.
#' SD remains positive. Supplied SDs are offsets and remove the ordinary
#' intercept, as for fixed SD formulas. The effect's process SD describes
#' variation in log SD, not observation error itself. There is no requirement
#' that the age curve rise at either tail; a quadratic fixed age term is also
#' a useful simpler alternative when that pattern is expected.
#'
#' One RW/AR1 term is initially allowed per observation SD surface. `mono()`
#' and `logistic()` remain catchability curves. Replication across years at each
#' age helps estimate residual spread; a few ages give little information about
#' AR1 persistence. Fitting warns about levels with fewer than three observations
#' and estimated processes with fewer than five observed levels per group.
#' These descriptive cautions exclude missing/zero observations and projections.
#' Observation-specific random SDs are scale mixtures, not the same variance
#' component as an additive Gaussian random mean; they are permitted with a
#' sparse-support warning. SD effects are integrated using RTMB's Laplace
#' approximation, so recovery checks are especially useful with sparse data.
#' This extension does not introduce formulas for latent N/F/M process SDs.
#' Estimated catch-SD age random walks also trigger an advisory: full-model
#' recovery simulations frequently failed when F was estimated jointly, despite
#' good isolated recovery. Compare IID age effects or a fixed age curve and
#' inspect numerical diagnostics. This caution does not forbid the model.
#'
#' @param x Column in the relevant observation table defining levels (IID) or
#'   ordered steps (RW/AR1).
#' @param by Optional column: numeric multiplier or categorical groups.
#' @param sd Optional positive fixed process SD. `NULL` estimates it.
#' @param phi Optional fixed AR1 correlation in `[0, 1)`. `NULL` estimates it.
#' @return Formula markers; direct calls raise an error.
#' @examples
#' ~ survey + iid(year, by = survey)
#' ~ survey + rw(year, by = survey)
#' ~ survey + ar1(year, by = survey)
#' ~ survey + (1 | vessel)
#' # F_settings = list(process = "iid", mu_form = ~ factor(age) + rw(year))
#' # M_settings = list(process = "off", mu_supplied = ~ I(0.2),
#' #                   mu_form = ~ 0 + ar1(year))
#' # M_settings = list(process = "ar1", mu_supplied = ~ I(0.3),
#' #   mu_form = ~ 0 + capelin + iid(M_age_group, by = capelin))
#' # catch_settings = list(sd_form = ~ iid(age))
#' # index_settings = list(q_form = ~ survey,
#' #                       sd_form = ~ survey + ar1(age, by = survey))
#' @seealso [catchability_curves], [prepare_tam()], [fit_tam()], [check_tam()]
#' @name formula_effects
#' @export
iid <- function(x, by = NULL, sd = NULL) .formula_marker_error("iid")

#' @rdname formula_effects
#' @export
rw <- function(x, by = NULL, sd = NULL) .formula_marker_error("rw")

#' @rdname formula_effects
#' @export
ar1 <- function(x, by = NULL, sd = NULL, phi = NULL) .formula_marker_error("ar1")

#' @rdname catchability_curves
#' @export
logistic <- function(x, by = NULL) .formula_marker_error("logistic")

.formula_marker_error <- function(term) {
  cli::cli_abort("{term}() is a formula marker. Use Gaussian terms in q_form, F/M mu_form, or catch/index sd_form; logistic() belongs in q_form or F mu_form.")
}

.formula_call_name <- function(x) {
  if (!is.call(x)) {
    return("")
  }
  fun <- x[[1L]]
  if (is.symbol(fun)) {
    return(as.character(fun))
  }
  if (is.call(fun) && as.character(fun[[1L]]) %in% c("::", ":::") &&
    identical(fun[[2L]], as.name("tinyAM"))) {
    return(as.character(fun[[3L]]))
  }
  ""
}

.parse_q_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) cli::cli_abort("{.arg q_form} must be a formula.")
  if (any(c("bh", "ricker") %in% all.names(formula))) {
    cli::cli_abort("bh() and ricker() are restricted to rec_settings$mu_form.")
  }
  special <- c("iid", "rw", "ar1", "logistic", "|")
  contains <- function(x) {
    is.call(x) && (.formula_call_name(x) %in% special ||
      any(vapply(as.list(x)[-1L], contains, logical(1))))
  }
  if (!contains(formula)) {
    return(.parse_mono_q_formula(formula, data))
  }
  if (length(formula) != 2L) cli::cli_abort("Structured {.arg q_form} must be one-sided.")
  specs <- list()
  column <- function(x, label) {
    if (!is.symbol(x) || !as.character(x) %in% names(data)) {
      cli::cli_abort("Structured {.arg {label}} must name an existing observation column.")
    }
    as.character(x)
  }
  fixed <- function(x, label, upper = Inf, allow_zero = FALSE) {
    if (is.null(x) || identical(x, quote(NULL))) {
      return(NULL)
    }
    value <- tryCatch(eval(x, envir = environment(formula)), error = function(e) NULL)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || (!allow_zero && value == 0) || value >= upper) {
      cli::cli_abort("{label} must be one finite {if (allow_zero) 'non-negative' else 'positive'} number below {upper}.")
    }
    value
  }
  strip <- function(x) {
    if (!contains(x)) {
      return(x)
    }
    type <- .formula_call_name(x)
    if (type %in% special) {
      if (type == "|") {
        if (length(x) != 3L || !identical(x[[2L]], 1)) {
          cli::cli_abort("Only random intercept syntax {.code (1 | group)} is supported.")
        }
        spec <- list(
          type = "iid", variable = column(x[[3L]], "group"), by = NULL,
          sd = NULL, phi = NULL
        )
      } else {
        fun <- switch(type,
          iid = iid,
          rw = rw,
          ar1 = ar1,
          logistic = logistic
        )
        x[[1L]] <- as.name(type)
        args <- tryCatch(as.list(match.call(fun, x))[-1L], error = function(e) {
          cli::cli_abort("Invalid {type}() arguments: {conditionMessage(e)}")
        })
        spec <- list(
          type = type, variable = column(args$x, "x"),
          by = if (is.null(args$by) || identical(args$by, quote(NULL))) NULL else column(args$by, "by"),
          sd = fixed(args$sd, "sd"), phi = fixed(args$phi, "phi", 1, TRUE)
        )
      }
      specs[[length(specs) + 1L]] <<- spec
      return(NULL)
    }
    if (type == "(") {
      return(strip(x[[2L]]))
    }
    if (type == "+") {
      parts <- Filter(Negate(is.null), lapply(as.list(x)[-1L], strip))
      if (!length(parts)) {
        return(NULL)
      }
      if (length(parts) == 1L) {
        return(parts[[1L]])
      }
      return(as.call(c(list(as.name("+")), parts)))
    }
    if (type == "-" && length(x) == 3L && !contains(x[[3L]])) {
      left <- strip(x[[2L]])
      return(call("-", if (is.null(left)) 1 else left, x[[3L]]))
    }
    cli::cli_abort("Structured effects must be additive, without interactions or transformations.")
  }
  ordinary <- formula
  rhs <- strip(formula[[2L]])
  ordinary[[2L]] <- if (is.null(rhs)) 1 else rhs
  design <- .parse_mono_q_formula(ordinary, data)
  terms <- lapply(specs, .make_q_term, data = data)
  ids <- vapply(terms, `[[`, character(1), "id")
  if (anyDuplicated(ids)) cli::cli_abort("Duplicate structured catchability terms are not identifiable.")
  design$q_terms <- stats::setNames(terms, ids)
  design
}

.make_q_term <- function(spec, data) {
  x <- data[[spec$variable]]
  if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x) ||
    (is.numeric(x) && any(!is.finite(x)))) {
    cli::cli_abort("Structured {.arg x} must have finite, non-missing values.")
  }
  temporal <- spec$type %in% c("rw", "ar1")
  if (temporal && !(is.ordered(x) ||
    (is.numeric(x) && all(x == floor(x))))) {
    cli::cli_abort("RW/AR1 coordinates must be integer numeric values or ordered factors.")
  }
  if (spec$type == "logistic" && !is.numeric(x)) {
    cli::cli_abort("logistic() {.arg x} must be numeric.")
  }
  by <- if (is.null(spec$by)) rep("all", nrow(data)) else data[[spec$by]]
  if (!is.atomic(by) || !is.null(dim(by)) || anyNA(by) ||
    (is.numeric(by) && any(!is.finite(by)))) {
    cli::cli_abort("Structured {.arg by} must have finite, non-missing values.")
  }
  numeric_by <- !is.null(spec$by) && is.numeric(by)
  if (numeric_by && spec$type == "logistic") {
    cli::cli_abort("logistic() only supports categorical {.arg by} groups.")
  }
  multiplier <- if (numeric_by) by else rep(1, nrow(data))
  if (numeric_by) by <- rep("all", nrow(data))
  by <- as.character(by)
  historical <- if (is.null(data$is_proj)) rep(TRUE, nrow(data)) else !data$is_proj
  observed <- historical
  if (!is.null(data$obs)) observed <- observed & is.finite(data$obs) & data$obs > 0
  id <- paste(c(spec$type, spec$variable, spec$by), collapse = "_")
  groups <- lapply(unique(by), function(group) {
    rows <- which(by == group)
    past <- rows[historical[rows]]
    support <- rows[observed[rows] & multiplier[rows] != 0]
    if (!length(support)) cli::cli_abort("{id}: group {group} has no informative historical observations.")
    if (spec$type == "logistic") {
      if (length(unique(x[support])) < 3L) {
        cli::cli_abort("logistic() needs at least three observed x values per group.")
      }
      return(list(
        group = group, rows = rows, midpoint = stats::median(x[support]),
        slope = 4 / diff(range(x[support]))
      ))
    }
    if (temporal && is.numeric(x)) {
      levels <- seq.int(min(x[past]), max(x[rows]))
      past_levels <- levels[levels <= max(x[past])]
    } else if (temporal) {
      levels <- levels(x)[seq_len(max(as.integer(x[rows])))]
      past_levels <- levels[seq_len(max(as.integer(x[past])))]
    } else {
      levels <- if (is.factor(x)) levels(droplevels(x[rows])) else sort(unique(x[rows]))
      past_levels <- levels[levels %in% x[past]]
      levels <- c(past_levels, setdiff(levels, past_levels))
    }
    n <- length(past_levels)
    if (n < 2L && is.null(spec$sd)) {
      cli::cli_abort("{id}: estimating an SD requires at least two historical levels; supply {.arg sd} or simplify.")
    }
    if (spec$type == "rw" && n < 2L) cli::cli_abort("rw() needs at least two historical steps.")
    if (spec$type == "ar1" && n < 3L && is.null(spec$phi)) {
      cli::cli_abort("Estimating AR1 correlation requires at least three historical steps; supply {.arg phi}.")
    }
    state_levels <- if (spec$type == "rw") levels[-1L] else levels
    list(
      group = group, rows = rows, levels = levels, n_historical = n,
      index = match(x[rows], levels), states = paste(group, state_levels, sep = ":"),
      is_proj = !state_levels %in% past_levels
    )
  })
  spec$id <- id
  spec$groups <- groups
  spec$multiplier <- multiplier
  spec$parameter <- if (spec$type == "logistic") NULL else paste0("eta_q_", id)
  spec$sd_parameter <- if (is.null(spec$sd) && spec$type != "logistic") paste0("log_sd_q_", id) else NULL
  spec$phi_parameter <- if (spec$type == "ar1" && is.null(spec$phi)) paste0("logit_phi_q_", id) else NULL
  if (spec$type == "logistic") {
    spec$a50_parameter <- paste0("q_a50_", id)
    spec$slope_parameter <- paste0("log_q_slope_", id)
  }
  spec
}

.q_term_parameters <- function(terms) {
  parameters <- list()
  for (term in terms) {
    labels <- vapply(term$groups, `[[`, character(1), "group")
    if (term$type == "logistic") {
      parameters[[term$a50_parameter]] <- stats::setNames(
        vapply(term$groups, `[[`, numeric(1), "midpoint"), labels
      )
      parameters[[term$slope_parameter]] <- stats::setNames(log(
        vapply(term$groups, `[[`, numeric(1), "slope")
      ), labels)
    } else {
      states <- unlist(lapply(term$groups, `[[`, "states"), use.names = FALSE)
      parameters[[term$parameter]] <- stats::setNames(numeric(length(states)), states)
      if (!is.null(term$sd_parameter)) parameters[[term$sd_parameter]] <- c(process = log(.1))
      if (!is.null(term$phi_parameter)) parameters[[term$phi_parameter]] <- c(process = stats::qlogis(.5))
    }
  }
  parameters
}

.q_random_parameters <- function(dat) {
  unlist(lapply(dat$q_terms, `[[`, "parameter"), use.names = FALSE)
}

.q_term_design <- function(term, n) {
  states <- sum(vapply(term$groups, function(g) length(g$states), integer(1)))
  z <- matrix(0, n, states)
  offset <- 0L
  for (g in term$groups) {
    index <- g$index - as.integer(term$type == "rw")
    rows <- g$rows[index > 0L]
    z[cbind(rows, offset + index[index > 0L])] <- term$multiplier[rows]
    offset <- offset + length(g$states)
  }
  z
}

.check_q_terms <- function(dat) {
  if (!length(dat$q_terms)) {
    return(invisible(NULL))
  }
  d <- dat$obs$index
  rows <- !d$is_proj & is.finite(d$obs) & d$obs > 0
  fixed <- dat$q_modmat[rows, , drop = FALSE]
  if (!is.null(dat$q_mono_modmat)) fixed <- cbind(fixed, dat$q_mono_modmat[rows, , drop = FALSE])
  gaussian <- Filter(function(term) term$type != "logistic", dat$q_terms)
  designs <- lapply(gaussian, .q_term_design, n = nrow(d))
  for (term in Filter(function(term) term$type == "logistic", dat$q_terms)) {
    jacobian <- matrix(0, nrow(d), 2L * length(term$groups))
    for (i in seq_along(term$groups)) {
      g <- term$groups[[i]]
      x <- d[[term$variable]][g$rows] - g$midpoint
      remaining <- 1 - stats::plogis(g$slope * x)
      jacobian[g$rows, 2L * i - 1L] <- -g$slope * remaining
      jacobian[g$rows, 2L * i] <- g$slope * x * remaining
    }
    combined <- cbind(fixed, jacobian[rows, , drop = FALSE])
    scale <- sqrt(colSums(combined^2))
    scale[scale == 0] <- 1
    if (qr(sweep(combined, 2, scale, "/"))$rank < ncol(combined)) {
      cli::cli_abort("{term$id} duplicates the ordinary or mono age/size curve. Remove redundant curve terms.")
    }
  }
  for (id in names(designs)) {
    term <- gaussian[[id]]
    z <- designs[[id]][rows, , drop = FALSE]
    informative <- colSums(abs(z)) > 0
    z <- z[, informative, drop = FALSE]
    if (!ncol(z)) cli::cli_abort("{id} has no informative effect states.")
    if (is.null(term$sd) && ncol(fixed) && qr(cbind(fixed, z))$rank == qr(fixed)$rank) {
      cli::cli_abort("{id} duplicates saturated fixed effects. Remove those fixed terms or supply {.arg sd}.")
    }
    independent <- term$type == "iid" || (term$type == "ar1" && identical(term$phi, 0))
    if (independent && is.null(term$sd) && all(colSums(z != 0) == 1L)) {
      diagonal <- rowSums(z^2)
      active <- diagonal > 0
      # Constant observation-specific variance is exactly the same component
      # as a freely estimated observation SD on those rows (for the log link).
      if (identical(dat$index_settings$q_link, "log") &&
        length(unique(diagonal[active])) == 1L) {
        s <- dat$sd_index_modmat[rows, , drop = FALSE]
        if (ncol(s) && qr(cbind(s, as.numeric(active)))$rank == qr(s)$rank) {
          noise <- if (is.null(dat$mean_component)) {
            "observation error"
          } else {
            paste("IID", dat$mean_component, "residual variance")
          }
          cli::cli_abort("{id} is indistinguishable from estimated {noise}. Use replicated levels or supply one SD.")
        }
      }
    }
  }
  independent <- vapply(gaussian, function(term) {
    term$type == "iid" ||
      (term$type == "ar1" && identical(term$phi, 0))
  }, logical(1))
  candidates <- names(designs)[independent & vapply(
    gaussian,
    function(term) is.null(term$sd), logical(1)
  )]
  if (length(candidates) > 1L) {
    partition <- function(z) {
      index <- max.col(abs(z), ties.method = "first")
      list(
        group = match(index, unique(index)),
        weight = z[cbind(seq_len(nrow(z)), index)]
      )
    }
    for (i in seq_len(length(candidates) - 1L)) {
      for (j in (i + 1L):length(candidates)) {
        a <- designs[[candidates[i]]][rows, , drop = FALSE]
        b <- designs[[candidates[j]]][rows, , drop = FALSE]
        if (identical(partition(a), partition(b))) {
          cli::cli_abort("{candidates[i]} and {candidates[j]} define the same IID variance component. Retain one term.")
        }
      }
    }
  }
  invisible(NULL)
}

.q_effects <- function(par, dat, simulate = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")
  contribution <- selectivity <- numeric(nrow(dat$obs$index))
  nll <- 0
  simulated <- list()
  rw_terms <- Filter(function(term) term$type == "rw", dat$q_terms)
  increments <- numeric(sum(vapply(rw_terms, function(term) {
    sum(vapply(term$groups, function(g) length(g$states), integer(1)))
  }, integer(1))))
  increment_offset <- 0L
  for (term in dat$q_terms) {
    if (term$type == "logistic") {
      a50 <- par[[term$a50_parameter]]
      slope <- exp(par[[term$slope_parameter]])
      for (i in seq_along(term$groups)) {
        rows <- term$groups[[i]]$rows
        predictor <- slope[i] * (dat$obs$index[[term$variable]][rows] - a50[i])
        selectivity[rows] <- selectivity[rows] - RTMB::logspace_add(0, -predictor)
      }
      next
    }
    sd <- if (is.null(term$sd)) exp(par[[term$sd_parameter]]) else term$sd
    phi <- if (is.null(term$phi)) {
      if (term$type == "ar1") plogis(par[[term$phi_parameter]]) else 0
    } else {
      term$phi
    }
    states <- par[[term$parameter]]
    offset <- 0L
    for (g in term$groups) {
      ii <- offset + seq_along(g$states)
      x <- states[ii]
      if (simulate) {
        if (term$type == "iid") x[] <- stats::rnorm(length(x), 0, sd)
        if (term$type == "rw") x[] <- cumsum(stats::rnorm(length(x), 0, sd))
        if (term$type == "ar1") {
          x[1L] <- stats::rnorm(1L, 0, sd / sqrt(1 - phi^2))
          if (length(x) > 1L) for (j in 2:length(x)) x[j] <- stats::rnorm(1L, phi * x[j - 1L], sd)
        }
        states[ii] <- x
      }
      if (term$type == "iid") nll <- nll - sum(RTMB::dnorm(x, 0, sd, log = TRUE))
      if (term$type == "rw") {
        delta <- x
        if (length(x) > 1L) delta[-1L] <- x[-1L] - x[-length(x)]
        increments[increment_offset + seq_along(delta)] <- delta
        increment_offset <- increment_offset + length(delta)
        nll <- nll - RTMB::dnorm(x[1L], 0, sd, log = TRUE)
        if (length(x) > 1L) {
          nll <- nll - sum(RTMB::dnorm(
            x[-1L] - x[-length(x)], 0, sd,
            log = TRUE
          ))
        }
        level <- numeric(length(x) + 1L)
        level[-1L] <- x
        x <- level
      }
      if (term$type == "ar1") {
        nll <- nll - RTMB::dnorm(x[1L], 0, sd / sqrt(1 - phi^2), log = TRUE)
        if (length(x) > 1L) {
          nll <- nll - sum(RTMB::dnorm(
            x[-1L], phi * x[-length(x)], sd,
            log = TRUE
          ))
        }
      }
      contribution[g$rows] <- contribution[g$rows] + x[g$index] * term$multiplier[g$rows]
      offset <- offset + length(g$states)
    }
    simulated[[term$parameter]] <- states
  }
  list(
    contribution = contribution, log_selectivity = selectivity,
    nll = nll, parameters = simulated, rw_increments = increments
  )
}

.q_level_metadata <- function(term) {
  all_levels <- unique(unlist(lapply(term$groups, `[[`, "levels"), use.names = FALSE))
  do.call(rbind, lapply(term$groups, function(g) {
    level <- as.character(g$levels)
    data.frame(
      term = term$id, process = term$type, variable = term$variable,
      by = if (is.null(term$by)) NA_character_ else term$by, group = g$group,
      level = level, coordinate = if (is.numeric(g$levels)) g$levels else match(g$levels, all_levels),
      is_proj = seq_along(level) > g$n_historical,
      is_fixed = term$type == "rw" & seq_along(level) == 1L
    )
  }))
}

.tidy_q_effects <- function(fit, interval = .95, increment_report = "eta_q_increments") {
  terms <- Filter(function(term) term$type != "logistic", fit$dat$q_terms)
  if (!length(terms)) {
    return(list())
  }
  estimate <- .tam_parameter_summary(fit, "Estimate")
  error <- .tam_parameter_summary(fit, "Std. Error")
  z <- stats::qnorm(.5 + interval / 2)
  report_errors <- if (is.null(fit[["sdrep"]])) {
    NULL
  } else {
    as.list(fit[["sdrep"]], "Std. Error", report = TRUE)
  }
  out <- list(levels = list(), increments = list(), contributions = list())
  offset <- 0L
  for (term in terms) {
    tab <- .q_level_metadata(term)
    tab$est <- tab$se <- numeric(nrow(tab))
    tab$est[!tab$is_fixed] <- estimate[[term$parameter]]
    tab$se[!tab$is_fixed] <- error[[term$parameter]]
    tab$lwr <- tab$est - z * tab$se
    tab$upr <- tab$est + z * tab$se
    tab$se_scale <- "reported"
    out$levels[[term$id]] <- tab
    if (term$type == "rw") {
      delta <- tab[!tab$is_fixed, , drop = FALSE]
      ii <- offset + seq_len(nrow(delta))
      delta$est <- fit$rep[[increment_report]][ii]
      delta$se <- if (is.null(report_errors)) NA_real_ else report_errors[[increment_report]][ii]
      delta$lwr <- delta$est - z * delta$se
      delta$upr <- delta$est + z * delta$se
      out$increments[[term$id]] <- delta
      offset <- offset + nrow(delta)
    }
    if (!is.null(term$by) && is.numeric(fit$dat$obs$index[[term$by]])) {
      d <- fit$dat$obs$index
      rows <- unlist(lapply(term$groups, `[[`, "rows"), use.names = FALSE)
      index <- unlist(lapply(term$groups, `[[`, "index"), use.names = FALSE)
      contribution <- cbind(
        d[rows, intersect(c("year", "age", "survey", "is_proj"), names(d)), drop = FALSE],
        tab[index, c("term", "process", "variable", "group", "level", "coordinate", "est", "lwr", "upr", "se", "se_scale")]
      )
      contribution$multiplier <- term$multiplier[rows]
      contribution$est <- contribution$est * contribution$multiplier
      limits <- cbind(contribution$lwr, contribution$upr) * contribution$multiplier
      contribution$lwr <- apply(limits, 1, min)
      contribution$upr <- apply(limits, 1, max)
      contribution$se <- contribution$se * abs(contribution$multiplier)
      out$contributions[[term$id]] <- unique(contribution)
    }
  }
  attr(out, "interval") <- interval
  out
}

#' Rising catchability and mean fishing-mortality curves
#'
#' @description
#' Use these terms in `index_settings$q_form` when catchability should increase
#' with age or size, or in `F_settings$mu_form` for a rising long-term mean F.
#' `mono()` fits flexible non-decreasing steps, including exact
#' plateaus. `logistic()` fits a smooth rising selectivity curve. Ordinary formula
#' terms supply baseline catchability; neither curve describes a dome with other
#' covariates held constant.
#'
#' @details
#' Both functions are additive formula markers, not numeric transformations.
#'
#' ## Flexible steps: mono()
#'
#' Numeric values are ordered
#' increasingly; factors (including ordered factors) use their declared levels.
#' Supply factor levels in the scientifically intended order.
#'
#' The first represented level has no increment. Subsequent link-scale levels add
#' cumulative non-negative increments `dq`, fitted directly with a lower bound
#' of zero by [fit_tam()]. These are increments on the log-q scale by default,
#' or logit-q when `index_settings$q_link = "logit"`, not absolute q.
#' Both links preserve non-decreasing q. Steps start at 0.05.
#' A zero step gives an exact plateau between separate
#' levels at a finite parameter value; pooled levels also retain identical q.
#' When optimizing [nll_fun()] directly, supply the same zero lower bounds for
#' `dq`. Standard errors describe local curvature on the increment scale;
#' inference at an active boundary need not follow a symmetric normal law.
#'
#' `by` gives each group independent steps, using only levels represented in
#' that group, in their declared order. Each group needs at least two levels.
#' Baselines come from ordinary formula terms: `~ mono(x, by = survey)` shares
#' an intercept, whereas `~ survey + mono(x, by = survey)` gives separate
#' survey baselines. Other ordinary covariates can modify q, so monotonicity
#' holds with those covariates held constant. Interactions involving `mono()`,
#' transformed arguments, and ordinary effects of the same `x` are unsupported.
#'
#' ## Smooth selectivity: logistic()
#'
#' For numeric age or size, the curve is
#' \eqn{S(x) = 1 / (1 + \exp\{-k(x-a_{50})\})}, with positive slope \eqn{k}.
#' The midpoint \eqn{a_{50}} is the age/size at half the asymptotic selectivity.
#' This curve multiplies catchability after the selected inverse link, rather
#' than adding a sigmoid on the link scale. Thus logit q remains below one.
#' Ordinary terms supply its baseline: `~ survey + logistic(age, by = survey)`
#' gives survey-specific baselines and curves. Each curve needs at least three
#' observed x values. Numeric `by` and interactions are unsupported.
#' Midpoint and positive slope are fixed effects, not Gaussian random effects.
#' Redundant curves, such as logistic plus an unrestricted factor of the same
#' age, are rejected. See [formula_effects] for variation around these curves.
#' Fitting warns if the midpoint lies outside the observed age/size range and
#' its 95% interval is wider than that range. Maximum catchability and curve
#' shape may then be difficult to distinguish. Uncertainty advisories require
#' a successful positive-definite uncertainty calculation; failures are reported
#' separately by the numerical convergence checks.
#'
#' ## Mean fishing mortality
#'
#' With `F_settings$process = "off"`, `"iid"` or `"ar1"`, use `~ mono(age)` or
#' `~ logistic(age)`. Ordinary terms provide the baseline F. One curve is allowed;
#' its numeric coordinate must be age or an age-derived `pop_covariates` column
#' constant across years at each age. Categorical `by` groups may have separate
#' curves. M mean curves are not enabled.
#'
#' For mono, \eqn{\log\mu_F = X\beta + C d}, where C contains cumulative
#' indicators and the directly fitted `dmu_f` increments satisfy \eqn{d\geq0}.
#' For logistic, \eqn{\log\mu_F = X\beta + \log S(x)}. The baseline therefore
#' gives the asymptotic F when other covariates are held constant. Midpoint
#' `F_a50_*` and positive slope `log_F_slope_*` are fitted separately from q
#' curves. With the residual process off, F equals this mean. Otherwise realized
#' F fluctuates around it and need not be monotonic.
#' `mu_F` is the exponential of expected log F, rather than an arithmetic
#' average over lognormal process draws.
#' Static curves are rejected with RW/correlated-RW residuals because they
#' cancel from the increments. Gaussian structured mean terms retain their
#' existing IID-residual or no-residual-process restriction, including with a curve.
#' `dmu_f` has the same zero bounds and boundary-inference limitations as `dq`.
#' [tidy_pop()] reports mean-F profiles with joint RTMB-derived uncertainty.
#'
#' @param x Name of an index column for q: numeric or factor for `mono()`,
#'   numeric age/size for `logistic()`. F mean curves require numeric age or
#'   an age-derived population covariate.
#' @param by Optional name of a categorical grouping column.
#' @return A formula marker; calling this function directly raises an error.
#' @seealso [formula_effects], [prepare_tam()], [make_par()], [tidy_obs_pred()], [tidy_par()], [tinyAM-model]
#' @examples
#' ~q_block # ordinary unconstrained q
#' ~ mono(q_block)
#' ~ survey + mono(q_block, by = survey)
#' ~ survey + logistic(age, by = survey)
#' # F_settings = list(process = "iid", mu_form = ~ mono(age))
#' # F_settings = list(process = "ar1", mu_form = ~ logistic(age))
#' @name catchability_curves
#' @export
mono <- function(x, by = NULL) {
  cli::cli_abort("mono() is only supported as an additive term in index_settings$q_form or F_settings$mu_form.")
}

.parse_mono_q_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg q_form} must be a formula.")
  }
  contains_mono <- function(x) {
    if (!is.call(x)) {
      return(FALSE)
    }
    if (identical(x[[1L]], as.name("mono"))) {
      return(TRUE)
    }
    if (as.character(x[[1L]])[1L] %in% c("::", ":::") &&
      identical(x[[3L]], as.name("mono"))) {
      return(TRUE)
    }
    any(vapply(as.list(x), contains_mono, logical(1)))
  }
  # Preserve the original model.matrix path, including contrasts and attributes.
  if (!contains_mono(formula)) {
    return(list(q_modmat = stats::model.matrix(formula, data = data)))
  }
  if (length(formula) != 2L) cli::cli_abort("Monotonic {.arg q_form} must be one-sided.")
  specs <- list()
  strip <- function(x) {
    if (!contains_mono(x)) {
      return(x)
    }
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
    if (is.call(x) && identical(x[[1L]], as.name("("))) {
      return(strip(x[[2L]]))
    }
    if (is.call(x) && identical(x[[1L]], as.name("-")) && length(x) == 3L &&
      !contains_mono(x[[3L]])) {
      left <- strip(x[[2L]])
      return(call("-", if (is.null(left)) 1 else left, x[[3L]]))
    }
    if (is.call(x) && identical(x[[1L]], as.name("+"))) {
      args <- lapply(as.list(x)[-1L], strip)
      args <- Filter(Negate(is.null), args)
      if (!length(args)) {
        return(NULL)
      }
      if (length(args) == 1L) {
        return(args[[1L]])
      }
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
  list(
    q_modmat = q_modmat,
    q_mono_modmat = do.call(cbind, lapply(design, `[[`, "matrix")),
    q_mono_steps = do.call(rbind, lapply(design, `[[`, "steps"))
  )
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
    label <- paste0(
      spec$variable,
      if (!is.null(spec$by)) paste0("[", spec$by, "=", groups[i], "]"),
      ":", utils::head(lev, -1L), "->", utils::tail(lev, -1L)
    )
    colnames(mat) <- label
    matrices[[i]] <- mat
    steps[[i]] <- data.frame(
      coef = label, variable = spec$variable,
      from_level = as.character(utils::head(lev, -1L)), to_level = as.character(utils::tail(lev, -1L)),
      by = if (is.null(spec$by)) NA_character_ else spec$by,
      by_level = if (is.null(spec$by)) NA_character_ else as.character(groups[i])
    )
  }
  list(matrix = do.call(cbind, matrices), steps = do.call(rbind, steps))
}


# Mean effects use the same unique-state compiler and densities as q effects.
.mean_formula_data <- function(dat, component) {
  d <- if (component == "F") dat$F_data else .population_data(dat)
  d$is_proj <- dat$is_proj[match(d$year, dat$years)]
  d$obs <- 1
  d
}

.parse_mean_formula <- function(formula, dat, component) {
  .population_data(dat, formula, paste0(component, "_settings$mu_form"))
  d <- .mean_formula_data(dat, component)
  design <- .parse_q_formula(formula, d)
  terms <- design$q_terms
  logistic <- Filter(function(x) x$type == "logistic", terms)
  mono <- design$q_mono_steps
  n_mono <- if (is.null(mono)) 0L else length(unique(mono$variable))
  has_curve <- n_mono + length(logistic) > 0L
  if (has_curve) {
    if (component != "F") cli::cli_abort("M_settings$mu_form does not support catchability curves: mono() and logistic() are enabled for F only.")
    if (n_mono + length(logistic) != 1L) cli::cli_abort("Use one mono() or logistic() curve in F_settings$mu_form.")
    if (!dat$F_settings$process %in% c("off", "iid", "ar1")) {
      cli::cli_abort("F mean curves require process = 'off', 'iid' or 'ar1'. Static age curves cancel from RW increments.")
    }
    variable <- if (n_mono) mono$variable[1L] else logistic[[1L]]$variable
    x <- d[[variable]]
    if (!is.numeric(x) || any(vapply(split(x, d$age), function(z) length(unique(z)) != 1L, logical(1)))) {
      cli::cli_abort("F mean curves must use numeric age or an age-derived population covariate constant across years at each age.")
    }
    by <- if (n_mono) mono$by[1L] else logistic[[1L]]$by
    if (!is.null(by) && !is.na(by) && is.numeric(d[[by]])) cli::cli_abort("F mean curves require categorical by groups.")
    if (!ncol(design$q_modmat)) cli::cli_abort("F mean curves require an ordinary baseline; retain the intercept or supply group baselines.")
    fixed <- design$q_modmat[!d$is_proj, , drop = FALSE]
    if (qr(cbind(fixed, 1))$rank > qr(fixed)$rank) cli::cli_abort("F mean curves require an ordinary baseline; retain the intercept or supply group baselines.")
    if (n_mono) {
      combined <- cbind(design$q_modmat, design$q_mono_modmat)
      combined <- combined[!d$is_proj, , drop = FALSE]
      if (qr(combined)$rank < ncol(combined)) cli::cli_abort("The mono F curve duplicates ordinary mean terms. Remove redundant age effects or supply group baselines.")
    }
  }
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term$id <- paste(component, term$id, sep = "_")
    if (term$type == "logistic") {
      term$a50_parameter <- paste0("F_a50_", term$id)
      term$slope_parameter <- paste0("log_F_slope_", term$id)
    } else {
      term$parameter <- paste0("eta_mu_", term$id)
    }
    if (!is.null(term$sd_parameter)) term$sd_parameter <- paste0("log_sd_mu_", term$id)
    if (!is.null(term$phi_parameter)) term$phi_parameter <- paste0("logit_phi_mu_", term$id)
    terms[[i]] <- term
  }
  if (length(terms)) names(terms) <- vapply(terms, `[[`, character(1), "id")
  list(
    matrix = design$q_modmat, terms = terms,
    mono_matrix = design$q_mono_modmat, mono_steps = mono
  )
}

# A slope has one age-defined coefficient reused across all modeled years.
.M_ar1_random_slope <- function(dat) {
  terms <- dat$M_terms
  if (!identical(dat$M_settings$process, "ar1") || length(terms) != 1L) return(FALSE)
  term <- terms[[1L]]
  if (term$type != "iid" || is.null(term$by)) return(FALSE)
  d <- .mean_formula_data(dat, "M")
  is.numeric(d[[term$by]]) && all(vapply(split(d[[term$variable]], d$age),
    function(x) length(unique(x)) == 1L, logical(1)))
}

.check_mean_terms <- function(dat, component) {
  terms <- dat[[paste0(component, "_terms")]]
  if (!length(terms)) {
    return(invisible(NULL))
  }
  settings <- dat[[paste0(component, "_settings")]]
  if (sum(vapply(terms, function(x) x$type %in% c("rw", "ar1"), logical(1))) > 1L) {
    cli::cli_abort(c("Use one temporal process in the {component} mean.",
      "i" = "Multiple RW/AR1 mean terms may explain the same changes in mortality; their separate contributions have not been validated.",
      "i" = "Retain one rw() or ar1() term and use IID residuals for independent fluctuations."
    ))
  }
  gaussian <- Filter(function(x) x$type != "logistic", terms)
  random_slope <- component == "M" && .M_ar1_random_slope(dat)
  if (length(gaussian) && !settings$process %in% c("iid", "off") && !random_slope) {
    alternative <- paste0(", or use process = 'off' to vary only the ", component, " mean")
    cli::cli_abort(c("Structured {component} means currently require IID residuals or process = off.",
      "i" = "Structured means and temporal residual processes may explain overlapping mortality variation; this combination has not been validated.",
      "i" = "Set process = 'iid'{alternative}.",
      if (component == "M") c("i" = "With AR1 M residuals, one iid(age_group, by = numeric_covariate) term is also supported. Its groups must be constant across years and respect M age blocks.")
    ))
  }
  d <- .mean_formula_data(dat, component)
  fixed <- dat[[paste0(component, "_modmat")]]
  proxy <- list(
    obs = list(index = d), q_terms = terms, q_modmat = fixed,
    q_mono_modmat = if (component == "F") dat$F_mono_modmat else NULL,
    mean_component = component,
    index_settings = list(q_link = "log"),
    sd_index_modmat = matrix(
      if (settings$process == "iid") 1 else numeric(), nrow(d),
      if (settings$process == "iid") 1L else 0L
    )
  )
  # M ages sharing one absolute process state supply one replicate, not many.
  if (component == "M" && settings$process == "iid") {
    for (term in terms) {
      if (term$type != "iid" || !is.null(term$sd)) next
      z <- .q_term_design(term, nrow(d))
      effective <- d$year %in% settings$years & d$age %in% settings$age_block_start
      outside <- !d$age %in% as.numeric(names(settings$age_blocks)) |
        !d$year %in% settings$years
      # A free mean variance must have replication beyond a single M state.
      active <- z[effective | outside, , drop = FALSE]
      if (all(colSums(active != 0) <= 1L)) {
        cli::cli_abort(c("Unreplicated IID M mean effects with IID residuals are not enabled with both SDs estimated.",
          "i" = "A shared M age block counts as one state, even when it contains several ages.",
          "i" = "Share mean effects across independent age blocks, supply sd in iid(), or use M_settings$process = 'off'."
        ))
      }
    }
  }
  .check_q_terms(proxy)
}

# These are cautions about supported designs, not proof of non-identifiability.
.mean_process_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  estimated <- function(terms) {
    any(vapply(terms, function(term) !is.null(term$sd_parameter), logical(1)))
  }
  if (.M_ar1_random_slope(dat)) {
    out <- rbind(out, data.frame(issue = "M_random_slopes_AR1", detail = paste(
      "Age-group random slopes with AR1 M residuals are an experimental combination.",
      "Persistent covariates and persistent residuals may explain similar changes in M; check slope and variance uncertainty and compare a fixed-slope model."
    )))
    term <- dat$M_terms[[1L]]
    if (!is.null(term$sd_parameter)) {
      d <- .mean_formula_data(dat, "M")
      rows <- !d$is_proj & d$year %in% dat$M_settings$years &
        d$age %in% dat$M_settings$age_block_start
      groups <- sum(colSums(abs(.q_term_design(term, nrow(d))[rows, , drop = FALSE])) > 0)
      if (groups < 5L) {
        out <- rbind(out, data.frame(issue = "M_random_slope_groups", detail = paste(
          "Only", groups, "independent slope groups are represented in fitted M states.",
          "Their slope SD may be weakly estimated; ages sharing an M state do not provide additional independent groups."
        )))
      }
    }
  }
  if (identical(dat$M_settings$process, "iid") && estimated(dat$M_terms)) {
    out <- rbind(out, data.frame(issue = "M_variance_separation", detail = paste(
      "Both M mean-process and IID residual SDs are estimated. Simulations showed that small shared M variation can be difficult to separate from residual variation.",
      "Start with M_settings$process = 'off', or supply sd in the mean term when external information supports it; check recovery before interpreting both components."
    )))
  }
  if (estimated(dat$F_terms) && estimated(dat$M_terms)) {
    out <- rbind(out, data.frame(issue = "joint_mortality_means", detail = paste(
      "F and M mean-process SDs are both estimated. The recovery study varied one mortality mean at a time; joint estimation has not been validated.",
      "Compare with models varying only one mean and examine sensitivity of F, M and abundance."
    )))
  }
  out
}

.mean_effects <- function(par, dat, component, simulate = FALSE) {
  proxy <- list(
    obs = list(index = .mean_formula_data(dat, component)),
    q_terms = dat[[paste0(component, "_terms")]]
  )
  out <- .q_effects(par, proxy, simulate = simulate)
  if (component == "F") {
    out$contribution <- out$contribution + out$log_selectivity
    if (!is.null(dat$F_mono_modmat)) {
      out$contribution <- out$contribution + drop(dat$F_mono_modmat %*% par$dmu_f)
    }
  }
  out
}

.formula_terms <- function(dat) {
  c(dat$q_terms, dat$F_terms, dat$M_terms, dat$sd_catch_terms, dat$sd_index_terms)
}

.formula_random_parameters <- function(dat) {
  unlist(lapply(.formula_terms(dat), `[[`, "parameter"), use.names = FALSE)
}

.tidy_formula_effects <- function(fit, interval = .95) {
  all <- list(q = .tidy_q_effects(fit, interval))
  for (component in c("F", "M", "sd_catch", "sd_index")) {
    if (!length(fit$dat[[paste0(component, "_terms")]])) next
    proxy <- fit
    proxy$dat$q_terms <- fit$dat[[paste0(component, "_terms")]]
    proxy$dat$obs$index <- if (startsWith(component, "sd_")) {
      fit$dat$obs[[sub("^sd_", "", component)]]
    } else {
      .mean_formula_data(fit$dat, component)
    }
    all[[component]] <- .tidy_q_effects(proxy, interval,
      increment_report = paste0(
        "eta_", if (startsWith(component, "sd_")) "" else "mu_",
        component, "_increments"
      )
    )
  }
  out <- list(levels = list(), increments = list(), contributions = list())
  for (component in names(all)) {
    for (part in names(out)) {
      tables <- lapply(all[[component]][[part]], function(tab) {
        tab$component <- component
        tab
      })
      out[[part]] <- c(out[[part]], tables)
    }
  }
  if (!length(out$levels)) {
    return(list())
  }
  attr(out, "interval") <- interval
  out
}


# Observation SD effects reuse the Gaussian formula-state compiler.
.parse_sd_formula <- function(formula, data, component) {
  design <- .parse_q_formula(formula, data)
  terms <- design$q_terms
  if (!is.null(design$q_mono_modmat) ||
    any(vapply(terms, function(x) x$type == "logistic", logical(1)))) {
    cli::cli_abort("{component}_settings$sd_form supports iid(), rw() and ar1(), not catchability curves.")
  }
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term$id <- paste("sd", component, term$id, sep = "_")
    term$parameter <- paste0("eta_", term$id)
    if (!is.null(term$sd_parameter)) term$sd_parameter <- paste0("log_", term$id)
    if (!is.null(term$phi_parameter)) term$phi_parameter <- paste0("logit_phi_", term$id)
    terms[[i]] <- term
  }
  if (length(terms)) names(terms) <- vapply(terms, `[[`, character(1), "id")
  list(matrix = design$q_modmat, terms = terms)
}

.sd_effects <- function(par, dat, component, simulate = FALSE) {
  .q_effects(par, list(
    obs = list(index = dat$obs[[component]]),
    q_terms = dat[[paste0("sd_", component, "_terms")]]
  ), simulate = simulate)
}

.check_sd_terms <- function(dat, component) {
  terms <- dat[[paste0("sd_", component, "_terms")]]
  if (!length(terms)) {
    return(invisible(NULL))
  }
  if (sum(vapply(terms, function(x) x$type %in% c("rw", "ar1"), logical(1))) > 1L) {
    cli::cli_abort("Use one RW/AR1 term in {component}_settings$sd_form; overlapping ordered SD processes have not been validated.")
  }
  # A random effect on SD is not additive observation variance on the mean scale.
  .check_q_terms(list(
    obs = list(index = dat$obs[[component]]), q_terms = terms,
    q_modmat = dat[[paste0("sd_", component, "_modmat")]],
    index_settings = list(q_link = NULL)
  ))
}

.sd_process_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  for (component in c("catch", "index")) {
    terms <- dat[[paste0("sd_", component, "_terms")]]
    if (!length(terms)) next
    d <- dat$obs[[component]]
    observed <- !d$is_proj & is.finite(d$obs) & d$obs > 0
    for (term in terms) {
      if (component == "catch" && term$type == "rw" && term$variable == "age" &&
        !is.null(term$sd_parameter)) {
        out <- rbind(out, data.frame(issue = "catch_SD_RW", detail = paste0(
          term$id,
          ": estimating an age random walk in catch SD can be numerically fragile when F is also estimated. Full assessment recovery simulations frequently failed despite repeated observations at each age. Compare IID age effects or a fixed age curve, and inspect check_tam(fit)."
        )))
      }
      for (group in term$groups) {
        rows <- group$rows[observed[group$rows] & term$multiplier[group$rows] != 0]
        support <- table(d[[term$variable]][rows])
        support <- support[support > 0]
        if (any(support < 3L)) {
          out <- rbind(out, data.frame(
            issue = "SD_effect_replication",
            detail = paste0(
              term$id, ", group ", group$group,
              ": some SD effect levels have fewer than three observations. Random SD effects describe residual spread, not mean trends; sparse levels may be weakly estimated. Compare a simpler SD formula and inspect uncertainty."
            )
          ))
        }
        if ((!is.null(term$sd_parameter) || !is.null(term$phi_parameter)) && length(support) < 5L) {
          out <- rbind(out, data.frame(issue = "SD_effect_support", detail = paste0(
            term$id,
            ", group ", group$group, ": only ", length(support),
            " observed levels inform the SD process. Its process SD or AR1 correlation may be imprecise; consider a simpler formula or externally supported sd/phi."
          )))
        }
      }
    }
  }
  out
}
