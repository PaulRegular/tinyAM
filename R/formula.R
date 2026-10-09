#' Random catchability effects in formulas
#'
#' @description
#' Allow catchability to vary among groups or through time without fitting an
#' unrelated fixed coefficient for every value. Use these markers as additive
#' terms in `index_settings$q_form`; ordinary formula terms supply the baseline.
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
#'
#' Forecast states follow the same normalized process density and are integrated
#' out when there are no observations. At the fitted mode, new IID levels have
#' zero effect, RW holds its last level, and AR1 decays toward zero. Simulation
#' includes future process variation. These are conditional predictions on the
#' link scale, not marginal response-scale means.
#'
#' @param x Column defining levels (IID) or ordered steps (RW/AR1).
#' @param by Optional column: numeric multiplier or categorical groups.
#' @param sd Optional positive fixed process SD. `NULL` estimates it.
#' @param phi Optional fixed AR1 correlation in `[0, 1)`. `NULL` estimates it.
#' @return Formula markers; direct calls raise an error.
#' @examples
#' ~ survey + iid(year, by = survey)
#' ~ survey + rw(year, by = survey)
#' ~ survey + ar1(year, by = survey)
#' ~ survey + (1 | vessel)
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
  cli::cli_abort("{term}() is only supported as an additive term in {.arg index_settings$q_form}.")
}

.formula_call_name <- function(x) {
  if (!is.call(x)) return("")
  fun <- x[[1L]]
  if (is.symbol(fun)) return(as.character(fun))
  if (is.call(fun) && as.character(fun[[1L]]) %in% c("::", ":::") &&
      identical(fun[[2L]], as.name("tinyAM"))) return(as.character(fun[[3L]]))
  ""
}

.parse_q_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) cli::cli_abort("{.arg q_form} must be a formula.")
  special <- c("iid", "rw", "ar1", "logistic", "|")
  contains <- function(x) {
    is.call(x) && (.formula_call_name(x) %in% special ||
      any(vapply(as.list(x)[-1L], contains, logical(1))))
  }
  if (!contains(formula)) return(.parse_mono_q_formula(formula, data))
  if (length(formula) != 2L) cli::cli_abort("Structured {.arg q_form} must be one-sided.")
  specs <- list()
  column <- function(x, label) {
    if (!is.symbol(x) || !as.character(x) %in% names(data)) {
      cli::cli_abort("Structured {.arg {label}} must name an existing index column.")
    }
    as.character(x)
  }
  fixed <- function(x, label, upper = Inf, allow_zero = FALSE) {
    if (is.null(x) || identical(x, quote(NULL))) return(NULL)
    value <- tryCatch(eval(x, envir = environment(formula)), error = function(e) NULL)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || (!allow_zero && value == 0) || value >= upper) {
      cli::cli_abort("{label} must be one finite {if (allow_zero) 'non-negative' else 'positive'} number below {upper}.")
    }
    value
  }
  strip <- function(x) {
    if (!contains(x)) return(x)
    type <- .formula_call_name(x)
    if (type %in% special) {
      if (type == "|") {
        if (length(x) != 3L || !identical(x[[2L]], 1)) {
          cli::cli_abort("Only random intercept syntax {.code (1 | group)} is supported.")
        }
        spec <- list(type = "iid", variable = column(x[[3L]], "group"), by = NULL,
                     sd = NULL, phi = NULL)
      } else {
        fun <- switch(type, iid = iid, rw = rw, ar1 = ar1, logistic = logistic)
        x[[1L]] <- as.name(type)
        args <- tryCatch(as.list(match.call(fun, x))[-1L], error = function(e)
          cli::cli_abort("Invalid {type}() arguments: {conditionMessage(e)}"))
        spec <- list(type = type, variable = column(args$x, "x"),
          by = if (is.null(args$by) || identical(args$by, quote(NULL))) NULL else column(args$by, "by"),
          sd = fixed(args$sd, "sd"), phi = fixed(args$phi, "phi", 1, TRUE))
      }
      specs[[length(specs) + 1L]] <<- spec
      return(NULL)
    }
    if (type == "(") return(strip(x[[2L]]))
    if (type == "+") {
      parts <- Filter(Negate(is.null), lapply(as.list(x)[-1L], strip))
      if (!length(parts)) return(NULL)
      if (length(parts) == 1L) return(parts[[1L]])
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
      return(list(group = group, rows = rows, midpoint = stats::median(x[support]),
                  slope = 4 / diff(range(x[support]))))
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
    list(group = group, rows = rows, levels = levels, n_historical = n,
         index = match(x[rows], levels), states = paste(group, state_levels, sep = ":"),
         is_proj = !state_levels %in% past_levels)
  })
  spec$id <- id
  spec$groups <- groups
  spec$multiplier <- multiplier
  spec$parameter <- if (spec$type == "logistic") NULL else paste0("eta_q_", id)
  spec$sd_parameter <- if (is.null(spec$sd) && spec$type != "logistic") paste0("log_sd_q_", id) else NULL
  spec$phi_parameter <- if (spec$type == "ar1" && is.null(spec$phi)) paste0("logit_phi_q_", id) else NULL
  spec
}

.q_term_parameters <- function(terms) {
  parameters <- list()
  for (term in terms) {
    labels <- vapply(term$groups, `[[`, character(1), "group")
    if (term$type == "logistic") {
      parameters[[paste0("q_a50_", term$id)]] <- stats::setNames(
        vapply(term$groups, `[[`, numeric(1), "midpoint"), labels)
      parameters[[paste0("log_q_slope_", term$id)]] <- stats::setNames(log(
        vapply(term$groups, `[[`, numeric(1), "slope")), labels)
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
  if (!length(dat$q_terms)) return(invisible(NULL))
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
          cli::cli_abort("{id} is indistinguishable from estimated observation error. Use replicated levels or supply one SD.")
        }
      }
    }
  }
  independent <- vapply(gaussian, function(term) term$type == "iid" ||
    (term$type == "ar1" && identical(term$phi, 0)), logical(1))
  candidates <- names(designs)[independent & vapply(gaussian,
    function(term) is.null(term$sd), logical(1))]
  if (length(candidates) > 1L) {
    partition <- function(z) {
      index <- max.col(abs(z), ties.method = "first")
      list(group = match(index, unique(index)),
           weight = z[cbind(seq_len(nrow(z)), index)])
    }
    for (i in seq_len(length(candidates) - 1L)) for (j in (i + 1L):length(candidates)) {
      a <- designs[[candidates[i]]][rows, , drop = FALSE]
      b <- designs[[candidates[j]]][rows, , drop = FALSE]
      if (identical(partition(a), partition(b))) {
        cli::cli_abort("{candidates[i]} and {candidates[j]} define the same IID variance component. Retain one term.")
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
  for (term in dat$q_terms) {
    if (term$type == "logistic") {
      a50 <- par[[paste0("q_a50_", term$id)]]
      slope <- exp(par[[paste0("log_q_slope_", term$id)]])
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
    } else term$phi
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
        nll <- nll - RTMB::dnorm(x[1L], 0, sd, log = TRUE)
        if (length(x) > 1L) nll <- nll - sum(RTMB::dnorm(
          x[-1L] - x[-length(x)], 0, sd, log = TRUE))
        level <- numeric(length(x) + 1L)
        level[-1L] <- x
        x <- level
      }
      if (term$type == "ar1") {
        nll <- nll - RTMB::dnorm(x[1L], 0, sd / sqrt(1 - phi^2), log = TRUE)
        if (length(x) > 1L) nll <- nll - sum(RTMB::dnorm(
          x[-1L], phi * x[-length(x)], sd, log = TRUE))
      }
      contribution[g$rows] <- contribution[g$rows] + x[g$index] * term$multiplier[g$rows]
      offset <- offset + length(g$states)
    }
    simulated[[term$parameter]] <- states
  }
  list(contribution = contribution, log_selectivity = selectivity,
       nll = nll, parameters = simulated)
}

#' Rising catchability curves in formulas
#'
#' @description
#' Use these terms in `index_settings$q_form` when catchability should increase
#' with age or size. `mono()` fits flexible non-decreasing steps, including exact
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
#'
#' @param x Name of a column in index observations: numeric or factor for
#'   `mono()`, numeric age/size for `logistic()`.
#' @param by Optional name of a categorical grouping column.
#' @return A formula marker; calling this function directly raises an error.
#' @seealso [formula_effects], [prepare_tam()], [make_par()], [tidy_obs_pred()], [tidy_par()], [tinyAM-model]
#' @examples
#' ~ q_block # ordinary unconstrained q
#' ~ mono(q_block)
#' ~ survey + mono(q_block, by = survey)
#' ~ survey + logistic(age, by = survey)
#' @name catchability_curves
#' @export
mono <- function(x, by = NULL) {
  cli::cli_abort("mono() is only supported as an additive term in {.arg index_settings$q_form}.")
}

.parse_mono_q_formula <- function(formula, data) {
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
