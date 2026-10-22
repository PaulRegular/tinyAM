#' Merge user-supplied starting values into parameter list
#'
#' @description
#' Internal helper that overlays user-provided starting values onto the default
#' parameter list from `make_par(dat)`. Only matching names are updated.
#'
#' @details
#' Scalars are replaced directly, vectors align by names (if present) or by
#' position, and matrices align by row/column names when available, otherwise
#' by shape. Non-matching elements remain unchanged. This allows reusing
#' estimated parameters from previous fits as starting values for new runs or
#' retros, ensuring consistent initialization without manual trimming.
#'
#' @param par0 Base parameter list (usually from `make_par(dat)`).
#' @param start List of user-supplied starting values with matching structure.
#'
#' @return A parameter list like `par0`, updated with matching values from
#'   `start`.
#'
#' @keywords internal
#' @noRd

.merge_start_par <- function(par0, start) {
  for (nm in intersect(names(par0), names(start))) {
    x0 <- par0[[nm]]
    xs <- start[[nm]]

    # Scalars
    if (is.numeric(x0) && length(dim(x0)) == 0 && length(x0) == 1L && is.null(names(x0))) {
      if (is.numeric(xs) && length(xs) == 1L) par0[[nm]] <- xs
      next
    }

    # Vectors with names
    if (is.null(dim(x0))) {
      if (!is.null(names(x0)) && !is.null(names(xs))) {
        ii <- intersect(names(x0), names(xs))
        x0[ii] <- xs[ii]
      } else if (length(xs) == length(x0)) {
        x0[] <- xs
      }
      par0[[nm]] <- x0
      next
    }

    # Matrices: align by dimnames when possible, else recycle shape if equal
    if (is.matrix(x0)) {
      rnm <- rownames(x0)
      cnm <- colnames(x0)
      if (!is.null(dimnames(x0)) && !is.null(dimnames(xs))) {
        rr <- intersect(rnm, rownames(xs))
        cc <- intersect(cnm, colnames(xs))
        if (length(rr) && length(cc)) x0[rr, cc] <- xs[rr, cc]
      } else if (all(dim(xs) == dim(x0))) {
        x0[, ] <- xs
      }
      par0[[nm]] <- x0
      next
    }
  }
  # A changed recruitment lag can move a year between fixed and random states.
  recruitment <- c(start$log_r_init, start$log_r)
  if (length(recruitment) && !is.null(names(recruitment))) {
    for (nm in c("log_r_init", "log_r")) {
      years <- intersect(names(par0[[nm]]), names(recruitment))
      if (length(years)) par0[[nm]][years] <- recruitment[years]
    }
  }
  par0
}

# Read the stored conditional mode without reevaluating or changing its objective.
.start_parameters <- function(start, dat) {
  if (inherits(start, "tam_fit")) {
    source <- start
    if (!is.null(source$parameter_values) && !is.null(source$obj)) {
      start <- source$obj$env$parList(par = source$parameter_values)
    } else if (!is.null(source$sdrep)) {
      start <- as.list(source$sdrep, "Estimate")
    } else {
      cli::cli_abort("Starting fit has no stored parameter estimates.")
    }
    if (!identical(source$dat$index_settings$q_link, dat$index_settings$q_link)) {
      start <- start[!grepl("^(log_q|logit_q$|dq$|eta_q_|log_sd_q_|logit_phi_q_|q_a50_|log_q_slope_)", names(start))]
      cli::cli_inform("Catchability starts were omitted because q_link changed.")
    }
  }
  if (!is.list(start) || (length(start) && (is.null(names(start)) ||
    anyNA(names(start)) || any(!nzchar(names(start))) || anyDuplicated(names(start))))) {
    cli::cli_abort("{.arg start_par} must be a named parameter list or a {.cls tam_fit} object.")
  }
  start
}


#' Fit a Tiny Assessment Model (TAM)
#'
#' @description
#' Estimates abundance, recruitment, and mortality from catch-at-age and survey
#' observations. Returns population trends, observation diagnostics, parameter
#' estimates, and uncertainty tables. Use `N_settings`, `F_settings`, and
#' `M_settings` to specify where unexplained biological variation is allowed.
#' See [prepare_tam()] for options and [tinyAM-model] for the model equations.
#'
#' @details
#' Builds data with [prepare_tam()] and initial parameters with [make_par()].
#' RTMB integrates random effects using the Laplace approximation; `nlminb`
#' optimizes the remaining fixed effects. [mono()] increments `dq` and `dmu_f` have a zero
#' lower bound, as do initial-abundance decrements `d_n0`. Their SEs remain on the increment scale, but symmetric Wald
#' intervals are only local approximations at a boundary.
#' Set `index_settings$q_link = "logit"` to restrict the full catchability
#' prediction to between zero and one; the default `"log"` allows q above one.
#' Formula coefficients are then named `logit_q`, and [mono()] increments
#' act on logit-q. Starting q coefficients must use the selected link scale:
#' old `log_q` starts are not reused as `logit_q`. See [prepare_tam()] and [tidy_par()].
#'
#' Random-effect blocks are chosen automatically from the model settings:
#'
#' - Always includes `log_f` and eligible `log_r` states. Early recruitment
#'   anchors are fixed; see [recruitment_formulas].
#' - Includes `log_n0` only if `N_settings$init == "latent"`.
#' - Includes `log_n` if `N_settings$process != "off"`.
#' - Includes `log_m` if `M_settings$process != "off"`.
#' - Includes the unique `eta_q_*` states of Gaussian catchability formula terms.
#' - Includes Gaussian F/M mean and catch/index log-SD formula states.
#'   Logistic midpoint/slope and process SD/correlation are fixed effects.
#'
#' A warning is issued if the number of random effects exceeds 1.5 times the
#' number of observed data points (rough identifiability check).
#'
#' For F random walks, [prepare_tam()] removes time-invariant mean columns
#' because they cancel from all increments. For M random walks, mean coefficients
#' with no effect on either increments or mortality outside the fitted process
#' states are held at their starting values.
#' An RW SD is also held fixed if the process has only one row (no increments).
#' These parameters cannot be estimated from the likelihood; no constraint is
#' imposed on the latent states themselves.
#' For AR1, the correlation for an axis containing a single row or column is
#' fixed at zero; otherwise it would be confounded with the SD.
#'
#' `is_converged` summarizes optimizer, gradient and uncertainty checks. Also inspect
#' `opt$convergence` and `opt$message`; a passing numerical check does not
#' establish that the data identify every biological component.
#' Fitting also warns about sparse catchability effects and wide formula-effect
#' uncertainty. These advisories are retained by [check_tam()] and do not change
#' the likelihood, parameter constraints or numerical convergence criteria.
#'
#' Omitted entries in a settings list use package defaults; [update()] instead
#' retains the original fit's entries. Use explicit `NULL` to remove an optional
#' formula or supplied offset. Unknown or duplicate entry names are errors.
#' Supply spawning timing in `data$maturity$spawn_time`; see [prepare_tam()].
#'
#' @param data A named list of tidy observation tables (e.g., `catch`, `index`,
#'   `weight`, `maturity`) and optional `pop_covariates`. A [tam_data][as_tam_data]
#'   object or plain list can be used. See [cod_obs] for an example.
#' @param interval Level in `(0, 1)` to use to generate confidence intervals,
#'   where applicable; default `0.95.`
#' @param silent Logical; if `TRUE`, disables tracing information.
#' @param add_osa_res Logical; add one-step-ahead residuals? Hard wired to
#'   apply the `"oneStepGaussianOffMode"` method. See [RTMB::oneStepPredict()]
#'   for details.
#' @param start_par Optional fitted `tam_fit` or list of starting parameter values, matching the
#'   structure returned by [make_par()]. Useful for warm starts or retrospective
#'   runs where parameter estimates from a previous fit are reused to improve
#'   convergence and speed. Non-matching or missing entries are ignored.
#' @param grad_tol Numeric tolerance passed to [check_tam()] when
#'   evaluating the fitted object's gradients. Defaults to `1e-2`.
#'
#' @inheritParams prepare_tam
#'
#' @return
#' A `tam_fit` list with components:
#'
#' - **call**: matched call.
#' - **dat**: data list returned by [prepare_tam()].
#' - **obj**: RTMB `ADFun` object.
#' - **opt**: `[stats::nlminb()]` optimization result.
#' - **rep**: list from `obj$report()`.
#' - **sdrep**: [RTMB::sdreport()] result.
#' - **fixed_par**: fixed parameter estimates in a tidy format (see [tidy_par()]).
#' - **random_par**: list of random parameter estimates in a tidy format (see
#'                   [tidy_par()]).
#' - **is_converged**: Did numerical checks pass? (see [check_tam()])
#' - **diagnostics**: Numerical, structural and advisory checks from [check_tam()].
#' - **obs_pred**: `obs$catch` and `obs$index` data augmented with
#'                  predicted values, parameter estimates, and
#'                  standardized residuals (see [tidy_obs_pred()]).
#' - **pop**: A collection of population summaries in tidy format (see
#'            [tidy_pop()]).
#' - **rec**: Recruitment predictions and residuals, with parent pairs and
#'            curve tables for stock-recruit models (see [tidy_recruitment()]).
#' - **formula_effects**: When Gaussian q, F/M mean or observation SD formula terms are present, signed effect levels,
#'   RW increments and numeric-by contributions for reporting and dashboards.
#'
#' @example inst/examples/example_fit_default.R
#' @examples
#' fit
#' head(fit$pop$ssb)
#'
#' ## Fit with projections (status quo F)
#' fit2 <- update(fit,
#'   proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)
#' )
#'
#' if (interactive()) {
#'   fits <- list("without_proj" = fit, "with_proj" = fit2)
#'   vis_tam(fits)
#' }
#'
#' @importFrom stats nlminb rnorm
#'
#' @seealso
#' [prepare_tam()], [make_par()], [sim_tam()], [fit_retro()],
#' [RTMB::MakeADFun()], [RTMB::sdreport()]
#' @eval .settings_usage("fit_tam")
#' @export
fit_tam <- function(
  data,
  years = NULL,
  ages = NULL,
  N_settings = list(process = "iid", init = "exp", mu_form = NULL, sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL, sd = NULL, phi = NULL),
  F_settings = list(process = "rw", mu_form = NULL, sd_form = ~1, mean_ages = NULL),
  M_settings = list(
    process = "off", mu_form = NULL, mu_supplied = ~ I(0.2), sd_form = ~1,
    age_breaks = NULL, first_dev_year = NULL, mean_ages = NULL
  ),
  catch_settings = list(sd_form = ~1, sd_supplied = NULL),
  index_settings = list(sd_form = ~1, sd_supplied = NULL, q_form = ~q_block, q_link = "log"),
  proj_settings = NULL,
  interval = 0.95,
  add_osa_res = FALSE,
  silent = FALSE,
  start_par = NULL,
  grad_tol = 1e-2
) {
  call <- match.call()

  if (!is.numeric(grad_tol) || length(grad_tol) != 1L ||
    !is.finite(grad_tol) || grad_tol <= 0) {
    cli::cli_abort("{.arg grad_tol} must be one positive finite number.")
  }

  refit_args <- .complete_settings(mget(names(formals(fit_tam))))
  list2env(refit_args[.settings_names], envir = environment())
  dat <- prepare_tam(
    data = data, years = years, ages = ages,
    N_settings = N_settings, rec_settings = rec_settings, F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    proj_settings = proj_settings
  )
  par <- make_par(dat)
  formula_advisories <- rbind(
    .mean_process_advisories(dat), .q_process_advisories(dat),
    .sd_process_advisories(dat), .rec_process_advisories(dat),
    .process_sd_advisories(dat)
  )
  if (nrow(formula_advisories)) {
    cli::cli_warn(c(
      "Formula-effect caution.",
      stats::setNames(formula_advisories$detail, rep("i", nrow(formula_advisories)))
    ))
  }
  if (!is.null(start_par)) {
    start_par <- .start_parameters(start_par, dat)
    refit_args$start_par <- start_par
    par <- .merge_start_par(par, start_par)
    if (!is.null(dat$rec$curve) && !any(c("log_sr_alpha", "log_sr_beta") %in% names(start_par))) {
      par <- .initialize_rec_curve(par, dat)
    }
  }
  for (name in c("dq", "dmu_f", "d_n0")) {
    if (!is.null(par[[name]]) && any(!is.finite(par[[name]]) | par[[name]] < 0)) {
      cli::cli_abort("Starting {.arg {name}} increments must be finite and non-negative.")
    }
  }

  ran <- c("log_f", "log_r", .formula_random_parameters(dat))
  if (dat$N_settings$init == "latent") {
    ran <- c(ran, "log_n0")
  }
  if (dat$N_settings$process != "off") {
    ran <- c(ran, "log_n")
  }
  if (dat$M_settings$process != "off") {
    ran <- c(ran, "log_m")
  }
  map <- list()
  for (process in c("n", "f", "m")) {
    settings <- dat[[paste0(toupper(process), "_settings")]]
    if (!settings$process %in% c("rw", "cor_rw")) next
    states <- par[[paste0("log_", process)]]
    if (nrow(states) < 2L) {
      design <- dat$process_sd[[toupper(process)]]
      map[[design$parameter]] <- factor(rep(NA, length(par[[design$parameter]])))
    }
    if (process == "n") next
    coef_name <- if (process == "f") "log_mu_f" else "mu_m"
    if (is.null(par[[coef_name]])) next
    design <- dat[[paste0(toupper(process), "_modmat")]]
    fixed <- vapply(seq_len(ncol(design)), function(j) {
      surface <- matrix(design[, j], length(dat$years), length(dat$ages),
        dimnames = list(dat$years, dat$ages)
      )
      if (process == "f") {
        field <- surface[!dat$is_proj, , drop = FALSE]
      } else {
        field <- surface[rownames(states), settings$age_block_start, drop = FALSE]
        surface[rownames(states), names(settings$age_blocks)] <- 0
        if (any(surface != 0)) {
          return(FALSE)
        } # Mean still controls M outside process states.
      }
      nrow(field) < 2L || all(field[-1, , drop = FALSE] == field[-nrow(field), , drop = FALSE])
    }, logical(1))
    if (any(fixed)) {
      indices <- seq_along(fixed)
      indices[fixed] <- NA_integer_
      map[[coef_name]] <- factor(indices)
    }
  }
  for (process in c("n", "f", "m")) {
    if (dat[[paste0(toupper(process), "_settings")]]$process != "ar1") next
    states <- par[[paste0("log_", process)]]
    singleton <- c(ncol(states), nrow(states)) == 1L # age, year
    if (any(singleton)) {
      nm <- paste0("logit_phi_", process)
      par[[nm]][singleton] <- qlogis(0)
      indices <- 1:2
      indices[singleton] <- NA_integer_
      map[[nm]] <- factor(indices)
    }
  }

  make_nll_fun <- function(f, d) function(p) f(p, d) # use closure to avoid global assignment of data
  structural <- .check_tam_structure(dat, map)
  if (any(structural$status == "fail")) {
    cli::cli_abort(c("Model formulas contain redundant fitted coefficients.",
      "x" = paste(structural$check[structural$status == "fail"], collapse = ", ")
    ))
  }
  obj <- RTMB::MakeADFun(
    make_nll_fun(nll_fun, dat),
    par,
    random = ran,
    map = map,
    silent = silent
  )

  lower <- rep(-Inf, length(obj$par))
  lower[names(obj$par) %in% c("dq", "dmu_f", "d_n0")] <- 0
  opt <- tryCatch(stats::nlminb(
    obj$par, obj$fn, obj$gr,
    lower = lower,
    control = list(eval.max = 1000, iter.max = 1000)
  ), error = function(e) cli::cli_abort("Optimization failed: {conditionMessage(e)}"))
  # Synchronize reports and the conditional random states with the returned solution.
  obj$fn(opt$par)
  gradient <- obj$gr(opt$par)
  parameter_values <- obj$env$last.par
  rep <- obj$report(parameter_values)
  sdreport_error <- NULL
  sdrep <- tryCatch(
    RTMB::sdreport(obj,
      par.fixed = opt$par,
      getJointPrecision = TRUE, getReportCovariance = FALSE
    ),
    error = function(e) {
      sdreport_error <<- conditionMessage(e)
      NULL
    }
  )

  out <- list(
    call = call,
    refit_args = refit_args,
    dat = dat,
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep,
    sdreport_error = sdreport_error,
    gradient = gradient,
    parameter_values = parameter_values,
    parameter_map = map,
    bounds = list(lower = lower, upper = rep(Inf, length(lower))),
    grad_tol = grad_tol
  )

  class(out) <- c("tam_fit", "list")

  par_tabs <- tidy_par(out, interval = interval)
  out$fixed_par <- par_tabs$fixed
  out$random_par <- par_tabs$random
  out$obs_pred <- tidy_obs_pred(out, add_osa_res = add_osa_res, interval = interval, trace = !silent)
  out$pop <- tidy_pop(out, interval = interval)
  out$rec <- tidy_recruitment(out, interval = interval)
  if (length(.formula_terms(dat))) out$formula_effects <- .tidy_formula_effects(out, interval = interval)
  out$diagnostics <- check_tam(out)
  out$is_converged <- out$diagnostics$is_converged
  if (!out$is_converged) {
    failed <- out$diagnostics$numerical
    failed <- failed$detail[failed$status != "pass"]
    cli::cli_warn(c("Numerical convergence checks did not pass; inspect {.code check_tam(fit)}.",
      stats::setNames(failed, rep("x", length(failed))),
      "i" = "Optimizer code zero alone is insufficient. A missed gradient tolerance alone does not establish a structural problem."
    ))
  }
  formula_uncertainty <- rbind(
    .formula_uncertainty_advisories(out), .rec_fit_advisories(out),
    .cor_rw_advisories(out), .process_sd_fit_advisories(out)
  )
  if (nrow(formula_uncertainty)) {
    cli::cli_warn(c(
      "Formula-effect uncertainty or support is limited.",
      stats::setNames(formula_uncertainty$detail, rep("i", nrow(formula_uncertainty)))
    ))
  }

  .new_tam_fit(out)
}


#' Update a fitted TAM model
#'
#' Refit a TAM model after changing one or more arguments. Arguments not
#' supplied in `...` are taken from the evaluated inputs stored in the
#' original fit, so updating does not depend on objects remaining available
#' in the original calling environment.
#'
#' Settings lists are updated by entry name: omitted entries retain their original
#' values. Formulas and vectors are replaced in full. Explicit `NULL` clears
#' optional entries; `proj_settings = NULL` disables projections. An empty list
#' changes nothing. Initial fits instead fill omitted settings from package defaults.
#'
#' @param object A fitted `tam_fit` object.
#' @param ... Named arguments to change in the original fit. For settings lists,
#'   supply only the entries to change, e.g. `N_settings = list(init = "free")`.
#' @param evaluate Logical; if `TRUE`, fit and return the updated model.
#'   If `FALSE`, return the updated call without fitting.
#'
#' @return An updated `tam_fit` object when `evaluate = TRUE`; otherwise,
#'   the updated call.
#'
#' @export
update.tam_fit <- function(object, ..., evaluate = TRUE) {
  object <- .require_tam_fit(object, arg = "object")

  changes <- list(...)

  if (length(changes) &&
    (is.null(names(changes)) || anyNA(names(changes)) ||
      any(names(changes) == "") || anyDuplicated(names(changes)))) {
    cli::cli_abort("All arguments supplied to {.fn update} must have unique, nonempty names.")
  }

  # Use requested inputs, before preparation changes offsets or designs.
  args <- object$refit_args
  change_expr <- as.list(match.call(expand.dots = FALSE)$...)
  new_call <- as.list(object$call)
  for (name in names(changes)) {
    value <- if (name %in% .settings_names) {
      .merge_settings(args[[name]], changes[[name]], name)
    } else {
      changes[[name]]
    }
    args[name] <- list(value)
    new_call[name] <- if (name %in% .settings_names) list(value) else change_expr[name]
  }
  new_call <- as.call(new_call)
  if (!evaluate) {
    return(new_call)
  }

  out <- do.call(fit_tam, args)

  # Preserve the concise user-facing call.
  out$call <- new_call

  out
}

