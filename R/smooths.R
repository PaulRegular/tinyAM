#' Penalized smooth F and catchability relationships
#'
#' Use `s()` or `t2()` in `F_settings$mu_form` or `index_settings$q_form`,
#' including inside [for_fleet()] and [for_survey()]. The basis is built by
#' mgcv; the assessment jointly estimates the relationship and its smoothness.
#' No separate GAM is fitted to catch or survey observations.
#'
#' @details Supported terms are univariate `s(x, k, bs, by)` with thin-plate
#' (`"tp"`) or cubic-regression (`"cr"`) bases, and `t2(age, year, k, bs, by)`
#' using these marginal bases. Predictors and categorical `by` groups must be
#' bare input-column names. Other mgcv options, numeric `by`, `te()`, `ti()`,
#' multidimensional `s()`, linked penalties, fixed smoothness, random-effect
#' bases and custom constraints are not supported. Choose a modest `k` supported
#' by distinct predictor values; it limits basis size, not the fitted EDF.
#'
#' Unpenalized basis directions are fixed coefficients. Each penalized block
#' has coefficients \eqn{u_j\sim N(0,\tau_j^2 I)} with its complete normalized
#' Gaussian density. Smoothing SDs are estimated through RTMB and coefficients
#' are integrated by Laplace approximation. Tensor terms can have several
#' independent penalties; factor-by smooths have separate penalties per level.
#' These penalties describe mean relationships, not F residual innovations.
#' Contributions enter the log F or selected q-link predictor, before any
#' multiplicative selectivity curve. Ordinary terms provide intercepts.
#'
#' F smooths support `process = "off"`, `"iid"` or `"ar1"`. Smooth means with
#' residual variation are experimental: both may describe the same pattern.
#' Inspect [check_tam()], uncertainty, and recovery simulations. Intrinsic
#' RW/correlated-RW F residuals are not supported with smooth means.
#'
#' Bases are constructed from historical covariates only and rebuilt in every
#' fold. New temporal smooth contributions hold the terminal training year in
#' projections; existing [rw()]/[ar1()] q-effect forecasts retain their process
#' definitions. Warm starts are retained only for compatible bases. Simulation
#' draws each penalized coefficient once per realization and applies the same
#' function throughout its grid. Conditional draws retain fitted coefficients.
#'
#' The conversion follows mgcv's `smoothCon()`/`smooth2random(type = 2)` and
#' the approach in sdmTMB. Prediction transforms are adapted from Simon Wood's
#' mgcv `smooth2random` example, also used by sdmTMB's `s2rPred()`.
#' Validation pins mgcv 1.9.4 and sdmTMB revision
#' `34740199c51fbebe86a6de7671f33f590db166bf`.
#' @references Wood, S. N. (2017). *Generalized Additive Models: An Introduction
#'   with R*, second edition. Chapman and Hall/CRC.
#'   <https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html>
#'   sdmTMB contributors, `R/smoothers.R` (GPL-3):
#'   <https://github.com/sdmTMB/sdmTMB/blob/34740199c51fbebe86a6de7671f33f590db166bf/R/smoothers.R>
#' @examples
#' ~ fleet + s(age, k = 5, by = fleet)
#' ~ fleet + t2(age, year, k = c(4, 5), bs = c("cr", "tp"), by = fleet)
#' ~ 0 + for_fleet("Trawl", ~ factor(year) + logistic(age)) +
#'   for_fleet("Gillnet", ~ factor(year) + s(age, k = 5))
#' ~ 0 + for_survey("RV", ~ s(age, k = 5)) +
#'   for_survey("Sentinel", ~ logistic(age))
#' @seealso [fit_tam()], [for_fleet()], [formula_effects], [double_logistic()]
#' @name smooth_formulas
NULL

.smooth_spec <- function(call, environment) {
  type <- .formula_call_name(call)
  fun <- if (type == "s") mgcv::s else mgcv::t2
  call[[1L]] <- as.name(type)
  args <- tryCatch(as.list(match.call(fun, call, expand.dots = FALSE))[-1L],
    error = function(e) cli::cli_abort("Invalid {type}() arguments: {conditionMessage(e)}"))
  unsupported <- setdiff(names(args), c("...", "k", "bs", "by"))
  if (length(unsupported)) cli::cli_abort("Unsupported {type}() options: {paste(unsupported, collapse = ', ')}. Supported options are k, bs and categorical by.")
  coordinates <- as.list(args$...)
  if (!length(coordinates) || any(nzchar(names(coordinates))) ||
      !all(vapply(coordinates, is.symbol, logical(1)))) {
    cli::cli_abort("{type}() predictors must be unnamed, bare column names.")
  }
  variables <- vapply(coordinates, as.character, character(1))
  if (type == "s" && length(variables) != 1L) cli::cli_abort("Only univariate s() is supported; use t2(age, year) for age-year smooths.")
  if (type == "t2" && (length(variables) != 2L || !setequal(variables,c("age","year")))) {
    cli::cli_abort("t2() currently supports age and year only.")
  }
  evaluate <- function(x) tryCatch(eval(x, environment),
    error = function(e) cli::cli_abort("Cannot evaluate {type}() option: {conditionMessage(e)}"))
  bs <- if (is.null(args$bs)) if (type == "s") "tp" else "cr" else evaluate(args$bs)
  if (!is.character(bs) || anyNA(bs) || any(!bs %in% c("tp","cr")) ||
      !length(bs) || !length(bs) %in% c(1L,length(variables))) {
    cli::cli_abort("Supported smooth bases are 'tp' and 'cr', one per predictor.")
  }
  if (!is.null(args$k)) {
    k <- evaluate(args$k)
    if (!is.numeric(k) || any(!is.finite(k)) || any(k < 3 | k != floor(k)) ||
        !length(k) %in% c(1L,length(variables))) cli::cli_abort("Smooth k must contain integer basis sizes of at least three, one per predictor.")
  }
  by <- NULL
  if (!is.null(args$by) && !identical(args$by,quote(NULL)) && !identical(args$by,quote(NA))) {
    if (!is.symbol(args$by)) cli::cli_abort("Smooth by must name a categorical input column.")
    by <- as.character(args$by)
  }
  if (is.null(by) && "by" %in% names(args)) {
    parts <- as.list(call)
    call <- as.call(parts[names(parts) != "by"])
  }
  eval_environment <- list2env(list(s=mgcv::s,t2=mgcv::t2),parent=environment)
  object <- tryCatch(eval(call,eval_environment),
    error = function(e) cli::cli_abort("Cannot construct {type}(): {conditionMessage(e)}"))
  list(type="smooth", variable=paste(variables,collapse=","), variables=variables,
    by=by, call=call, object=object, label=object$label)
}

# Adapted from mgcv's smooth2random example (Simon N. Wood, GPL >= 2),
# also used by sdmTMB::s2rPred at the revision cited above (GPL-3).
.smooth_prediction <- function(smooth, random, data) {
  X <- mgcv::PredictMat(smooth,data)
  if (!is.null(random$trans.U)) X <- X %*% random$trans.U
  X <- t(t(X)*random$trans.D)
  X[,random$rind] <- X[,random$pen.ind != 0]
  index <- random$pen.ind
  index[random$rind] <- index[index > 0]
  list(fixed=X[,random$pen.ind == 0,drop=FALSE],
    random=lapply(seq_along(random$rand),function(i)X[,index == i,drop=FALSE]))
}

.compile_smooths <- function(specs, data) {
  fixed <- matrix(numeric(),nrow(data),0)
  terms <- list()
  historical <- if (is.null(data$is_proj)) rep(TRUE,nrow(data)) else !data$is_proj
  for (spec in specs) {
    columns <- c(spec$variables,spec$by)
    missing <- setdiff(columns,names(data))
    if (length(missing)) cli::cli_abort("{spec$label} is missing input columns: {paste(missing,collapse=', ')}.")
    for (variable in spec$variables) if (!is.numeric(data[[variable]]) || any(!is.finite(data[[variable]]))) {
      cli::cli_abort("{spec$label} needs complete finite numeric {variable} values, including rows with missing observations.")
    }
    d <- as.data.frame(data)
    if (!is.null(spec$by)) {
      x <- d[[spec$by]]
      if (!(is.factor(x) || is.character(x) || is.logical(x)) || anyNA(x)) cli::cli_abort("Smooth by must be a complete categorical column; numeric by is not supported.")
      if (!is.factor(x)) x <- factor(x)
      d[[spec$by]] <- x
    }
    training <- .drop_unused_levels(d[historical,,drop=FALSE])
    if ("year" %in% spec$variables && any(!historical)) {
      d$year[!historical] <- max(training$year)
    }
    bases <- tryCatch(mgcv::smoothCon(spec$object,training,absorb.cons=TRUE,modCon=3,
      diagonal.penalty=FALSE), error=function(e)cli::cli_abort(
        "Cannot build {spec$label}: {conditionMessage(e)}. Check distinct training values and reduce k if needed.",parent=e))
    for (basis in bases) {
      random <- mgcv::smooth2random(basis,names(training),type=2)
      design <- .smooth_prediction(basis,random,d)
      id <- paste0("smooth_",basis$label)
      labels <- paste0(id,":null",seq_len(ncol(design$fixed)))
      colnames(design$fixed) <- labels
      fixed <- cbind(fixed,design$fixed)
      rows <- if (is.null(spec$by)) seq_len(nrow(d)) else which(as.character(d[[spec$by]]) == basis$by.level)
      support <- rows[historical[rows]]
      if (!is.null(d$obs)) support <- support[is.finite(d$obs[support]) & d$obs[support] > 0]
      if (!length(support)) cli::cli_abort("{basis$label} has no informative historical observations.")
      reference <- unique(training[columns])
      reference <- reference[do.call(order,reference),,drop=FALSE]
      rownames(reference) <- NULL
      signature <- .smooth_prediction(basis,random,reference)
      for (j in seq_along(design$random)) {
        key <- paste0(id,":penalty",j)
        Z <- design$random[[j]]
        states <- paste0("basis",seq_len(ncol(Z)))
        colnames(Z) <- states
        term <- list(type="smooth",id=key,variable=spec$variable,variables=spec$variables,
          by=spec$by,sd=NULL,phi=NULL,parameter=paste0("eta_q_",key),
          sd_parameter=paste0("log_sd_q_",key),phi_parameter=NULL,
          multiplier=rep(1,nrow(d)),design=Z,
          groups=list(list(group=if (is.null(spec$by)) "all" else as.character(basis$by.level),
            rows=rows,levels=states,states=states,n_historical=length(states),
            index=seq_along(states),is_proj=rep(FALSE,length(states)))),
          smooth_id=id,smooth_fixed=design$fixed,smooth=basis,random_conversion=random,
          penalty=j,basis_signature=list(grid=reference,design=signature))
        terms[[key]] <- term
      }
    }
  }
  list(fixed=fixed,terms=terms)
}

# A smooth2random basis can change eigenvector signs (or rotate an eigenspace)
# without changing its function or Gaussian penalty. Transfer coefficients only
# when both the function and, for random coefficients, the penalty are preserved.
.smooth_design_transform <- function(old, new, orthogonal = TRUE) {
  tolerance <- 1e-10
  if (!is.matrix(old) || !is.matrix(new) || !identical(dim(old), dim(new)) ||
      any(!is.finite(old)) || any(!is.finite(new))) return(NULL)
  n <- ncol(old)
  if (!n) return(matrix(numeric(), 0L, 0L))
  if (identical(unname(old), unname(new))) return(diag(n))
  transform <- tryCatch(qr.solve(new, old, tol = tolerance),
    error = function(e) NULL)
  if (is.null(transform) || any(!is.finite(transform))) return(NULL)
  scale <- max(abs(old))
  if (!scale || max(abs(new %*% transform - old)) > tolerance * scale) return(NULL)
  if (orthogonal && max(abs(crossprod(transform) - diag(n))) > tolerance) return(NULL)
  if (!orthogonal && qr(transform, tol = tolerance)$rank != n) return(NULL)
  transform
}

.smooth_basis_transform <- function(old, new) {
  fields <- c("type", "smooth_id", "variables", "by", "scope")
  if (!identical(old$type, "smooth") || !identical(new$type, "smooth") ||
      !identical(old[fields], new[fields]) ||
      !identical(class(old$smooth), class(new$smooth)) ||
      !identical(old$smooth$bs.dim, new$smooth$bs.dim) ||
      !isTRUE(all.equal(old$basis_signature$grid, new$basis_signature$grid,
        tolerance = 0))) return(NULL)
  before <- old$basis_signature$design
  after <- new$basis_signature$design
  if (length(before$random) != length(after$random)) return(NULL)
  fixed <- .smooth_design_transform(before$fixed, after$fixed, orthogonal = FALSE)
  random <- Map(.smooth_design_transform, before$random, after$random)
  if (is.null(fixed) || any(vapply(random, is.null, logical(1)))) return(NULL)
  list(fixed = fixed, random = random)
}

.smooth_basis_compatible <- function(old, new) {
  !is.null(.smooth_basis_transform(old, new))
}

.smooth_warm_starts <- function(start, old, new) {
  designs <- list(F = old$F_modmat, q = old$q_modmat)
  for (component in c("q", "F")) {
    field <- if (component == "q") "q_terms" else "F_terms"
    coefficient <- if (component == "F") "log_mu_f" else
      if (new$index_settings$q_link == "logit") "logit_q" else "log_q"
    terms <- Filter(function(x) x$type == "smooth", new[[field]])
    ids <- unique(vapply(terms, `[[`, character(1), "smooth_id"))
    for (id in ids) {
      pieces <- Filter(function(x) identical(x$smooth_id, id), terms)
      previous <- Filter(function(x) identical(x$type, "smooth") &&
        identical(x$smooth_id, id), old[[field]])
      transform <- if (identical(names(previous), names(pieces)))
        .smooth_basis_transform(previous[[1L]], pieces[[1L]]) else NULL
      columns <- colnames(pieces[[1L]]$smooth_fixed)
      if (is.null(transform)) {
        parameters <- unlist(lapply(pieces, function(x) c(x$parameter, x$sd_parameter)))
        start[parameters] <- NULL
        if (!is.null(start[[coefficient]])) {
          start[[coefficient]] <- start[[coefficient]][!names(start[[coefficient]]) %in% columns]
        }
        next
      }
      for (piece in pieces) {
        parameter <- piece$parameter
        if (!is.null(start[[parameter]])) {
          value <- drop(transform$random[[piece$penalty]] %*% start[[parameter]])
          names(value) <- names(start[[parameter]])
          start[[parameter]] <- value
        }
      }
      if (length(columns)) {
        if (!is.null(start[[coefficient]]) && all(columns %in% names(start[[coefficient]]))) {
          start[[coefficient]][columns] <- drop(transform$fixed %*% start[[coefficient]][columns])
        }
        # Express the old design in the new coordinates before checking the
        # ordinary fixed-effect design, including the smooth null space.
        designs[[component]][, columns] <- designs[[component]][, columns, drop = FALSE] %*%
          solve(transform$fixed)
      }
    }
  }
  list(parameters = start, designs = designs)
}

.smooth_function_values <- function(par, dat, component, effects) {
  terms <- Filter(function(x)x$type == "smooth",dat[[if(component == "q") "q_terms" else "F_terms"]])
  n <- if (component == "q") nrow(dat$obs$index) else nrow(dat$F_data)
  ids <- unique(vapply(terms,`[[`,character(1),"smooth_id"))
  out <- matrix(0,n,length(ids),dimnames=list(row=seq_len(n),smooth=ids))
  if (!length(ids)) return(out)
  p <- par
  p[names(effects$parameters)] <- effects$parameters
  coefficient <- if (component == "F") "log_mu_f" else
    if (dat$index_settings$q_link == "logit") "logit_q" else "log_q"
  "[<-" <- RTMB::ADoverload("[<-")
  for (id in ids) {
    pieces <- Filter(function(x)x$smooth_id == id,terms)
    X <- pieces[[1L]]$smooth_fixed
    value <- drop(X %*% p[[coefficient]][colnames(X)])
    for (term in pieces) value <- value + drop(term$design %*% p[[term$parameter]])
    out[,id] <- value
  }
  out
}

.tidy_smooth_functions <- function(fit, interval, component) {
  terms <- Filter(function(x)x$type == "smooth",fit$dat[[if(component == "q")"q_terms" else "F_terms"]])
  if (!length(terms)) return(list())
  field <- if (component == "q") "eta_q_smooth" else "eta_mu_F_smooth"
  estimate <- fit$rep[[field]]
  if (is.null(estimate)) return(list())
  error <- if (is.null(fit[["sdrep"]])) NULL else as.list(fit$sdrep,"Std. Error",report=TRUE)[[field]]
  data <- if (component == "q") fit$dat$obs$index else .mean_formula_data(fit$dat,"F")
  ids <- unique(vapply(terms,`[[`,character(1),"smooth_id"))
  z <- stats::qnorm(.5+interval/2)
  out <- lapply(seq_along(ids),function(j) {
    term <- Filter(function(x)x$smooth_id == ids[j],terms)[[1L]]
    rows <- term$groups[[1L]]$rows
    columns <- unique(c(intersect(c("year","age","fleet","survey","is_proj"),names(data)),term$variables,term$by))
    tab <- data[rows,columns,drop=FALSE]
    tab$term <- ids[j]; tab$component <- component; tab$process <- "smooth"
    tab$variable <- term$variable; tab$group <- term$groups[[1L]]$group
    tab$by <- if (is.null(term$by)) NA_character_ else term$by
    tab$est <- as.numeric(estimate[rows,j])
    tab$se <- if (is.null(error)) NA_real_ else as.numeric(error[rows,j])
    tab$lwr <- tab$est-z*tab$se; tab$upr <- tab$est+z*tab$se
    tab$se_scale <- "reported"
    tab$unit <- if (component == "F") "log F contribution" else
      if (fit$dat$index_settings$q_link == "logit")"logit q contribution" else "log q contribution"
    tab
  })
  link <- if (component == "q") fit$dat$index_settings$q_link else "log"
  stats::setNames(out,paste0(component,":",link,":",ids))
}

.smooth_design_advisories <- function(dat) {
  out <- data.frame(issue=character(),detail=character())
  F <- Filter(function(x)x$type == "smooth",dat$F_terms)
  if (length(F) && dat$F_settings$process %in% c("iid","ar1")) out <- rbind(out,data.frame(
    issue="smooth_F_residuals",detail=paste(
      "Smooth F means with IID/AR1 residuals are experimental. Both components may explain the same age/year patterns;",
      "inspect smoothing and residual SD uncertainty, their correlation, and recovery simulations. No component is simplified automatically.")))
  for (component in c("q","F")) for (term in Filter(function(x)x$type == "smooth",
      dat[[if(component == "q") "q_terms" else "F_terms"]])) {
    rows <- term$groups[[1L]]$rows
    d <- if (component == "q") dat$obs$index else .mean_formula_data(dat,"F")
    eligible <- !d$is_proj
    if (component == "q") eligible <- eligible & is.finite(d$obs) & d$obs > 0 else
      eligible <- eligible & .f_state_rows(dat)
    rows <- rows[eligible[rows]]
    cells <- nrow(unique(d[rows,term$variables,drop=FALSE]))
    if (cells < 5L) out <- rbind(out,data.frame(issue="smooth_support",
      detail=paste0(term$id,": only ",cells," distinct historical predictor cells support the smoothing SD. Repeated ages are not independent function coordinates; inspect uncertainty.")))
  }
  out
}

.smooth_fit_advisories <- function(fit) {
  out <- data.frame(issue=character(),detail=character())
  sdr <- fit[["sdrep"]]
  if (is.null(sdr) || !isTRUE(sdr$pdHess) || is.null(sdr$cov.fixed)) return(out)
  terms <- Filter(function(x)x$type == "smooth",fit$dat$F_terms)
  if (!length(terms) || !fit$dat$F_settings$process %in% c("iid","ar1")) return(out)
  smooth <- unlist(lapply(terms,`[[`,"sd_parameter"),use.names=FALSE)
  residual <- fit$dat$process_sd$F$parameter
  i <- which(names(sdr$par.fixed) %in% smooth)
  j <- which(names(sdr$par.fixed) %in% residual)
  variance <- diag(sdr$cov.fixed)
  if (!length(i) || !length(j) || any(!is.finite(variance[c(i,j)]) | variance[c(i,j)] <= 0)) return(out)
  correlation <- sdr$cov.fixed[i,j,drop=FALSE]/sqrt(outer(variance[i],variance[j]))
  if (any(abs(correlation) > .9,na.rm=TRUE)) out <- rbind(out,data.frame(
    issue="smooth_F_separation",detail="Smoothing and F residual SD estimates have absolute correlation above 0.9. Their separate contributions may be weakly identified; compare supported simpler hypotheses and recovery simulations."))
  out
}
