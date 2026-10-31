#' Apply a formula within one fleet or survey
#'
#' Give a selected fishery its own F relationship, or a selected survey its own
#' catchability relationship, while keeping one formula. Ordinary terms outside
#' the marker can share effects across groups.
#' @param fleet,survey One exact name from the catch or index table.
#' @param formula A one-sided formula evaluated on rows of that fleet/survey.
#' @details Use `for_fleet()` in `F_settings$mu_form` and `for_survey()` in
#' `index_settings$q_form`. Each inner formula has its own intercept. Use an
#' outer `~ 0 + ...` when the inner formulas supply all group baselines.
#' Scoped contributions are zero outside their selected rows; this means no
#' contribution to log F or the q predictor, rather than zero fishing mortality.
#' Joint design checks reject redundant intercepts and age relationships.
#' @return Formula markers; direct calls raise an error.
#' @examples
#' ~ 0 + for_fleet("Trawl", ~ factor(year) + logistic(age)) +
#'   for_fleet("Gillnet", ~ factor(year) + factor(age))
#' ~ 0 + for_survey("RV", ~ logistic(age)) +
#'   for_survey("Sentinel", ~ factor(age))
#' @seealso [fit_tam()], [catchability_curves], [formula_effects]
#' @export
for_fleet <- function(fleet, formula) {
  cli::cli_abort("for_fleet() is a formula marker for F_settings$mu_form.")
}

#' @rdname for_fleet
#' @export
for_survey <- function(survey, formula) {
  cli::cli_abort("for_survey() is a formula marker for index_settings$q_form.")
}

.contains_scope <- function(x) is.call(x) && (
  .formula_call_name(x) %in% c("for_fleet", "for_survey") ||
    any(vapply(as.list(x)[-1L], .contains_scope, logical(1))))

.formula_has_call <- function(formula, names) {
  walk <- function(x, env) {
    if (!is.call(x)) return(FALSE)
    name <- .formula_call_name(x)
    if (name %in% names) return(TRUE)
    if (name %in% c("for_fleet","for_survey")) {
      spec <- .scope_spec(x,env)
      return(walk(spec$formula,environment(spec$formula)))
    }
    any(vapply(as.list(x)[-1L],walk,logical(1),env=env))
  }
  walk(formula,if(inherits(formula,"formula")) environment(formula) else parent.frame())
}

# Input columns, resolving scoped formula objects in their own environments.
.formula_input_variables <- function(formula) {
  walk <- function(x, env) {
    if (is.symbol(x)) return(as.character(x))
    if (!is.call(x)) return(character())
    name <- .formula_call_name(x)
    if (name %in% c("for_fleet", "for_survey")) {
      spec <- .scope_spec(x, env)
      return(walk(spec$formula[[2L]], environment(spec$formula)))
    }
    parts <- as.list(x)[-1L]
    if (name %in% c("iid", "rw", "ar1", "warped_normal", "double_logistic")) {
      parts <- as.list(match.call(get(name, mode="function"), x))[-1L]
      parts <- parts[!names(parts) %in% c("sd", "phi", "b")]
    }
    unique(unlist(lapply(parts, walk, env=env),use.names=FALSE))
  }
  if (is.null(formula)) character() else walk(formula[[2L]], environment(formula))
}

.scope_spec <- function(x, environment) {
  type <- sub("^for_", "", .formula_call_name(x))
  marker <- if (type == "fleet") for_fleet else for_survey
  args <- tryCatch(as.list(match.call(marker, x))[-1L],
    error = function(e) cli::cli_abort("Invalid for_{type}() arguments: {conditionMessage(e)}"))
  if (is.null(args[[type]]) || is.null(args$formula)) cli::cli_abort("for_{type}() needs a name and a one-sided formula.")
  name <- eval(args[[type]], environment)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    cli::cli_abort("for_{type}() needs one nonempty {type} name.")
  }
  formula <- eval(args$formula, environment)
  if (!inherits(formula, "formula") || length(formula) != 2L) cli::cli_abort("for_{type}() needs a one-sided formula.")
  if (.contains_scope(formula)) cli::cli_abort("Do not nest fleet/survey scopes.")
  list(type = type, name = name, formula = formula,
    id = paste0(type, "[", name, "]"))
}

.split_scopes <- function(formula, data, allowed) {
  scopes <- list()
  walk <- function(x) {
    if (!.contains_scope(x)) return(x)
    type <- .formula_call_name(x)
    if (type %in% c("for_fleet", "for_survey")) {
      spec <- .scope_spec(x, environment(formula))
      if (is.null(allowed) || spec$type != allowed) {
        cli::cli_abort("for_{spec$type}() is only supported in {if (spec$type == 'fleet') 'F_settings$mu_form' else 'index_settings$q_form'}.")
      }
      if (!spec$type %in% names(data)) cli::cli_abort("for_{spec$type}() requires {spec$type} labels in its input table.")
      available <- unique(as.character(data[[spec$type]]))
      if (!spec$name %in% available) cli::cli_abort("Unknown {spec$type} {.val {spec$name}}. Available names: {paste(available, collapse = ', ')}.")
      spec$rows <- which(as.character(data[[spec$type]]) == spec$name)
      scopes[[length(scopes) + 1L]] <<- spec
      return(NULL)
    }
    if (type == "(") return(walk(x[[2L]]))
    if (type == "+") {
      parts <- Filter(Negate(is.null), lapply(as.list(x)[-1L], walk))
      if (!length(parts)) return(NULL)
      if (length(parts) == 1L) return(parts[[1L]])
      return(Reduce(function(a,b) call("+",a,b),parts))
    }
    if (type == "-" && length(x) == 3L && !.contains_scope(x[[3L]])) {
      left <- walk(x[[2L]])
      return(call("-",if(is.null(left)) 1 else left,x[[3L]]))
    }
    cli::cli_abort("Fleet/survey scopes must be additive terms, without interactions or transformations around a scope.")
  }
  ordinary <- formula
  rhs <- walk(formula[[2L]])
  ordinary[[2L]] <- if (is.null(rhs)) 1 else rhs
  list(formula = ordinary, scopes = scopes)
}

.scoped_q_design <- function(formula, data, fixed_design, allowed) {
  split <- .split_scopes(formula, data, allowed)
  global <- .parse_q_formula(split$formula, data, fixed_design)
  out <- global
  out$q_scope_components <- list()
  n <- nrow(data)
  padding <- function(x, rows, id) {
    z <- matrix(0, n, ncol(x))
    z[rows, ] <- x
    colnames(z) <- if(ncol(x)) paste0(id, ":", colnames(x)) else character()
    z
  }
  has_curve <- function(design) !is.null(design$q_mono_modmat) ||
    any(vapply(design$q_terms, .is_selectivity_curve, logical(1)))
  curve_rows <- if (has_curve(global)) list(seq_len(n)) else list()
  for (spec in split$scopes) {
    rows <- spec$rows
    piece <- tryCatch(.parse_q_formula(spec$formula, data[rows, , drop = FALSE], fixed_design),
      error=function(e) cli::cli_abort("Cannot compile {spec$id}: {conditionMessage(e)}",parent=e))
    if (has_curve(piece)) {
      if (any(vapply(curve_rows, function(x) length(intersect(x,rows)) > 0L, logical(1)))) {
        cli::cli_abort("Overlapping structural curves in {spec$id}. Use one age/size curve on each applicable row.")
      }
      curve_rows[[length(curve_rows) + 1L]] <- rows
    }
    out$q_modmat <- cbind(out$q_modmat, padding(piece$q_modmat, rows, spec$id))
    if (!is.null(piece$q_mono_modmat)) {
      out$q_mono_modmat <- cbind(out$q_mono_modmat, padding(piece$q_mono_modmat, rows, spec$id))
      steps <- piece$q_mono_steps
      steps$coef <- paste0(spec$id, ":", steps$coef)
      steps$scope <- spec$name
      out$q_mono_steps <- .bind_rows(list(out$q_mono_steps, steps))
    }
    for (term in piece$q_terms) {
      term$id <- paste0(spec$id, ":", term$id)
      term$scope <- spec[c("type","name","id")]
      if (.is_selectivity_curve(term)) term <- .curve_parameter_names(term)
      prefixes <- c(parameter="eta_q_",sd_parameter="log_sd_q_",phi_parameter="logit_phi_q_",
        a50_parameter="q_a50_",slope_parameter="log_q_slope_")
      for (field in c("parameter","sd_parameter","phi_parameter")) {
        if (!is.null(term[[field]])) {
          prefix <- prefixes[[field]]
          term[[field]] <- paste0(prefix,spec$id,":",substring(term[[field]],nchar(prefix)+1L))
        }
      }
      for (i in seq_along(term$groups)) term$groups[[i]]$rows <- rows[term$groups[[i]]$rows]
      multiplier <- numeric(n); multiplier[rows] <- term$multiplier; term$multiplier <- multiplier
      if (term$id %in% names(out$q_terms)) cli::cli_abort("Duplicate scoped term {term$id}.")
      out$q_terms[[term$id]] <- term
    }
    if (!is.null(piece$q_power)) {
      if (is.null(out$q_power)) out$q_power <- list(by=NULL,levels=character(),index=integer(n))
      active <- piece$q_power$index > 0L
      if (any(out$q_power$index[rows[active]] > 0L)) cli::cli_abort("Overlapping power(N) terms in {spec$id}.")
      offset <- length(out$q_power$levels)
      out$q_power$levels <- c(out$q_power$levels,paste0(spec$id, ":",piece$q_power$levels))
      out$q_power$index[rows[active]] <- piece$q_power$index[active] + offset
    }
    out$q_scope_components[[length(out$q_scope_components) + 1L]] <- spec
  }
  observed <- if (is.null(data$is_proj)) rep(TRUE,n) else !data$is_proj
  if (allowed == "survey") observed <- observed & is.finite(data$obs) & data$obs > 0
  fixed <- out$q_modmat[observed,,drop=FALSE]
  if (!is.null(out$q_mono_modmat)) fixed <- cbind(fixed,out$q_mono_modmat[observed,,drop=FALSE])
  if (ncol(fixed) && qr(fixed)$rank < ncol(fixed)) {
    cli::cli_abort("Scoped formula has redundant or unsupported fixed terms. Remove overlapping intercepts or age effects within the applicable fleet/survey.")
  }
  out
}

.check_f_curve_design <- function(design, data, terms) {
  specifications <- lapply(Filter(.is_selectivity_curve,terms),function(x) {
    list(variable=x$variable, by=x$by, rows=unique(unlist(lapply(x$groups,`[[`,"rows"))))
  })
  mono <- design$q_mono_steps
  if (!is.null(mono)) {
    if (!"scope" %in% names(mono)) mono$scope <- NA_character_
    for (i in which(!duplicated(mono[c("variable","by","scope")]))) {
      rows <- seq_len(nrow(data))
      if (!is.na(mono$scope[i])) {
        scope <- Filter(function(x) x$name == mono$scope[i],design$q_scope_components)[[1L]]
        rows <- scope$rows
      }
      specifications[[length(specifications)+1L]] <- list(variable=mono$variable[i],by=mono$by[i],rows=rows)
    }
  }
  historical <- !data$is_proj
  fixed <- design$q_modmat[historical,,drop=FALSE]
  for (spec in specifications) {
    d <- data[spec$rows,,drop=FALSE]; x <- d[[spec$variable]]
    group <- interaction(d$age,if("fleet" %in% names(d)) d$fleet else 1,drop=TRUE)
    if (!is.numeric(x) || any(vapply(split(x,group),function(z) length(unique(z)) != 1L,logical(1)))) {
      cli::cli_abort("F mean curves must use numeric age or an age-derived catch covariate constant across years at each age within its scope.")
    }
    if (!is.null(spec$by) && !is.na(spec$by) && is.numeric(d[[spec$by]])) cli::cli_abort("F mean curves require categorical by groups.")
    indicator <- as.numeric(seq_len(nrow(data)) %in% spec$rows)[historical]
    if (!ncol(fixed) || qr(cbind(fixed,indicator))$rank > qr(fixed)$rank) {
      cli::cli_abort("F mean curves require an ordinary baseline within their applicable fleet; retain an inner intercept or supply group baselines.")
    }
  }
  invisible(NULL)
}
