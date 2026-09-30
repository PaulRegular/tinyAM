#' Translate SAM assumptions into tinyAM fitting settings
#'
#' Generate a transparent approximation using existing tinyAM options. Review
#' [sam_to_tam_audit()] before fitting; translating settings does not reproduce
#' the SAM likelihood and does not guarantee numerical convergence.
#'
#' @param sam_fit A fitted SAM object.
#' @param overrides Named list of `fit_tam()` arguments. Nested settings are
#'   merged by name; explicit `NULL` values replace defaults. `obs` is supplied
#'   separately and cannot be overridden here.
#' @return A named list suitable for
#'   `do.call(fit_tam, c(list(obs = sam_to_tam_obs(sam_fit)), settings))`.
#' @details
#' The baseline uses IID survival errors, free initial abundance, independent
#' F random-walk increments, supplied input M, and original weight/maturity.
#' q and observation-SD blocks retain SAM's keys. Observation errors are
#' independent lognormal; catch scaling, biological process models, correlated
#' innovations and additional variance groups are not reproduced.
#' Missing biological values are never filled. Choose a complete input period
#' explicitly and subset the translated observations before calling `fit_tam()`.
#' No fitted SAM estimates constrain the resulting tinyAM fit.
#' @export
sam_to_tam_settings <- function(sam_fit, overrides = list()) {
  x <- .sam_source(sam_fit)
  obs <- sam_to_tam_obs(sam_fit)
  if (anyNA(obs$index$q_key) || anyNA(obs$index$sd_key) || anyNA(obs$catch$sd_key) ||
      any(obs$index$sd_key < 0) || any(obs$catch$sd_key < 0)) {
    cli::cli_abort("SAM active observation cells require complete q and observation-SD keys.")
  }
  q_form <- if (any(obs$index$q_key == -1)) {
    cols <- grep("^q_key_", names(obs$index), value = TRUE)
    if (length(cols)) stats::reformulate(cols, intercept = FALSE) else ~ 0
  } else if (nlevels(obs$index$q_block) == 1L) ~ 1 else ~ 0 + q_block
  catch_sd <- if (nlevels(obs$catch$sd_block) == 1L) ~ 1 else ~ 0 + sd_block
  index_sd <- if (nlevels(obs$index$sd_block) == 1L) ~ 1 else ~ 0 + sd_block
  if (length(x$conf$fbarRange) != 2L || anyNA(x$conf$fbarRange) ||
      any(!x$conf$fbarRange %in% seq.int(x$conf$minAge, x$conf$maxAge))) {
    cli::cli_abort("SAM fbarRange must identify two modeled age bounds.")
  }
  settings <- list(years = x$years, ages = seq.int(x$conf$minAge, x$conf$maxAge),
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "rw", mu_form = NULL,
                      mean_ages = seq.int(min(x$conf$fbarRange), max(x$conf$fbarRange))),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption,
                      age_breaks = NULL, first_dev_year = NULL),
    catch_settings = list(sd_form = catch_sd, sd_supplied = NULL, fill_missing = FALSE),
    index_settings = list(q_form = q_form, sd_form = index_sd,
                          sd_supplied = NULL, fill_missing = FALSE),
    proj_settings = NULL)
  overrides <- .validate_named_list(overrides, "overrides", allow_empty = TRUE)
  allowed <- c(setdiff(names(formals(make_dat)), "obs"),
               setdiff(names(formals(fit_tam)), c("obs", "...")))
  unknown <- setdiff(names(overrides), allowed)
  if (length(unknown)) cli::cli_abort("Unknown settings override(s): {paste(unknown, collapse = ', ')}.")
  utils::modifyList(settings, overrides, keep.null = TRUE)
}

.sam_design_matches <- function(form, data, key) {
  if (is.null(form) || anyNA(key)) return(FALSE)
  if (any(grepl("mono\\s*\\(", deparse(form)))) return(FALSE)
  design <- tryCatch(stats::model.matrix(form, data), error = function(e) NULL)
  if (is.null(design) || nrow(design) != nrow(data) || any(!is.finite(design))) return(FALSE)
  levels <- sort(unique(key[key >= 0]))
  expected <- matrix(0, length(key), length(levels))
  for (i in seq_along(levels)) expected[, i] <- as.numeric(key == levels[i])
  rank <- function(m) if (!ncol(m)) 0L else qr(m)$rank
  rank(design) == rank(expected) && rank(cbind(design, expected)) == rank(expected)
}

.sam_audit_settings <- function(out, sam_fit, settings) {
  obs <- sam_to_tam_obs(sam_fit)
  describe <- function(x) paste(deparse(x, width.cutoff = 80L), collapse = " ")
  out$tam_setting <- "See mapping"
  set <- function(field, setting, exact = TRUE, note = "") {
    i <- match(field, out$sam_setting)
    if (is.na(i)) return(invisible(NULL))
    out$tam_setting[i] <<- describe(setting)
    if (!exact && out$tam_status[i] != "not_checked") out$tam_status[i] <<- "unsupported"
    if (nzchar(note)) out$notes[i] <<- paste(out$notes[i], note)
  }
  set("minAge/maxAge", settings$ages, identical(as.integer(settings$ages), seq.int(sam_fit$conf$minAge, sam_fit$conf$maxAge)))
  set("corFlag", settings$F_settings, identical(settings$F_settings$process, "rw") &&
        is.null(settings$F_settings$mu_form),
      "The applied F model is shown in tam_setting; no correlated innovations are introduced.")
  set("keyVarF", settings$F_settings$process, identical(settings$F_settings$process, "rw"))
  set("keyVarLogN", settings$N_settings, identical(settings$N_settings$process, "iid"))
  set("initState", settings$N_settings$init, identical(settings$N_settings$init, "free"))
  set("keyLogFpar", settings$index_settings$q_form,
      .sam_design_matches(settings$index_settings$q_form, obs$index, obs$index$q_key))
  set("keyVarObs", list(catch = settings$catch_settings$sd_form, index = settings$index_settings$sd_form),
      .sam_design_matches(settings$catch_settings$sd_form, obs$catch, obs$catch$sd_key) &&
        .sam_design_matches(settings$index_settings$sd_form, obs$index, obs$index$sd_key) &&
        is.null(settings$catch_settings$sd_supplied) && is.null(settings$index_settings$sd_supplied))
  supplied <- tryCatch(.log_supplied(settings$M_settings$mu_supplied, obs$weight, "M mu_supplied"), error = function(e) NULL)
  set("nm.dat", settings$M_settings,
      identical(settings$M_settings$process, "off") && is.null(settings$M_settings$mu_form) &&
        !is.null(supplied) && isTRUE(all.equal(exp(supplied), obs$weight$M_assumption, check.attributes = FALSE)))
  set("fbarRange", settings$F_settings$mean_ages)
  selected <- lapply(obs, function(d) d[d$year %in% settings$years & d$age %in% settings$ages, ])
  complete <- length(settings$years) > 0 && length(settings$ages) > 0 &&
    all(settings$years %in% sam_fit$data$years) &&
    all(settings$ages %in% seq.int(sam_fit$conf$minAge, sam_fit$conf$maxAge)) &&
    nrow(selected$weight) == length(settings$years) * length(settings$ages) &&
    all(is.finite(selected$weight$obs)) && all(is.finite(selected$maturity$obs)) &&
    all(is.finite(selected$weight$M_assumption) & selected$weight$M_assumption > 0)
  extra <- data.frame(component = "input period", sam_setting = "selected years",
    sam_value = paste(range(sam_fit$data$years), collapse = ":"),
    tam_status = if (!complete) "unsupported" else if (identical(as.numeric(settings$years), as.numeric(sam_fit$data$years))) "supported" else "partially_supported",
    tam_mapping = "Original biological inputs; no automatic imputation",
    notes = if (complete) "Comparison must use the common fitted period." else "Selected original biological inputs contain missing values; explicitly choose a complete period before fitting.",
    tam_setting = paste(range(settings$years), collapse = ":"))
  rbind(out, extra)
}
