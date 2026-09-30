.sam_estimates <- function(log_est, log_se = NULL, interval = 0.95) {
  if (is.null(log_se)) log_se <- rep(NA_real_, length(log_est))
  z <- stats::qnorm(0.5 + interval / 2)
  data.frame(est = exp(as.vector(log_est)), lwr = exp(as.vector(log_est) - z * as.vector(log_se)),
    upr = exp(as.vector(log_est) + z * as.vector(log_se)), se = as.vector(log_se), se_scale = "log")
}

.sam_comparison_tables <- function(sam_fit, interval) {
  if (!requireNamespace("stockassessment", quietly = TRUE)) cli::cli_abort("Install stockassessment to extract SAM comparison tables.")
  x <- .sam_source(sam_fit)
  years <- x$years
  ages <- seq.int(x$conf$minAge, x$conf$maxAge)
  N <- stockassessment::ntable(sam_fit)
  F <- stockassessment::faytable(sam_fit)
  state <- function(m, se = NULL) {
    d <- .sam_long(m)
    names(d)[names(d) == "obs"] <- "value"
    out <- cbind(d[c("year", "age")], .sam_estimates(log(d$value), se, interval))
    out$is_proj <- FALSE
    out
  }
  N_tab <- state(N, if (!is.null(sam_fit$plsd$logN)) t(sam_fit$plsd$logN) else NULL)
  fkeys <- sam_fit$conf$keyLogFsta[which(sam_fit$data$fleetTypes == 0), ]
  fse <- if (!is.null(sam_fit$plsd$logF) && all(fkeys >= 0)) t(sam_fit$plsd$logF[fkeys + 1L, , drop = FALSE]) else NULL
  F_tab <- state(F, fse)
  M <- x$data$nm[as.character(years), as.character(ages), drop = FALSE]
  if (isTRUE(sam_fit$conf$mortalityModel > 0) && length(sam_fit$pl$logNM)) {
    M <- exp(sam_fit$pl$logNM[seq_along(years), seq_along(ages), drop = FALSE])
    dimnames(M) <- dimnames(N)
  }
  trend <- function(nm) {
    values <- sam_fit$sdrep$value
    i <- which(names(values) == nm)
    if (length(i) != length(years)) return(NULL)
    cbind(data.frame(year = years), .sam_estimates(values[i], sam_fit$sdrep$sd[i], interval), is_proj = FALSE)
  }
  plain <- function(values) data.frame(year = years, est = values, lwr = NA_real_, upr = NA_real_,
                                      se = NA_real_, se_scale = NA_character_, is_proj = FALSE)
  pop <- list(N = N_tab, F = F_tab, M = state(M), ssb = trend("logssb"),
              recruitment = trend("logR"), F_bar = trend("logfbar"), biomass = trend("logtsb"),
              abundance = plain(rowSums(N)))
  # These two summaries are unambiguous state transformations, without invented SEs.
  pop$S <- state(F / apply(F, 1, max))
  pop$M_bar <- plain(rowSums(M * N) / rowSums(N))
  obs <- sam_to_tam_obs(sam_fit)
  predictions <- x$raw
  predictions$pred <- if (length(sam_fit$rep$predObs) == nrow(predictions)) exp(sam_fit$rep$predObs) else NA_real_
  predictions$sd <- NA_real_
  for (f in x$fleets$fleet_id) {
    cov <- sam_fit$rep$obsCov[[f]]
    rows <- which(predictions$fleet_id == f)
    if (!is.matrix(cov)) next
    sd <- sqrt(diag(cov))[predictions$age[rows] - x$fleets$min_age[f] + 1L]
    w <- x$observation_weights[rows]
    if (length(w)) {
      weighted <- !is.na(w)
      flag <- sam_fit$conf$fixVarToWeight
      if (length(flag)) sd[weighted] <- if (rep(flag, length.out = nrow(x$fleets))[f] == 1)
        sqrt(w[weighted]) else sd[weighted] / sqrt(w[weighted])
    }
    # Prediction-dependent SD cannot be recovered from a static covariance table.
    link <- sam_fit$conf$predVarObsLink
    if (is.matrix(link)) sd[link[cbind(f, predictions$age[rows] - min(ages) + 1L)] >= 0 &
                               !is.na(link[cbind(f, predictions$age[rows] - min(ages) + 1L)])] <- NA_real_
    predictions$sd[rows] <- sd
  }
  key <- function(d) paste(d$year, d$fleet_id, d$age, sep = ":")
  obs_pred <- lapply(obs[c("catch", "index")], function(d) {
    i <- match(key(d), key(predictions))
    d$pred <- predictions$pred[i]
    d$sd <- predictions$sd[i]
    d$std_res <- ifelse(is.finite(d$obs) & d$obs > 0 & d$sd > 0,
                       (log(d$obs) - log(d$pred)) / d$sd, NA_real_)
    d$is_proj <- FALSE
    d
  })
  # tinyAM's dashboard yield definition uses stock weight, including for SAM.
  cw <- obs$weight$obs[match(paste(obs_pred$catch$year, obs_pred$catch$age),
                            paste(obs$weight$year, obs$weight$age))]
  yield <- function(v) plain(vapply(years, function(y) {
    i <- obs_pred$catch$year == y
    if (anyNA(v[i])) NA_real_ else sum(v[i] * cw[i])
  }, numeric(1)))
  pop$total_yield <- yield(obs_pred$catch$obs)
  pop$total_yield_pred <- yield(obs_pred$catch$pred)
  q <- obs_pred$index$q_key
  obs_pred$index$q <- ifelse(q == -1, 1, exp(sam_fit$pl$logFpar[pmax(q + 1L, 1L)]))
  # SAM conditional residuals retain SAM covariance/likelihood interpretation.
  for (f in x$fleets$fleet_id) {
    if (!identical(as.character(sam_fit$conf$obsLikelihoodFlag[f]), "LN") ||
        (!is.null(sam_fit$conf$fracMixObs) && sam_fit$conf$fracMixObs[f] != 0)) {
      for (nm in names(obs_pred)) obs_pred[[nm]]$std_res[obs_pred[[nm]]$fleet_id == f] <- NA_real_
    }
  }
  fixed <- data.frame(par = character(), coef = character(), est = numeric(), lwr = numeric(),
                      upr = numeric(), se = numeric(), se_scale = character())
  if (length(sam_fit$pl$logFpar)) {
    fixed <- cbind(data.frame(par = "q", coef = paste0("q_block", seq_along(sam_fit$pl$logFpar) - 1L)),
                   .sam_estimates(sam_fit$pl$logFpar, sam_fit$plsd$logFpar, interval))
  }
  random <- list(log_f = transform(F_tab, par = "f"), log_r = transform(N_tab[N_tab$age == min(ages), ], par = "r"))
  pop <- Filter(Negate(is.null), pop)
  attr(pop, "interval") <- attr(fixed, "interval") <- attr(random, "interval") <- interval
  list(pop = pop, obs_pred = obs_pred, fixed_par = fixed, random_par = random)
}

#' Prepare a SAM fit for comparison in tinyAM dashboards
#'
#' Display existing SAM outputs alongside tinyAM fits without refitting SAM or
#' manufacturing a tinyAM optimizer. This is a reporting object, not a `tam_fit`;
#' it cannot be simulated, updated or used for retrospective fitting.
#'
#' @param sam_fit A fitted SAM object.
#' @param interval Confidence level in `(0, 1)`. Log-scale uncertainty is
#'   transformed to natural-scale intervals, with `se` retained on the log scale.
#' @return A `tam_comparison` list containing `dat`, `pop`, `obs_pred`,
#'   `fixed_par`, `random_par`, source fit and reporting notes. Population tables
#'   retain SAM's native definitions; abundance/selectivity/Mbar are explicitly
#'   calculated from fitted states without newly fabricated uncertainty.
#' @details
#' Requires suggested package \pkg{stockassessment}. SAM's arithmetic Fbar and
#' spawning-time SSB must not be confused with tinyAM's native definitions.
#' Missing reports and unavailable uncertainty remain missing. Observation
#' residuals are conditional standardized log residuals, not one-step residuals.
#' @seealso [sam_to_tam_obs()], [sam_to_tam_settings()], [sam_to_tam_audit()], [vis_tam()]
#' @export
sam_to_tam_comparison <- function(sam_fit, interval = 0.95) {
  if (length(interval) != 1L || !is.finite(interval) || interval <= 0 || interval >= 1) cli::cli_abort("interval must be in (0, 1).")
  tabs <- .sam_comparison_tables(sam_fit, interval)
  x <- .sam_source(sam_fit)
  notes <- c("SAM reference: native arithmetic Fbar, spawning-time SSB and modeled biological surfaces.",
             "Original biological inputs are shown in the input tables; missing values are retained.",
             "State-based abundance, selectivity and Mbar have no newly calculated uncertainty.",
             "Only q is displayed in the fixed-parameter panel; other SAM parameters have different meanings.",
             "Dashboard yield uses original stock weights for both models; this is not SAM's native catch biomass.",
             "Unavailable reports are omitted from mixed panels; missing uncertainty is never filled.")
  structure(c(list(call = match.call(), dat = list(obs = sam_to_tam_obs(sam_fit), years = x$years,
    ages = seq.int(x$conf$minAge, x$conf$maxAge), is_proj = rep(FALSE, length(x$years))),
    source_fit = sam_fit, reporting_notes = notes), tabs), class = c("tam_comparison", "list"))
}

#' @export
update.tam_comparison <- function(object, ...) {
  cli::cli_abort("A tam_comparison is a reporting object and cannot be updated or fitted.")
}
