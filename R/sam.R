.sam_long <- function(x) {
  data.frame(year = rep(as.integer(rownames(x)), times = ncol(x)),
             age = rep(as.integer(colnames(x)), each = nrow(x)), obs = as.vector(x))
}

.sam_key <- function(conf, field, fleet, ages) {
  key <- conf[[field]]
  if (is.null(key)) return(rep(NA_real_, length(ages)))
  expected <- conf$maxAge - conf$minAge + 1L
  if (!is.matrix(key) || ncol(key) != expected || nrow(key) < max(fleet)) {
    cli::cli_abort("SAM {field} dimensions do not match fleet and modeled ages.")
  }
  if (any(key != trunc(key), na.rm = TRUE) || any(key < -1, na.rm = TRUE)) {
    cli::cli_abort("SAM {field} keys must be integers, with -1 for inactive cells.")
  }
  key[cbind(fleet, ages - conf$minAge + 1L)]
}

.sam_block <- function(key) {
  key[!is.na(key) & key < 0] <- NA_real_
  factor(key, levels = sort(unique(key[!is.na(key)])))
}

.sam_formula <- function(block, setting) {
  if (anyNA(block) || !length(block)) return("No complete estimated block mapping")
  paste0(setting, " = ", if (length(unique(block)) == 1L) "~ 1" else
    paste0("~ 0 + ", if (grepl("q_form", setting, fixed = TRUE)) "q_block" else "sd_block"))
}

.sam_q_formula <- function(key) {
  if (!length(key) || anyNA(key)) return("No complete q-key mapping")
  if (all(key >= 0)) return(.sam_formula(.sam_block(key), "index_settings$q_form"))
  estimated <- sort(unique(key[key >= 0]))
  columns <- if (length(estimated)) paste0("q_key_", estimated) else character()
  paste0("index_settings$q_form = ~ 0", if (length(columns))
    paste0(" + ", paste(columns, collapse = " + ")) else "")
}

#' Translate standard SAM observations for tinyAM
#'
#' Use the same catch, survey, and biological inputs in tinyAM. Successful
#' conversion does not imply that the two assessment models are equivalent;
#' inspect [sam_to_tam_audit()] before constructing a model.
#'
#' @param sam_fit A fitted SAM object, for example from `stockassessment::fitfromweb()`.
#' @return A standard `obs` list with `catch`, `index`, `weight`, and `maturity`.
#'   Source fleet IDs, q/SD keys and factor blocks, plus-group flags, catch weight,
#'   natural mortality (`M_assumption`), and spawning timing are retained.
#' @details
#' Exactly one type-0 catch fleet and type-2 age-specific surveys are supported.
#' Multiple catches, biomass surveys, and other fleet types abort rather than
#' silently changing population or observation equations.
#'
#' A full catch grid includes `NA` for unobserved years, including terminal survey
#' years. Biological values must exist for each modeled year and age; they are
#' never extrapolated or filled. Missing biological inputs require an explicit
#' complete fitting period. SAM has already removed nonpositive observations;
#' their original values and effort denominators cannot be recovered. `q_block` and `sd_block` use the original global SAM
#' keys, preserving sharing across surveys. Separate tinyAM catch/index SD
#' parameter vectors cannot enforce sharing between those two components.
#' When q keys include `-1` (fixed q = 1), additional `q_key_0`, `q_key_1`, etc.
#' numeric indicator columns permit an exact no-intercept formula with zeros at
#' fixed-q cells. If every q key is `-1`, `q_form = ~ 0` fixes all q to 1.
#'
#' `propF` and `propM` are metadata only. Stock mean weight supplies `weight$obs`;
#' catch mean weight is retained as `catch_weight` and does not replace it.
#' @export
sam_to_tam_obs <- function(sam_fit) {
  x <- .sam_source(sam_fit)
  fleets <- x$fleets
  if (sum(fleets$fleet_type == 0) != 1L) {
    cli::cli_abort("Conversion requires one SAM catch fleet; multiple fleets cannot be silently combined.")
  }
  if (any(!fleets$fleet_type %in% c(0, 2))) {
    cli::cli_abort("Unsupported SAM fleet types: {paste(unique(fleets$fleet_type[!fleets$fleet_type %in% c(0, 2)]), collapse = ', ')}. Only catch type 0 and age-specific survey type 2 can be converted.")
  }
  if (!any(fleets$fleet_type == 2)) cli::cli_abort("At least one age-specific survey is required for tinyAM conversion.")
  conf <- x$conf
  bounds <- c(conf$minAge, conf$maxAge)
  if (!length(bounds)) bounds <- range(as.integer(colnames(x$data$catch[[1]])))
  ages <- seq.int(bounds[1], bounds[2])
  if (any(fleets$min_age < min(ages) | fleets$max_age > max(ages))) {
    cli::cli_abort("Source fleet ages fall outside the configuration; explicit age reduction is required before conversion.")
  }
  pg <- conf$maxAgePlusGroup
  if (length(pg) == nrow(fleets) && any(pg[fleets$fleet_type == 2] == 1 &
                                      fleets$max_age[fleets$fleet_type == 2] < max(ages))) {
    cli::cli_abort("Survey plus groups below the modeled terminal age cannot be translated faithfully.")
  }
  if (length(pg) && isTRUE(pg[which(fleets$fleet_type == 0)] == 0)) {
    cli::cli_abort("tinyAM requires a terminal catch plus group; this SAM structure cannot be translated faithfully.")
  }
  years <- x$years
  grid <- expand.grid(year = years, age = ages)
  grid <- grid[order(grid$year, grid$age), ]
  rownames(grid) <- NULL
  surface <- function(name, optional = FALSE) {
    m <- x$data[[name]]
    if (is.null(m) && optional) return(rep(NA_real_, nrow(grid)))
    if (is.matrix(m) && optional) {
      tab <- .sam_long(m)
      return(tab$obs[match(paste(grid$year, grid$age), paste(tab$year, tab$age))])
    }
    if (!is.matrix(m) || !all(as.character(years) %in% rownames(m)) ||
        !all(as.character(ages) %in% colnames(m))) {
      cli::cli_abort("SAM {name} must supply the complete modeled biological grid; no extrapolation is performed.")
    }
    as.vector(m[cbind(as.character(grid$year), as.character(grid$age))])
  }
  add_keys <- function(tab) {
    if (length(conf$minAge) && length(conf$maxAge)) {
      tab$sd_key <- .sam_key(conf, "keyVarObs", tab$fleet_id, tab$age)
      tab$sd_block <- .sam_block(tab$sd_key)
    }
    pg <- conf$maxAgePlusGroup
    tab$is_plus_group <- if (!length(pg)) NA else
      if (length(pg) == 1L) tab$age == max(ages) & pg == 1 else
        tab$age == fleets$max_age[tab$fleet_id] & pg[tab$fleet_id] == 1
    tab
  }
  catch <- grid
  source_catch <- .sam_long(x$data$catch[[1]])
  key <- function(tab) paste(tab$year, tab$age, sep = ":")
  catch$obs <- source_catch$obs[match(key(catch), key(source_catch))]
  catch$obs[catch$obs < 0 & !is.na(catch$obs)] <- NA_real_
  catch$fleet_id <- fleets$fleet_id[fleets$fleet_type == 0]
  catch$fleet_name <- fleets$fleet_name[catch$fleet_id]
  catch <- add_keys(catch)
  index <- do.call(rbind, lapply(seq_along(x$data$surveys), function(i) {
    m <- x$data$surveys[[i]]
    tab <- .sam_long(m)
    tab$fleet_id <- fleets$fleet_id[fleets$fleet_type == 2][i]
    tab$fleet_name <- fleets$fleet_name[tab$fleet_id]
    tab$survey <- tab$fleet_name
    tab$samp_time <- fleets$samp_time[tab$fleet_id]
    tab$effort <- attr(m, "effort")[match(tab$year, as.integer(rownames(m)))]
    tab <- add_keys(tab)
    if (length(conf$minAge) && length(conf$maxAge)) {
      tab$q_key <- .sam_key(conf, "keyLogFpar", tab$fleet_id, tab$age)
    }
    tab
  }))
  if ("q_key" %in% names(index)) {
    index$q_block <- .sam_block(index$q_key)
    if (any(index$q_key == -1, na.rm = TRUE)) for (k in sort(unique(index$q_key[index$q_key >= 0]))) {
      index[[paste0("q_key_", k)]] <- as.numeric(index$q_key == k)
    }
  }
  if ("sd_key" %in% names(index)) index$sd_block <- .sam_block(index$sd_key)
  index <- index[order(index$year, index$fleet_id, index$age), ]
  rownames(index) <- NULL
  weight <- transform(grid, obs = surface("sw"), catch_weight = surface("cw", TRUE),
                      M_assumption = surface("nm"), propF = surface("pf", TRUE),
                      propM = surface("pm", TRUE))
  maturity <- transform(grid, obs = surface("mo"))
  list(catch = catch, index = index, weight = weight, maturity = maturity)
}

#' Audit SAM assumptions against tinyAM
#'
#' Distinguish exact input/formula mappings from differences in latent states,
#' covariance, and biological summaries. An input conversion is not an assessment
#' replication. No fitting or automatic approximation is performed.
#'
#' @param sam_fit A fitted SAM object.
#' @param settings Named `fit_tam()` arguments, normally from [sam_to_tam_settings()].
#'   `NULL` audits the generated baseline.
#' @return A data frame with `component`, `sam_setting`, `sam_value`, `tam_status`,
#'   `tam_mapping`, and `notes`. Statuses are `supported`, `partially_supported`,
#'   `unsupported`, or `not_checked`. Missing and unreviewed fields remain
#'   `not_checked`; configuration defaults are never silently inferred.
#' @details
#' Exact q and observation-SD mappings use factor columns built directly from SAM
#' keys, with one coefficient per block (or an intercept for a single block).
#' Formulas do not couple latent F states or create process/observation covariance.
#' Shared catch/index SD keys require one parameter across two tinyAM parameter
#' vectors and therefore are only partially supported.
#' SAM integrates initial abundance states even without an initial-state density;
#' tinyAM's free initializer estimates initial abundance and recruitment as fixed
#' parameters. Their conditional state equations agree but marginal likelihoods
#' differ, so this boundary treatment is only partially supported.
#'
#' The audit targets standard Gaussian SAM models. It checks enabled mixtures,
#' density-dependent q, scaling, initial-state priors, biological process models,
#' and supplied observation attributes separately. Unknown fields are reported,
#' so this table is not a universal SAM compatibility certificate.
#' @seealso [sam_to_tam_obs()], [sam_to_tam_comparison()]
#' @export
sam_to_tam_audit <- function(sam_fit, settings = NULL) {
  x <- .sam_source(sam_fit)
  if (is.null(settings)) settings <- sam_to_tam_settings(sam_fit)
  conf <- x$conf
  fleets <- x$fleets
  nc <- sum(fleets$fleet_type == 0)
  ns <- sum(fleets$fleet_type == 2)
  rows <- list()
  value <- function(v) {
    if (is.null(v)) return("missing")
    txt <- as.character(v)
    if (length(txt) > 24L) txt <- c(txt[seq_len(24)], paste0("... (", length(v), " values)"))
    paste(txt, collapse = ", ")
  }
  add <- function(component, setting, status, mapping = "", notes = "", v = conf[[setting]]) {
    if (is.null(v)) status <- "not_checked"
    rows[[length(rows) + 1L]] <<- data.frame(component = component, sam_setting = setting,
      sam_value = value(v), tam_status = status, tam_mapping = mapping, notes = notes)
  }
  same <- function(v) length(unique(v[!is.na(v) & v >= 0])) <= 1L
  add("age range", "minAge/maxAge", "supported", "ages = minAge:maxAge",
      "Conversion rejects source ages outside this range rather than aggregating silently.",
      if (length(conf$minAge) && length(conf$maxAge)) c(conf$minAge, conf$maxAge) else NULL)
  pg <- conf$maxAgePlusGroup
  pg_ok <- length(pg) > 0 && pg[1] == 1 &&
    (length(pg) == 1L || all(pg[fleets$fleet_type == 2] == as.integer(fleets$max_age[fleets$fleet_type == 2] == conf$maxAge)))
  add("plus groups", "maxAgePlusGroup", if (isTRUE(pg_ok)) "supported" else "unsupported",
      "tinyAM terminal plus group", "tinyAM always has a terminal plus group; lower-age survey plus groups cannot be translated as ordinary age observations.")
  add("catch fleets", "fleetTypes (catch)", if (nc == 1L) "supported" else "unsupported",
      "one type-0 catch table", "Multiple catch fleets cannot share tinyAM's single F surface without changing assumptions.", fleets$fleet_type[fleets$fleet_type %in% c(0, 1, 7)])
  add("survey fleets", "fleetTypes (survey)", if (ns > 0 && all(fleets$fleet_type %in% c(0, 2))) "supported" else "unsupported",
      "index$survey and index$samp_time", "Type 2 only; stored timing and already normalized observations are retained.", fleets$fleet_type[fleets$fleet_type != 0])
  fkeys <- conf$keyLogFsta
  frow <- if (is.matrix(fkeys) && nc == 1L) fkeys[fleets$fleet_type == 0, ] else numeric()
  f_ok <- length(frow) > 0 && !anyNA(frow) && all(frow >= 0) && !anyDuplicated(frow)
  add("F states", "keyLogFsta", if (f_ok) "supported" else "unsupported",
      "one independent latent log_f column per modeled age",
      "Repeated keys share the SAME state in SAM; mu_form shares only expected/mean structure and cannot enforce this equality.")
  add("F innovation correlation", "corFlag", if (isTRUE(all(conf$corFlag == 0))) "supported" else "unsupported",
      'F_settings = list(process = "rw", mu_form = NULL)',
      "SAM 0=independent, 1=compound symmetry, 2=AR1 across innovations, 3=separable, 4=0.99 correlation. tinyAM rw has independent increments; stationary ar1 is not equivalent.")
  qkeys <- conf$keyLogFpar
  q <- if (is.matrix(qkeys)) unlist(lapply(which(fleets$fleet_type == 2), function(i)
    qkeys[i, seq.int(fleets$min_age[i], fleets$max_age[i]) - conf$minAge + 1L])) else NULL
  q_ok <- !is.null(q) && length(q) && !anyNA(q) && all(q >= -1) && all(q == trunc(q))
  add("catchability", "keyLogFpar", if (q_ok) "supported" else "unsupported",
      if (q_ok) .sam_q_formula(q) else "",
      "q_block uses original global keys, including sharing across surveys. Fixed q=1 cells (-1) use zero rows in numeric indicator columns, without an intercept.")
  fv <- if (is.matrix(conf$keyVarF) && length(frow)) conf$keyVarF[fleets$fleet_type == 0, match(unique(frow[frow >= 0]), frow)] else numeric()
  fv_ok <- length(fv) && !anyNA(fv) && all(fv >= 0) && same(fv)
  add("F process variance", "keyVarF", if (nc == 1L && fv_ok) "supported" else "unsupported",
      "one tinyAM log_sd_f", "SAM selects variance for each unique F state from its first fleet-age occurrence; different active keys need different process SDs.")
  nkeys <- conf$keyVarLogN
  n_ok <- length(nkeys) > 1L && !anyNA(nkeys) && all(nkeys >= 0) && same(nkeys[-1]) && !nkeys[1] %in% nkeys[-1]
  add("N process variance", "keyVarLogN", if (n_ok) "supported" else "unsupported",
      'N_settings = list(process = "iid", init = "free"); separate sd_r and sd_n',
      "Exact variance grouping needs one survival SD and a separate recruitment SD; tinyAM cannot tie those two parameters or fit age-specific survival SDs.")
  obskeys <- conf$keyVarObs
  ckeys <- if (is.matrix(obskeys)) obskeys[fleets$fleet_type == 0, , drop = FALSE] else NULL
  ikeys <- if (is.matrix(obskeys)) obskeys[fleets$fleet_type == 2, , drop = FALSE] else NULL
  shared <- intersect(ckeys[!is.na(ckeys) & ckeys >= 0], ikeys[!is.na(ikeys) & ikeys >= 0])
  ov <- if (is.matrix(obskeys)) unlist(lapply(which(fleets$fleet_type %in% c(0, 2)), function(i)
    obskeys[i, seq.int(fleets$min_age[i], fleets$max_age[i]) - conf$minAge + 1L])) else numeric()
  ov_ok <- length(ov) && !anyNA(ov) && all(ov >= 0)
  add("observation variance", "keyVarObs", if (!ov_ok) "unsupported" else if (length(shared)) "partially_supported" else "supported",
      paste(.sam_formula(.sam_block(as.vector(ckeys)[as.vector(ckeys) >= 0]), "catch_settings$sd_form"),
            .sam_formula(.sam_block(as.vector(ikeys)[as.vector(ikeys) >= 0]), "index_settings$sd_form"), sep = "; "),
      "sd_block retains source keys. Sharing is exact within each table, but shared catch/index keys cannot be enforced by separate SD parameter vectors.")
  obs_fleets <- which(fleets$fleet_type %in% c(0, 2, 3, 6))
  independent <- isTRUE(all(as.character(conf$obsCorStruct)[obs_fleets] == "ID"))
  add("observation covariance", "obsCorStruct", if (independent) "supported" else "unsupported",
      "independent observation errors", "AR and US covariance cannot be represented by sd_form.")
  add("observation correlations", "keyCorObs", if (independent) "supported" else "unsupported",
      "inactive when obsCorStruct=ID", "Correlation keys matter only for enabled AR observation covariance.")
  add("observation likelihood", "obsLikelihoodFlag", if (isTRUE(all(as.character(conf$obsLikelihoodFlag)[obs_fleets] == "LN"))) "supported" else "unsupported",
      "independent lognormal likelihood", "ALN models proportions plus totals; Gaussian-mixture and observation-weight options are audited separately.")
  add("recruitment", "stockRecruitmentModelCode", if (identical(as.numeric(conf$stockRecruitmentModelCode), 0)) "supported" else "unsupported",
      "random walk on log recruitment", "tinyAM's recruitment is always a log random walk; N_settings$process controls survival residuals.")
  add("natural mortality", "nm.dat", if (all(is.finite(x$data$nm) & x$data$nm > 0)) "supported" else "unsupported", 'M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption)',
      "Supplied year-age natural mortality, with no estimated mean intercept. tinyAM's supplied log-M surface requires positive M; zero-M cases need separate treatment.", x$data$nm)
  add("maturity", "mo.dat", "supported", "maturity$obs", "Unchanged proportions; modeled maturity is a separate option.", x$data$mo)
  add("stock weight", "sw.dat", "supported", "weight$obs", "Unchanged stock mean weight, with original source units.", x$data$sw)
  add("catch weight", "cw.dat", "partially_supported", "weight$catch_weight metadata",
      "Catch numbers use the same Baranov equation, but tinyAM yield uses stock weight; compare catch biomass using retained catch weight separately.", x$data$cw)
  zero_timing <- !is.null(x$data$pf) && !is.null(x$data$pm) && isTRUE(all(x$data$pf == 0)) && isTRUE(all(x$data$pm == 0))
  add("spawning timing", "propF/propM", if (zero_timing) "supported" else "unsupported",
      "weight$propF and weight$propM metadata", "SAM SSB includes exp(-F*propF-M*propM); tinyAM SSB is at the beginning of the year.",
      if (is.null(x$data$pf) || is.null(x$data$pm)) NULL else c(x$data$pf, x$data$pm))
  add("Fbar", "fbarRange", "partially_supported", "F_settings$mean_ages = seq(min(fbarRange), max(fbarRange)); compare rowMeans(F[, ages])",
      "The age range maps exactly; SAM's arithmetic Fbar differs from tinyAM's abundance-weighted reported F_bar.")
  add("catch scaling", "noScaledYears", if (length(conf$noScaledYears) && conf$noScaledYears == 0) "supported" else "unsupported",
      "no estimated catch scaling", "SAM subtracts logScale[keyParScaledYA] from catch predictions in keyScaledYears; tinyAM has no matching parameter.")
  add("catch scaling years", "keyScaledYears", if (!length(conf$keyScaledYears)) "supported" else "unsupported",
      "inactive only when catch scaling is off", "Scaling years must not be replaced by fixed data adjustments to unknown fitted scales.")
  add("catch scaling keys", "keyParScaledYA", if (!length(conf$keyParScaledYA)) "supported" else "unsupported",
      "inactive only when catch scaling is off", "Keys share estimated scaling parameters across years/ages.")
  obs_attributes <- if (any(!is.na(x$observation_weights))) "weight" else character()
  add("observation weights", "fixVarToWeight", if (!length(obs_attributes)) "supported" else "partially_supported",
      "supplied SD = sqrt(weight) when flag=1; otherwise SD offset = 1/sqrt(weight)",
      "Inactive when no observation-weight attributes are supplied. Otherwise exact only with known positive weights, compatible SD groups and independent LN errors. Supplied observation weights are retained for audit; the baseline omits their likelihood effects.")
  add("observation attributes", "supplied observation attributes", if (!length(obs_attributes)) "supported" else "not_checked",
      "no supplied covariance/weight/partition attributes", "Any supplied weights require separate validation; they must not be confused with keyVarObs.",
      if (length(obs_attributes)) obs_attributes else "none")
  add("initial abundance", "initState", if (identical(as.numeric(conf$initState), 0)) "partially_supported" else "unsupported",
      'N_settings$init = "free"; no starting-state density',
      "Both omit the initial-state density, but SAM integrates all initial logN states while tinyAM estimates log_r0 and free log_n0 as fixed parameters. Conditional equations agree; marginal likelihoods differ.")
  for (nm in c("fracMixF", "fracMixN", "fracMixObs", "stockWeightModel", "catchWeightModel", "matureModel", "mortalityModel", "logNMeanAssumption")) {
    add("additional model option", nm, if (isTRUE(all(conf[[nm]] == 0))) "supported" else "unsupported",
        "zero/off option",
        "Only the zero/off option is matched; nonzero values require a separate mathematical mapping.")
  }
  add("density-dependent q", "keyQpow", if (isTRUE(all(conf$keyQpow < 0))) "supported" else "unsupported",
      "fixed linear q*N relationship", "Nonnegative keys estimate powers of predicted abundance; ordinary q_form does not do this.")
  add("prediction-dependent observation SD", "predVarObsLink", if (all(is.na(conf$predVarObsLink) | conf$predVarObsLink < 0)) "supported" else "unsupported",
      "inactive links", "sd_form can use known covariates, not an estimated power of the model prediction.")
  add("extra observation SD", "keyXtraSd", if (!length(conf$keyXtraSd)) "supported" else "not_checked",
      "no extra year-fleet-age SD keys", "Nonempty keys require inspection of the specific sharing and supplied-weight interactions.")
  for (nm in c("keyVarLogP", "keyBiomassTreat", "constRecBreaks")) {
    inactive <- switch(nm, keyVarLogP = !length(conf[[nm]]) && !any(fleets$fleet_type == 6),
      keyBiomassTreat = isTRUE(all(conf[[nm]] == -1)) && !any(fleets$fleet_type == 3),
      constRecBreaks = !length(conf[[nm]]) && identical(as.numeric(conf$stockRecruitmentModelCode), 0))
    add("additional model option", nm, if (inactive) "supported" else "not_checked",
        "inactive for the selected fleets/recruitment model", "Active configurations require a separate mapping.")
  }
  bio_keys <- c(keyStockWeightMean = "stockWeightModel", keyStockWeightObsVar = "stockWeightModel",
    keyCatchWeightMean = "catchWeightModel", keyCatchWeightObsVar = "catchWeightModel",
    keyMatureMean = "matureModel", keyMortalityMean = "mortalityModel", keyMortalityObsVar = "mortalityModel")
  for (nm in names(bio_keys)) {
    inactive <- identical(as.numeric(conf[[bio_keys[[nm]]]]), 0) && all(is.na(conf[[nm]]))
    add("biological process keys", nm, if (inactive) "supported" else "not_checked",
        paste0("inactive when ", bio_keys[[nm]], " = 0"), "Keys do not estimate biology when the corresponding process model is off.")
  }
  checked <- unlist(strsplit(vapply(rows, function(r) r$sam_setting, character(1)), "/", fixed = TRUE))
  for (nm in setdiff(names(conf), checked)) {
    add("unreviewed configuration", nm, "not_checked", notes = "Retained without inferring compatibility; inspect if active.")
  }
  out <- do.call(rbind, rows)
  .sam_audit_settings(out, sam_fit, settings)
}
