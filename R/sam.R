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

# Ordinary values from a SAM fit; no parsing or serialized optimizer calls.
.sam_source <- function(sam_fit) {
  if (!is.list(sam_fit) || !all(c("data", "conf", "pl", "opt") %in% names(sam_fit)) ||
      !is.list(sam_fit$opt) || length(sam_fit$opt$objective) != 1L ||
      !is.finite(sam_fit$opt$objective)) cli::cli_abort("Supply a fitted SAM object with data, conf, pl and a finite optimizer objective.")
  dat <- sam_fit$data
  conf <- sam_fit$conf
  bounds <- c(conf$minAge, conf$maxAge)
  if (length(bounds) != 2L || anyNA(bounds) || any(bounds != trunc(bounds)) ||
      bounds[1] < 0 || bounds[2] <= bounds[1]) cli::cli_abort("Invalid SAM modeled age range.")
  ages <- seq.int(bounds[1], bounds[2])
  years <- dat$years
  if (length(years) < 2L || anyNA(years) || any(years != trunc(years)) || any(diff(years) != 1)) cli::cli_abort("SAM years must be consecutive.")
  aux <- dat$aux
  if (!is.matrix(aux) || !all(c("year", "fleet", "age") %in% colnames(aux)) ||
      nrow(aux) != length(dat$logobs)) cli::cli_abort("SAM observation rows and logobs do not align.")
  raw <- data.frame(year = aux[, "year"], fleet_id = aux[, "fleet"], age = aux[, "age"], obs = exp(dat$logobs))
  if (anyDuplicated(raw[c("year", "fleet_id", "age")])) cli::cli_abort("Duplicate SAM year/fleet/age observation rows.")
  nf <- length(dat$fleetTypes)
  if (!nf || anyNA(raw[c("year", "fleet_id", "age")]) ||
      any(!raw$fleet_id %in% seq_len(nf)) || any(!raw$year %in% years)) cli::cli_abort("Invalid SAM observation identifiers.")
  fleet_names <- attr(dat, "fleetNames")
  if (length(fleet_names) != nf) fleet_names <- paste("Fleet", seq_len(nf))
  if (length(dat$sampleTimes) != nf) cli::cli_abort("SAM sampling times do not match fleets.")
  if (any(!raw$age %in% ages)) cli::cli_abort("SAM observation ages fall outside the modeled age range.")
  if (any(!seq_len(nf) %in% raw$fleet_id)) cli::cli_abort("Every SAM fleet must contain observation rows.")
  fleet_ages <- lapply(seq_len(nf), function(i) range(raw$age[raw$fleet_id == i]))
  fleets <- data.frame(fleet_id = seq_len(nf), fleet_name = fleet_names, fleet_type = dat$fleetTypes,
                       min_age = vapply(fleet_ages, min, numeric(1)), max_age = vapply(fleet_ages, max, numeric(1)),
                       samp_time = dat$sampleTimes)
  mats <- lapply(seq_len(nf), function(i) {
    d <- raw[raw$fleet_id == i, ]
    ys <- seq.int(min(d$year), max(d$year))
    aa <- seq.int(fleets$min_age[i], fleets$max_age[i])
    m <- matrix(NA_real_, length(ys), length(aa), dimnames = list(year = ys, age = aa))
    m[cbind(match(d$year, ys), match(d$age, aa))] <- d$obs
    # Effort was already applied by SAM and cannot be recovered from logobs.
    attr(m, "effort") <- rep(NA_real_, length(ys))
    m
  })
  names(mats) <- fleet_names
  data <- list(catch = mats[dat$fleetTypes == 0], surveys = mats[dat$fleetTypes == 2])
  fields <- c(sw = "stockMeanWeight", cw = "catchMeanWeight", mo = "propMat",
              nm = "natMor", pf = "propF", pm = "propM")
  for (nm in names(fields)) {
    m <- dat[[fields[[nm]]]]
    if (length(dim(m)) == 3L && dim(m)[3] == 1L) m <- m[, , 1]
    if (is.matrix(m)) {
      if (is.null(rownames(m))) rownames(m) <- utils::head(years, nrow(m))
      if (is.null(colnames(m))) colnames(m) <- utils::head(ages, ncol(m))
    }
    data[[nm]] <- m
  }
  list(data = data, fleets = fleets, conf = conf, years = years, raw = raw,
       observation_weights = dat$weight)
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
#' @seealso [sam_to_tam_obs()], [sam_to_tam_list()]
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

#' Convert a fitted SAM object to a list for comparison with tinyAM
#'
#' Convert a fitted SAM object into a list similar to a `tam_fit` object to ease
#' comparisons, especially through a dashboard made with [vis_tam()]. The list
#' contains SAM's fitted values and available uncertainty, arranged for use with
#' [tidy_tam()] and [vis_tam()]. The supplied SAM fit is left unchanged.
#'
#' @param sam_fit A fitted SAM object.
#' @param interval Confidence level in `(0, 1)`. Log-scale uncertainty is
#'   transformed to natural-scale intervals, with `se` retained on the log scale.
#' @return A list of class `tam_list` containing `dat`, `pop`, `obs_pred`,
#'   `fixed_par`, `random_par`, the source fit and notes explaining comparison
#'   differences. It resembles a `tam_fit` for reporting but cannot be fitted,
#'   updated, simulated or used for retrospective fitting.
#' @details
#' Requires suggested package \pkg{stockassessment}. SAM's own definitions are
#' retained: average fishing mortality gives equal weight to each age, and
#' spawning biomass accounts for mortality before spawning. tinyAM's reported
#' average fishing mortality is weighted by abundance, and its spawning biomass
#' is measured at the beginning of the year. Biological inputs may also differ
#' from SAM's fitted biological values.
#'
#' Missing values and unavailable uncertainty remain missing. Total abundance,
#' selectivity and average natural mortality are calculated from fitted states
#' without adding uncertainty estimates. Observation residuals describe the fit
#' to individual observations; they are not one-step-ahead residuals.
#' @seealso [sam_to_tam_obs()], [sam_to_tam_settings()], [sam_to_tam_audit()], [vis_tam()]
#' @export
sam_to_tam_list <- function(sam_fit, interval = 0.95) {
  if (length(interval) != 1L || !is.finite(interval) || interval <= 0 || interval >= 1) cli::cli_abort("interval must be in (0, 1).")
  tabs <- .sam_comparison_tables(sam_fit, interval)
  x <- .sam_source(sam_fit)
  notes <- c("Average fishing mortality (Fbar): SAM gives each age equal weight; tinyAM weights ages by their estimated abundance.",
             "Spawning biomass (SSB): SAM accounts for mortality before spawning; tinyAM reports biomass at the beginning of the year.",
             "Weight, maturity and natural mortality: SAM may estimate these values. The input tables show the original supplied values, including any missing values.",
             "Catch biomass (yield): this dashboard uses the original stock weights for both models. SAM's own catch biomass uses catch weights, which may differ.",
             "Uncertainty: total abundance, selectivity and average natural mortality are calculated from SAM's fitted states without new confidence intervals. Missing reports or uncertainty are left out or shown as missing.",
             "Parameters: only survey catchability (q) is shown among SAM's fixed parameters, because the other parameters do not correspond directly to tinyAM's.",
             "Residuals: the SAM residuals shown here describe the fit to individual observations; they are not one-step-ahead residuals.")
  structure(c(list(call = match.call(), dat = list(obs = sam_to_tam_obs(sam_fit), years = x$years,
                                                   ages = seq.int(x$conf$minAge, x$conf$maxAge), is_proj = rep(FALSE, length(x$years))),
                   source_fit = sam_fit, reporting_notes = notes), tabs), class = c("tam_list", "list"))
}

#' @export
update.tam_list <- function(object, ...) {
  cli::cli_abort("A tam_list is a reporting object and cannot be updated or fitted.")
}


