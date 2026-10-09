# Mean effects use the same unique-state compiler and densities as q effects.
.mean_formula_data <- function(dat, component) {
  d <- dat$obs[[if (component == "F") "catch" else "weight"]]
  d$obs <- 1
  d
}

.parse_mean_formula <- function(formula, dat, component) {
  d <- .mean_formula_data(dat, component)
  design <- .parse_q_formula(formula, d)
  terms <- design$q_terms
  if (!is.null(design$q_mono_modmat) ||
      any(vapply(terms, function(x) x$type == "logistic", logical(1)))) {
    cli::cli_abort("{component}_settings$mu_form supports iid(), rw() and ar1(), not catchability curves.")
  }
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term$id <- paste(component, term$id, sep = "_")
    term$parameter <- paste0("eta_mu_", term$id)
    if (!is.null(term$sd_parameter)) term$sd_parameter <- paste0("log_sd_mu_", term$id)
    if (!is.null(term$phi_parameter)) term$phi_parameter <- paste0("logit_phi_mu_", term$id)
    terms[[i]] <- term
  }
  if (length(terms)) names(terms) <- vapply(terms, `[[`, character(1), "id")
  list(matrix = design$q_modmat, terms = terms)
}

.check_mean_terms <- function(dat, component) {
  terms <- dat[[paste0(component, "_terms")]]
  if (!length(terms)) return(invisible(NULL))
  settings <- dat[[paste0(component, "_settings")]]
  if (sum(vapply(terms, function(x) x$type %in% c("rw", "ar1"), logical(1))) > 1L) {
    cli::cli_abort(c("Use one temporal process in the {component} mean.",
      "i" = "Multiple RW/AR1 mean terms may explain the same changes in mortality; their separate contributions have not been validated.",
      "i" = "Retain one rw() or ar1() term and use IID residuals for independent fluctuations."))
  }
  if (!settings$process %in% c("iid", "off")) {
    alternative <- if (component == "M") ", or use process = 'off' to vary only the M mean" else ""
    cli::cli_abort(c("Structured {component} means currently require IID residuals{if (component == 'M') ' or process = off' else ''}.",
      "i" = "Structured means and temporal residual processes may explain overlapping mortality variation; this combination has not been validated.",
      "i" = "Set process = 'iid'{alternative}."))
  }
  d <- .mean_formula_data(dat, component)
  fixed <- dat[[paste0(component, "_modmat")]]
  proxy <- list(obs = list(index = d), q_terms = terms, q_modmat = fixed,
    mean_component = component,
    index_settings = list(q_link = "log"),
    sd_index_modmat = matrix(if (settings$process == "iid") 1 else numeric(), nrow(d),
                            if (settings$process == "iid") 1L else 0L))
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
          "i" = "Share mean effects across independent age blocks, supply sd in iid(), or use M_settings$process = 'off'."))
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
  if (identical(dat$M_settings$process, "iid") && estimated(dat$M_terms)) {
    out <- rbind(out, data.frame(issue = "M_variance_separation", detail = paste(
      "Both M mean-process and IID residual SDs are estimated. Simulations showed that small shared M variation can be difficult to separate from residual variation.",
      "Start with M_settings$process = 'off', or supply sd in the mean term when external information supports it; check recovery before interpreting both components.")))
  }
  if (estimated(dat$F_terms) && estimated(dat$M_terms)) {
    out <- rbind(out, data.frame(issue = "joint_mortality_means", detail = paste(
      "F and M mean-process SDs are both estimated. The recovery study varied one mortality mean at a time; joint estimation has not been validated.",
      "Compare with models varying only one mean and examine sensitivity of F, M and abundance.")))
  }
  out
}

.mean_uncertainty_advisories <- function(fit) {
  out <- data.frame(issue = character(), detail = character())
  terms <- c(fit$dat$F_terms, fit$dat$M_terms)
  sdr <- fit[["sdrep"]]
  if (!length(terms) || !is.list(sdr) || !isTRUE(sdr$pdHess)) return(out)
  p <- sdr$par.fixed
  covariance <- sdr$cov.fixed
  if (!is.matrix(covariance) || nrow(covariance) != length(p)) return(out)
  variance <- diag(covariance)
  variance[!is.finite(variance) | variance < 0] <- NA_real_
  se <- sqrt(variance)
  half_width <- stats::qnorm(.975) * se
  wide_sd <- character()
  for (component in c("F", "M")) {
    component_terms <- fit$dat[[paste0(component, "_terms")]]
    if (!length(component_terms)) next
    sd_names <- c(unlist(lapply(component_terms, `[[`, "sd_parameter")),
      if (identical(fit$dat[[paste0(component, "_settings")]]$process, "iid"))
        paste0("log_sd_", tolower(component)))
    for (term in component_terms) {
      i <- match(term$phi_parameter, names(p))
      if (!length(i) || is.na(i) || !is.finite(half_width[i]) || !is.finite(p[i])) next
      width <- diff(stats::plogis(p[i] + c(-1, 1) * half_width[i]))
      if (width > .5) out <- rbind(out, data.frame(issue = "mortality_AR1_uncertainty",
        detail = paste0(term$id, ": the 95% AR1 correlation interval spans more than 0.5. ",
          "Persistence is weakly estimated; compare with IID/RW means or a scientifically supported fixed phi before interpreting it.")))
    }
    i <- match(sd_names, names(p))
    i <- i[!is.na(i)]
    wide <- is.finite(half_width[i]) & is.finite(p[i]) & 2 * half_width[i] > log(10)
    wide_sd <- c(wide_sd, names(p)[i[wide]])
  }
  if (length(wide_sd)) out <- rbind(out, data.frame(issue = "mortality_SD_uncertainty",
    detail = paste0("95% intervals span more than tenfold for ", paste(unique(wide_sd), collapse = ", "),
      ". These variance components are weakly estimated; compare simpler mean/residual structures or use externally supported mean-term SDs.")))
  out
}

.mean_effects <- function(par, dat, component, simulate = FALSE) {
  proxy <- list(obs = list(index = .mean_formula_data(dat, component)),
                q_terms = dat[[paste0(component, "_terms")]])
  .q_effects(par, proxy, simulate = simulate)
}

.formula_terms <- function(dat) c(dat$q_terms, dat$F_terms, dat$M_terms)

.formula_random_parameters <- function(dat) {
  unlist(lapply(.formula_terms(dat), `[[`, "parameter"), use.names = FALSE)
}

.tidy_formula_effects <- function(fit, interval = .95) {
  all <- list(q = .tidy_q_effects(fit, interval))
  for (component in c("F", "M")) {
    if (!length(fit$dat[[paste0(component, "_terms")]])) next
    proxy <- fit
    proxy$dat$q_terms <- fit$dat[[paste0(component, "_terms")]]
    proxy$dat$obs$index <- .mean_formula_data(fit$dat, component)
    all[[component]] <- .tidy_q_effects(proxy, interval,
      increment_report = paste0("eta_mu_", component, "_increments"))
  }
  out <- list(levels = list(), increments = list(), contributions = list())
  for (component in names(all)) for (part in names(out)) {
    tables <- lapply(all[[component]][[part]], function(tab) {
      tab$component <- component
      tab
    })
    out[[part]] <- c(out[[part]], tables)
  }
  if (!length(out$levels)) return(list())
  attr(out, "interval") <- interval
  out
}
