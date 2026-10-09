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
    cli::cli_abort("Use one temporal process in the {component} mean. Multiple temporal mean processes need a separate identifiability review.")
  }
  if (!settings$process %in% c("iid", "off")) {
    cli::cli_abort("Structured {component} means currently require IID residuals{if (component == 'M') ' or process = off' else ''}. Combining temporal mean and temporal residual processes is not enabled.")
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
        cli::cli_abort("Unreplicated IID M mean effects with IID residuals are not enabled with both SDs estimated. Share effects across independent age blocks or supply one SD.")
      }
    }
  }
  .check_q_terms(proxy)
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
