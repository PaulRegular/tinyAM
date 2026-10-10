# Observation SD effects reuse the Gaussian formula-state compiler.
.parse_sd_formula <- function(formula, data, component) {
  design <- .parse_q_formula(formula, data)
  terms <- design$q_terms
  if (!is.null(design$q_mono_modmat) ||
      any(vapply(terms, function(x) x$type == "logistic", logical(1)))) {
    cli::cli_abort("{component}_settings$sd_form supports iid(), rw() and ar1(), not catchability curves.")
  }
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    term$id <- paste("sd", component, term$id, sep = "_")
    term$parameter <- paste0("eta_", term$id)
    if (!is.null(term$sd_parameter)) term$sd_parameter <- paste0("log_", term$id)
    if (!is.null(term$phi_parameter)) term$phi_parameter <- paste0("logit_phi_", term$id)
    terms[[i]] <- term
  }
  if (length(terms)) names(terms) <- vapply(terms, `[[`, character(1), "id")
  list(matrix = design$q_modmat, terms = terms)
}

.sd_effects <- function(par, dat, component, simulate = FALSE) {
  .q_effects(par, list(obs = list(index = dat$obs[[component]]),
    q_terms = dat[[paste0("sd_", component, "_terms")]]), simulate = simulate)
}

.check_sd_terms <- function(dat, component) {
  terms <- dat[[paste0("sd_", component, "_terms")]]
  if (!length(terms)) return(invisible(NULL))
  if (sum(vapply(terms, function(x) x$type %in% c("rw", "ar1"), logical(1))) > 1L) {
    cli::cli_abort("Use one RW/AR1 term in {component}_settings$sd_form; overlapping ordered SD processes have not been validated.")
  }
  # A random effect on SD is not additive observation variance on the mean scale.
  .check_q_terms(list(obs = list(index = dat$obs[[component]]), q_terms = terms,
    q_modmat = dat[[paste0("sd_", component, "_modmat")]],
    index_settings = list(q_link = NULL)))
}

.sd_process_advisories <- function(dat) {
  out <- data.frame(issue = character(), detail = character())
  for (component in c("catch", "index")) {
    terms <- dat[[paste0("sd_", component, "_terms")]]
    if (!length(terms)) next
    d <- dat$obs[[component]]
    observed <- !d$is_proj & is.finite(d$obs) & d$obs > 0
    for (term in terms) {
      for (group in term$groups) {
        rows <- group$rows[observed[group$rows] & term$multiplier[group$rows] != 0]
        support <- table(d[[term$variable]][rows])
        support <- support[support > 0]
        if (any(support < 3L)) out <- rbind(out, data.frame(issue = "SD_effect_replication",
          detail = paste0(term$id, ", group ", group$group,
            ": some SD effect levels have fewer than three observations. Random SD effects describe residual spread, not mean trends; sparse levels may be weakly estimated. Compare a simpler SD formula and inspect uncertainty.")))
        if ((!is.null(term$sd_parameter) || !is.null(term$phi_parameter)) && length(support) < 5L) {
          out <- rbind(out, data.frame(issue = "SD_effect_support", detail = paste0(term$id,
            ", group ", group$group, ": only ", length(support),
            " observed levels inform the SD process. Its process SD or AR1 correlation may be imprecise; consider a simpler formula or externally supported sd/phi.")))
        }
      }
    }
  }
  out
}
