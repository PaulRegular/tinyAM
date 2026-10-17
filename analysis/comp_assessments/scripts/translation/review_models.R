source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/fit_default_first.R")
source("analysis/comp_assessments/tests/helper_stock.R")

review_database <- function() {
  .assessment_database(function(name) {
    file <- tempfile(fileext = ".csv")
    on.exit(unlink(file))
    path <- paste0("86e44bf:analysis/comp_assessments/database/", name)
    status <- system2("git", c("show", shQuote(path)), stdout = file)
    if (status != 0) cli::cli_abort("Could not read the pinned baseline database.")
    read.csv(file, na.strings = c("", "NA"), check.names = FALSE)
  }, source = "committed", commit = "86e44bf")
}

## Controlled model trials ----
review_model <- function(assessment_id, candidate, changes, rationale,
                         database = review_database(), observation_change = NULL) {
  record <- read_assessment(assessment_id, database)
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script))
  path <- file.path("analysis/comp_assessments/scripts/translation/stocks",
                    paste0(assessment_id, ".R"))
  status <- system2("git", c("show", shQuote(paste0("86e44bf:", path))), stdout = script)
  if (status != 0) cli::cli_abort("Could not read the pinned baseline stock script.")
  spec <- .test_stock(record, script = script)
  if (!is.null(observation_change)) spec$obs <- observation_change(spec$obs)
  settings <- spec$settings
  for (name in names(changes)) {
    base <- if (is.null(settings[[name]])) list() else settings[[name]]
    settings[[name]] <- utils::modifyList(base, changes[[name]], keep.null = TRUE)
  }
  directory <- file.path("analysis/comp_assessments/results/cache/model_review_20261016",
                         assessment_id, candidate)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(revision = database$commit, rationale = rationale,
               changes = changes, years = spec$years, ages = spec$ages,
               settings = settings), file.path(directory, "specification.rds"))
  model <- function(...) do.call(tinyAM::fit_tam,
    c(list(data = spec$obs, years = spec$years, ages = spec$ages,
           silent = TRUE), settings, list(...)))

  # Reuse only the recipe's existing conditional source-start construction.
  retries <- list()
  walk <- function(x) {
    if (missing(x) || !is.call(x)) return(invisible(NULL))
    if (identical(x[[1]], quote(fit_default_first))) retries[[length(retries)+1L]] <<- x$retry_start
    for (z in as.list(x)[-1]) walk(z)
  }
  for (x in parse(script)) walk(x)
  retry_start <- function() {
    starts <- NULL
    if (length(retries)) {
      e <- list2env(spec, parent = globalenv())
      starts <- eval(retries[[1]], e)()
    } else if (!is.null(spec$start_par)) {
      starts <- spec$start_par
    }
    if (is.null(starts)) cli::cli_abort("The recipe has no additional starting-value construction.")
    par <- tinyAM::make_par(do.call(tinyAM::prepare_tam,
      c(list(data = spec$obs, years = spec$years, ages = spec$ages), settings)))
    starts[intersect(names(starts), names(par)[lengths(par) > 0L])]
  }
  begun <- Sys.time()
  fit <- tryCatch(fit_default_first(model, retry_start,
    file.path(directory, "initialization_attempts.rds")), error = identity)
  if (inherits(fit, "error")) {
    saveRDS(list(error = conditionMessage(fit)), file.path(directory, "result.rds"))
    return(invisible(NULL))
  }
  saveRDS(fit, file.path(directory, "fit.rds"))
  diagnostics <- .assessment_diagnostics(assessment_id, database,
    if (isTRUE(fit$is_converged)) "converged" else "not_converged", fit,
    elapsed = as.numeric(difftime(Sys.time(), begun, units = "secs")))
  context <- .assessment_context(record, spec)
  comparison <- tryCatch(.assessment_comparison(fit, context), error = identity)
  result <- list(diagnostics = diagnostics, rationale = rationale, changes = changes)
  if (inherits(comparison, "error")) {
    result$comparison_error <- conditionMessage(comparison)
  } else {
    result$summary <- comparison$summary
    result$ref <- comparison$ref
  }
  saveRDS(result, file.path(directory, "result.rds"))
  print(diagnostics)
  if (!is.null(result$summary)) print(result$summary[c("metric",
    "median_absolute_percent_difference", "trend_correlation")])
  invisible(result)
}

## Reproduce a Pacific trial ----
review_candidate <- function(assessment_id, candidate) {
  if (candidate == "source_M_prior_BH" && assessment_id != "joint_hake_coastal_2025") {
    cli::cli_abort("The age-2 BH trial is specified only for coastal hake.")
  }
  change_data <- NULL
  if (candidate == "correlated_F") {
    changes <- list(F_settings=list(process="cor_rw",mu_form=NULL))
    rationale <- "Age-correlated annual F increments approximate shared fishing intensity/selectivity changes."
  } else if (candidate == "age_dependent_catch_sd") {
    changes <- list(catch_settings=list(sd_form=~age+I(age^2)))
    rationale <- "Estimate an age-dependent log-SD to approximate unequal relative composition precision."
  } else if (candidate == "annual_F_selectivity") {
    ascending <- assessment_id %in% c("afsc_yellowfin_sole_ebs_2024",
      "afsc_northern_rock_sole_ebs_2024", "afsc_alaska_plaice_bsai_2024", "dfo_pop_bc_2023")
    mean <- if (ascending) ~F_year+logistic(age) else ~F_year+factor(age)
    changes <- list(F_settings=list(process="iid",mu_form=mean))
    change_data <- function(obs) {
      obs$pop_covariates$F_year <- factor(obs$pop_covariates$year)
      obs
    }
    rationale <- "Shared annual fishing levels and an age-selectivity mean with IID residual F; terminal F-year covariates carry into projections."
  } else if (candidate == "source_M_q_priors" && assessment_id == "afsc_northern_rock_sole_ebs_2024") {
    changes <- list(M_settings=list(process="off",mu_form=~0+iid(M_prior_group,sd=.2),mu_supplied=~I(.15)),
      index_settings=list(q_form=~0+logistic(q_age)+iid(q_prior_group,sd=.2)))
    change_data <- function(obs) {
      obs$pop_covariates$M_prior_group <- factor("stock")
      obs$index$source_obs <- obs$index$obs
      obs$index$obs <- obs$index$obs/1.5
      obs$index$q_age <- pmin(obs$index$age,17)
      obs$index$q_prior_group <- factor("stock")
      obs
    }
    rationale <- "Native M and q log-penalties. Index/q are in units relative to the fixed prior center 1.5; log-likelihood residuals are unchanged. Original index is retained in source_obs, and the age-17 selectivity plateau is explicit."
  } else if (candidate %in% c("source_M_prior", "regularized_initial_abundance", "source_M_prior_BH")) {
    changes <- list()
    if (candidate == "regularized_initial_abundance") {
      changes$N_settings <- list(init="random")
      rationale <- "Regularize initial ages; adjacent-age survival through first-year Z differs from the source initialization."
    }
    if (candidate == "source_M_prior_BH") changes$rec_settings <- list(mu_form=~bh(ssb))
    if (candidate %in% c("source_M_prior", "source_M_prior_BH") || assessment_id == "afsc_northern_rock_sole_ebs_2024") {
      prior <- switch(assessment_id,
        joint_hake_coastal_2025=c(center=exp(-1.60944),sd=.1),
        afsc_northern_rock_sole_ebs_2024=c(center=.15,sd=.2),
        afsc_pollock_ai_2024=c(center=.2,sd=.2),
        afsc_pop_bsai_2024=c(center=.05*exp(-.05^2/2),sd=.05),
        dfo_pop_bc_2023=c(center=.06/sqrt(1+.3^2),sd=sqrt(log1p(.3^2))))
      if (is.null(prior)) cli::cli_abort("No verified M penalty is configured for this case.")
      changes$M_settings <- list(process="off",
        mu_form=stats::as.formula(paste0("~0+iid(M_prior_group,sd=",format(prior["sd"],digits=16),")")),
        mu_supplied=stats::as.formula(paste0("~I(",format(prior["center"],digits=16),")")))
      change_data <- function(obs) {
        obs$pop_covariates$M_prior_group <- factor("stock")
        obs
      }
      if (candidate == "source_M_prior_BH") rationale <- "Native M penalty plus a lag-two age-2 BH recruitment curve. Omitted juvenile fishing, h prior and bias/terminal rules remain approximations."
      if (candidate == "source_M_prior") rationale <- if (assessment_id == "dfo_pop_bc_2023") {
        "Positive lognormal moment approximation to the source normal M prior; bounds and exact density differ."
      } else "One shared log-M level with the native fixed-SD Gaussian penalty; Laplace integration and native bounds differ."
    }
  } else {
    cli::cli_abort("Unknown Pacific review candidate.")
  }
  review_model(assessment_id,candidate,changes,rationale,observation_change=change_data)
}
