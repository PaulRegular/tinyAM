#' Audit source-assessment assumptions against current tinyAM structure
#'
#' This analysis-local audit classifies documented source features without
#' selecting tinyAM settings. A partial match means the broad feature exists
#' but its likelihood, process, or boundary treatment differs.
#'
#' @param assessment_id One assessment identifier from `assessments.csv`.
#' @param assumptions Canonical `assumptions.csv` data.
#'
#' @return A long data frame retaining the source setting and adding a tinyAM
#'   support class and plain-language explanation.
audit_assumptions <- function(assessment_id, assumptions) {
  required <- c("assessment_id", "component", "setting", "value",
                "source_reference", "notes")
  missing <- setdiff(required, names(assumptions))
  if (length(missing)) {
    .translation_abort(paste("assumptions is missing required columns:",
                             paste(missing, collapse = ", ")))
  }
  if (length(assessment_id) != 1L || is.na(assessment_id) ||
      !nzchar(assessment_id)) {
    .translation_abort("assessment_id must be one non-empty value.")
  }
  x <- assumptions[!is.na(assumptions$assessment_id) &
                     assumptions$assessment_id == assessment_id, , drop = FALSE]
  if (!nrow(x)) .translation_abort(paste("No assumptions recorded for", assessment_id))

  for (column in c("fleet", "survey", "sex", "region", "season")) {
    if (!column %in% names(x)) x[[column]] <- NA_character_
  }
  x$tinyam_support <- "not_checked"
  x$audit_notes <- "Review this source description against the chosen tinyAM model before fitting."
  key <- paste(x$component, x$setting, sep = "|")
  value <- as.character(x$value)
  known <- !is.na(value) & nzchar(value) & !tolower(value) %in% c("unknown", "not reported")

  mark <- function(rows, support, note) {
    rows[is.na(rows)] <- FALSE
    rows <- rows & known
    x$tinyam_support[rows] <<- support
    x$audit_notes[rows] <<- note
  }

  mark(key == "population|age_range", "supported",
       "tinyAM accepts explicit integer model ages; population dynamics still need comparison.")
  mark(key == "population|plus_group", "partially_supported",
       "tinyAM supports a plus age, but its construction and the source age range must match.")
  mark(key == "N|cohort_survival", "partially_supported",
       "tinyAM uses cohort survival with Z = F + M; recruitment and initial-age treatments may differ.")
  mark(key == "N|initial_abundance", "partially_supported",
       "tinyAM offers exp, free, and random initial abundance; compare the source initialization directly.")
  mark(key == "recruitment|process", "partially_supported",
       "tinyAM has recruitment and cohort processes, but these do not automatically reproduce the source process.")
  mark(key == "recruitment|deviation_distribution", "partially_supported",
       "tinyAM supports IID, RW and stationary AR1 recruitment residuals. An AR1 recruitment rate multiplied by SSB is a different equation; verify the source recruitment age, lag and bias correction.")
  recruitment_curve <- x$component == "recruitment" &
    x$setting %in% c("stock_recruit", "stock_recruitment")
  mark(recruitment_curve & grepl("Beverton.Holt|Ricker", value, ignore.case = TRUE),
       "partially_supported",
       "BH and Ricker median-recruitment curves are available with IID or AR1 residuals. Check age, parent-SSB timing, active parameters and source bias correction; a fixed steepness of one is not a free BH curve.")
  mark(recruitment_curve & grepl("hockey|segmented", value, ignore.case = TRUE),
       "unsupported", "tinyAM does not implement the source hockey-stick recruitment curve; BH and Ricker are not equivalent replacements.")
  mark(key == "F|fishing_mortality", "partially_supported",
       "tinyAM estimates an age-year F surface; fleet-specific fully recruited F and selectivity are not identical.")
  mark(key == "F|fishery_selectivity", "partially_supported",
       "tinyAM can estimate flexible F-at-age but does not reproduce this period-specific selectivity structure directly.")
  mark(key == "M|process", "partially_supported",
       "tinyAM supports age-blocked IID, RW and AR1 M states and structured mean formulas. Check the source density, first-state treatment, fixed process SDs and priors separately.")
  mark(key == "M|process_sd", "partially_supported",
       "M residual SDs can be estimated with ordinary age-sharing formulas, constant within fitted M state blocks. A mean-formula RW with supplied SD and process off can represent fixed increments; verify initial-state treatment and source priors separately.")
  mark(key == "M|initial_priors", "unsupported",
       "tinyAM has no matching prior distribution for the starting M levels; reported source means may only be used as starting values.")
  mark(key == "M|natural_mortality", "partially_supported",
       "Fixed M can be supplied directly; estimated source M requires an explicit simplified tinyAM process choice.")
  mark(key == "M|age_groups", "partially_supported",
       "tinyAM can share M states over age blocks; verify the selected ages and plus-group boundary.")
  mark(key == "index|sampling_time", "supported",
       "tinyAM uses a fraction of year for index timing; an approximate season midpoint must remain labelled as approximate.")
  approximate_time <- key == "index|sampling_time" & grepl("approx", value, ignore.case = TRUE)
  approximate_time[is.na(approximate_time)] <- FALSE
  x$tinyam_support[approximate_time] <- "partially_supported"
  x$audit_notes[approximate_time] <- "The source month is known, but tinyAM uses an approximate within-month timing value."
  unknown_time <- key == "index|sampling_time" &
    grepl("unknown|not reported", value, ignore.case = TRUE)
  unknown_time[is.na(unknown_time)] <- FALSE
  x$tinyam_support[unknown_time] <- "unsupported"
  x$audit_notes[unknown_time] <- "No sampling time is available; a documented timing choice is required before this index can be fitted."
  mark(key == "index|age_aggregation", "partially_supported",
       "tinyAM expects age-specific indices; reconstructing them from a total and age composition changes the observation model.")
  mark(key == "index|likelihood", "partially_supported",
       "tinyAM fits age-specific lognormal observations; aggregate-index and composition likelihoods are different.")
  mark(key == "index|catchability_age_structure", "partially_supported",
       "tinyAM supports fixed sharing groups, monotone/logistic age curves and IID/RW/AR1 catchability effects. Verify the documented keys, normalization and q link before selecting a formula.")
  mark(key == "index|age_composition_likelihood", "partially_supported",
       "tinyAM does not fit a multivariate composition likelihood; any abundance-at-age conversion is an approximation.")
  mark(key == "catch|catch_likelihood", "partially_supported",
       "tinyAM uses lognormal catch-at-age observations, unlike a source model that fits totals and age composition separately.")
  observation_correlation <- x$component %in% c("catch", "index", "observation") &
    grepl("correlation", x$setting)
  mark(observation_correlation & grepl("independent|^ID$", value, ignore.case = TRUE),
       "supported", "Independent lognormal observation errors are available; their SD can have fixed or structured formula effects.")
  mark(observation_correlation & grepl("AR|covariance|correlated", value, ignore.case = TRUE) &
         !grepl("independent", value, ignore.case = TRUE),
       "unsupported", "Structured mean or SD effects do not reproduce correlation between observation residuals; tinyAM has no matching multivariate observation likelihood.")
  mark(x$component %in% c("catch", "index", "observation") &
         x$setting %in% c("observation_process", "observation_likelihood") &
         grepl("AR|covariance|correlated", value, ignore.case = TRUE) &
         !grepl("independent", value, ignore.case = TRUE),
       "partially_supported", "Lognormal observations are available, but the source residual covariance is not reproduced by mean or SD formula effects.")
  mark(x$component == "F" & grepl("correlation|process", x$setting) &
         grepl("correlated.*increment|increment.*correlation|age.correlated.*random.walk", value, ignore.case = TRUE),
       "partially_supported", "F process cor_rw supports AR1 age-correlated annual RW innovations with marginal age-specific SDs. Verify the source covariance, SD sharing and first-state treatment; arbitrary covariance and correlated N/M RWs are not supported.")
  mark(key == "biology|weight_at_age", "supported",
       "tinyAM accepts annual weight-at-age when the required model grid is complete.")
  mark(key == "weight|weight_at_age", "supported",
       "tinyAM accepts annual weight-at-age when the required model grid is complete.")
  maturity <- key == "biological|maturity_at_age"
  cohort_indexed_maturity <- grepl(
    "birth[ _-]cohort|cohort[ _-]year|cohort-indexed",
    paste(value, x$notes), ignore.case = TRUE
  )
  mark(maturity, "supported",
       "The reported maturity values are indexed by calendar year and age; a cohort effect in the estimation model does not change the table indexing.")
  mark(maturity & cohort_indexed_maturity, "partially_supported",
       "Maturity values indexed by birth cohort need a defensible mapping to calendar year before tinyAM use.")

  x[c("component", "fleet", "survey", "setting", "value", "tinyam_support",
      "source_reference", "notes", "audit_notes")]
}
