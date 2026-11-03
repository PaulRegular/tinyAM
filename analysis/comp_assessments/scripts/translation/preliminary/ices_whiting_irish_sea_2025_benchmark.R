source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1980:2023
ages <- 0:6
assessment_id <- source$assessment$assessment_id[[1]]
surveys <- c("NIGFS-Q1", "NIGFS-Q4", "NIMIK")

obs <- database_to_tam_obs(
  assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("NIGFS-Q1" = 0.23, "NIGFS-Q4" = 0.855, "NIMIK" = 0.48),
  assumptions = source$assumptions
)

obs$index$q_age <- ifelse(
  obs$index$survey == "NIGFS-Q1",
  ifelse(obs$index$age == 1, "age1", ifelse(obs$index$age <= 5, "ages2_5", "age6")),
  ifelse(
    obs$index$survey == "NIGFS-Q4",
    ifelse(obs$index$age == 0, "age0", ifelse(obs$index$age <= 3, "ages1_3", paste0("age", obs$index$age))),
    "age0"
  )
)
obs$index$q_key <- interaction(obs$index$survey, obs$index$q_age,
  drop = TRUE, lex.order = TRUE
)
obs$index$relative_sd <- 1
has_weight <- is.finite(obs$index$relative_precision_weight)
obs$index$relative_sd[has_weight] <-
  1 / sqrt(obs$index$relative_precision_weight[has_weight])

## Comparison ----
age_plus_group <- 6L
comparison_outputs <- source$outputs
comparison_aggregates <- c("abundance", "biomass", "ssb")

background <- c(
  "### Irish Sea whiting: 2025 ICES SAM benchmark",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | SAM fits 1980–2023 at ages 0–6+, with age 6 as the plus group for catch and two NIGFS indices. | Fit the same years and ages with age 6 as the comparison plus group. | The age-0-only NIMIK series is retained as observed. |",
  "| Recruitment and N | Recruitment is at age 0; SAM separates age-0 process variance from the shared ages 1–6 variance. | Use a random-walk recruitment process and an IID N process with exponential initialization. | The N transition covariance and initial-state treatment are simplified. |",
  "| F and M | SAM estimates age-specific F with AR(1) correlation across ages and reports Fbar for ages 1–3; M is fixed at ages 0–6. | Use correlated random walks across F ages and retain the source M surface. | tinyAM does not reproduce the full SAM correlation and variance structure. |",
  "| Catch | Total commercial catch-at-age is fitted with one shared lognormal SD. | Keep catch-at-age and use one shared catch SD. | No extra catch-age SD groups are introduced. |",
  "| Surveys and q | NIGFS-Q1, NIGFS-Q4 and NIMIK are fitted at timing 0.23, 0.855 and 0.48; q is shared by source age groups. NIGFS residuals are AR-correlated across ages. | Keep all three surveys, timing, q groups and NIGFS relative precision weights. | tinyAM uses independent index residuals; the two missing NIMIK observations (2020 and 2022) remain missing. |",
  "| Biology and SSB | SAM supplies annual stock/catch/landings/discard weights, maturity, fixed M and spawning-time fractions. | Use the recorded stock weight, maturity and fixed M surfaces. | Source SSB is retained with its native confidence interval; tinyAM uses its own fitted population and shared biology. |",
  "",
  "The fit tries tinyAM's normal initialization first. If numerical checks fail, it makes one retry with source-aligned starting values; these values initialize the optimizer and do not constrain estimates. The 2025 SAG summary series is stored separately and is not used as the benchmark reference."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) {
    if (is.null(start_par)) {
      tinyAM::fit_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = list(process = "iid", init = "exp"),
        rec_settings = list(process = "rw", mu_form = NULL),
        F_settings = list(process = "cor_rw", mean_ages = 1:3),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(sd_form = ~1),
        index_settings = list(
          q_form = ~0 + q_key,
          sd_form = ~0 + survey,
          sd_supplied = ~relative_sd
        ),
        silent = silent
      )
    } else {
      tinyAM::fit_tam(
        data = obs,
        years = years,
        ages = ages,
        N_settings = list(process = "iid", init = "exp"),
        rec_settings = list(process = "rw", mu_form = NULL),
        F_settings = list(process = "cor_rw", mean_ages = 1:3),
        M_settings = list(
          process = "off", mu_form = NULL,
          mu_supplied = ~M_assumption
        ),
        catch_settings = list(sd_form = ~1),
        index_settings = list(
          q_form = ~0 + q_key,
          sd_form = ~0 + survey,
          sd_supplied = ~relative_sd
        ),
        silent = silent,
        start_par = start_par
      )
    }
  }

  retry_start <- function() {
    source_N <- .translation_start_surface(
      source$outputs, "population", "numbers_at_age", years, ages
    )
    source_F <- .translation_start_surface(
      source$outputs, "mortality", "fishing_mortality_at_age", years, ages
    )
    source_M_rows <- .translation_surface(
      source$inputs[source$inputs$type == "M" &
        source$inputs$measure == "natural_mortality_at_age", , drop = FALSE],
      years, ages, "Official fixed M"
    )
    source_M <- matrix(
      source_M_rows$obs, nrow = length(years), ncol = length(ages),
      dimnames = list(as.character(years), as.character(ages))
    )
    names(dimnames(source_N)) <- names(dimnames(source_F)) <- c("year", "age")
    q_design <- stats::model.matrix(~0 + q_key, data = obs$index)
    survey_year <- match(as.character(obs$index$year), rownames(source_N))
    survey_age <- match(as.character(obs$index$age), colnames(source_N))
    survey_rows <- cbind(survey_year, survey_age)
    n_at_survey <- source_N[survey_rows] * exp(
      -(source_F[survey_rows] + source_M[survey_rows]) * obs$index$samp_time
    )
    q_start <- vapply(seq_len(ncol(q_design)), function(j) {
      selected <- q_design[, j] == 1
      stats::median(log(obs$index$obs[selected] / n_at_survey[selected]))
    }, numeric(1))
    names(q_start) <- colnames(q_design)

    list(
      log_r0 = log(source_N[1, "0"]),
      log_n0 = stats::setNames(log(source_N[1, as.character(ages[-1])]), as.character(ages[-1])),
      log_r = stats::setNames(log(source_N[-1, "0"]), as.character(years[-1])),
      log_n = log(source_N[-1, as.character(ages[-1]), drop = FALSE]),
      log_f = log(source_F),
      log_q = q_start
    )
  }

  fit <- fit_default_first(
    fit_model,
    retry_start = retry_start,
    attempt_file = file.path(
      "analysis", "comp_assessments", "results", "cache",
      assessment_id, "initialization_attempts.rds"
    )
  )
}
