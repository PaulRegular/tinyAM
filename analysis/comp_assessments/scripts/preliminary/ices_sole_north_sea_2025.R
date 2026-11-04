## Observations ----
years <- 1957:2024
ages <- 1:10
assessment_id <- source$assessment$assessment_id[[1]]
inputs <- as.data.frame(source$inputs)

source_catch <- inputs[
  inputs$assessment_id == assessment_id &
    inputs$type == "catch" & inputs$measure == "numbers_at_age",
  , drop = FALSE
]
catch_sum <- stats::aggregate(
  value ~ year + age, source_catch, sum
)
catch_rows <- source_catch[match(
  paste(catch_sum$year, catch_sum$age),
  paste(source_catch$year, source_catch$age)
), , drop = FALSE]
catch_rows$value <- catch_sum$value
catch_rows$fleet <- ""
catch_rows$source_type <- "translation_assumption"
catch_rows$source_reference <- paste(
  unique(source_catch$source_reference), collapse = "; "
)
catch_rows$transformation <- paste(
  "Summed landings and discards by year and age where each source series is available.",
  "The discard series starts in 2002; earlier rows contain reported landings only."
)
catch_rows$notes <- "Analysis-only combined catch stream; canonical fleet-specific records are unchanged."
inputs <- rbind(
  inputs[!(inputs$assessment_id == assessment_id &
    inputs$type == "catch" & inputs$measure == "numbers_at_age"), ],
  catch_rows
)

surveys <- c("BTS Q3", "COAST (SNS + DYFS)")
inputs <- inputs[!(inputs$type == "index" &
  inputs$measure == "relative_biomass_index" & inputs$survey %in% surveys), ]
index_rows <- inputs$type == "index" & inputs$measure == "numbers_at_age" &
  inputs$survey %in% surveys
inputs$unit[index_rows] <- "native survey index (unit unresolved)"
obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = surveys,
  sampling_times = c("BTS Q3" = 0.625, "COAST (SNS + DYFS)" = 0.625),
  assumptions = source$assumptions
)
obs$index$q_block <- factor(
  ifelse(obs$index$age >= 5, "5+", as.character(obs$index$age)),
  levels = c(as.character(1:4), "5+")
)
obs$index$q_key <- interaction(
  obs$index$survey, obs$index$q_block,
  drop = TRUE, lex.order = TRUE
)
obs$index$sd_block <- factor(obs$index$survey)
obs$maturity$spawn_time <- 0

## Fit and comparison ----
age_plus_group <- 10L
comparison_outputs <- source$outputs
comparison_definitions <- list(
  recruitment = list(
    status = "matched",
    definition = "Official same-year N at age 1; native recruitment at age 0 is retained separately.",
    reason = "The tinyAM population starts at age 1 because catch and F are fully available from that age."
  ),
  ssb = list(
    status = "matched",
    definition = "Start-of-year mature biomass over modeled ages 1-10 using official N and the shared source weights and maturity.",
    reason = "The common-definition comparison excludes age 0; native aggregate SSB remains available for context."
  )
)

background <- c(
  "### North Sea sole: 2025 ICES Stock Synthesis assessment",
  "",
  print_sources(source$assessment),
  "",
  "| Component | Official assessment | tinyAM representation | Main difference |",
  "|---|---|---|---|",
  "| Years and ages | The model spans 1957-2024 and ages 0-10+, with recruitment at age 0. | Fit 1957-2024 at ages 1-10+, retaining the age-10 plus group. | Age 0 is omitted because no catch or F-at-age is reported there; recruitment is compared using same-year N at age 1. |",
  "| Recruitment and N | Beverton-Holt recruitment with fixed steepness 0.6646 and the Stock Synthesis population structure. | Use a random walk for log recruitment, an IID cohort-survival process, and exponential initial abundance. | tinyAM cannot hold steepness fixed; estimating both Beverton-Holt curve parameters produced a non-positive-definite Hessian, so the fit uses the simpler recruitment process instead. |",
  "| F | Landings and discards have separate selectivity; landings selectivity varies through time, while both are constant from age 5 onward. Fbar covers ages 2-6. | Use one random-walk F process across ages, with Fbar ages 2-6. | Fleet-specific selectivity and the source's age-sharing structure are not represented. |",
  "| M | Fixed age-specific Lorenzen mortality. | Use the same fixed age-specific M values. | M is an input in both representations, not an estimated tinyAM process. |",
  "| Catch | Separate landings and discard catch-at-age streams; discards begin in 2002. | Sum the available catch streams by age, using reported landings alone before the discard series begins. Fit one shared catch observation SD. | The fleets and their separate likelihoods are combined; no missing discard values are added to the database. |",
  "| Surveys and q | BTS and COAST (SNS + DYFS) age-specific Q3 indices, with q constant from age 5 onward; aggregate biomass indices are also used. | Fit the age-specific series, with a separate q for ages 1-4 and a shared q for ages 5-10 in each survey. | Aggregate indices are omitted to avoid fitting closely related aggregate and age-specific series together. The source index multiplier is not reported, so its scale is retained. |",
  "| Biology and SSB | Annual stock weights, fixed maturity, and start-of-year SSB. | Retain source weights, maturity and start-of-year timing. | The common-definition SSB comparison uses ages 1-10; source native SSB includes age 0. |",
  "",
  "The 2026 ICES record is summary-only. The 2025 report is the latest detailed assessment used here; later summaries are not mixed into the fitted period."
)

fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "rw", mu_form = NULL),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(
      process = "off", mu_form = NULL,
      mu_supplied = ~M_assumption
    ),
    catch_settings = list(sd_form = ~1),
    index_settings = list(
      q_form = ~0 + q_key,
      sd_form = ~0 + sd_block
    ),
    silent = silent
  )
}
