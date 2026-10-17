source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1978:2018
ages <- 1:15
obs <- database_to_tam_obs(source$assessment$assessment_id, source$inputs,
  years = years, ages = ages, maturity_multiplier = .5,
  weight_survey = "AI_bottom_trawl", sampling_times = c(AI_bottom_trawl = .5),
  assumptions = source$assumptions)
obs$maturity$spawn_time <- 2 / 12
obs$index$q_block <- factor(pmin(obs$index$age, 8))

## Background and comparisons ----
age_plus_group <- 15
comparison_aggregates <- "ssb"
comparison_definitions <- list(ssb = list(status = "approximate",
  definition = "Native female spawning biomass versus tinyAM female biomass using annual survey weights",
  reason = "The separate official constant population-weight vector is unavailable; source SSB is retained, not reconstructed under invented biology."))
background <- c(
  "### Aleutian Islands pollock: official 2024 Model 15.1", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1978–2024, ages 1–15+. | 1978–2018, same ages. | Last matching fishery/survey age compositions are from 2018. Later aggregate-only observations cannot enter the age-index approximation without invented compositions. Official historical estimates also use those later data. |",
  "| N and recruitment | Deterministic survival; equilibrium initialization; Beverton–Holt recruitment with fixed steepness 0.7 and source variability setting 0.6. | Deterministic survival, free initial ages and a Beverton–Holt mean with estimated curve and IID log SD. | Source equilibrium, fixed steepness, variability scale and bias correction are not reproduced. The mean curve in tinyAM describes median recruitment. |",
  "| F | Annual intensity and penalized age selectivity, changing in four periods; ages 8+ share selectivity. | Independent temporal F random walks by age. | Source temporal blocks, smoothness penalties and terminal-age sharing are omitted. |",
  "| M | One estimated constant under a prior with mean 0.2 and CV 0.2. | One estimated constant without that prior. | M and q may be weakly identified; the official M estimate is not treated as a supplied input. |",
  "| Catch | Best-estimate total biomass and separate multinomial number compositions. | Numbers reconstructed from the same totals, proportions and annual catch weights, with estimated lognormal SD. | Source likelihood, ageing error and composition weights are not reproduced. |",
  "| Index and q | NRA survey biomass plus design-based age proportions; scalar q fixed at one under ages 5–12 normalization; penalized age selectivity. | Reconstructed age indices using original survey weights; q estimated by age, shared at ages 8+; estimated common log SD; time 0.5. | Exact configured timing is unresolved; 0.5 is a summer approximation. The fixed-q normalization, supplied total-index CVs and separate composition likelihood are omitted. |",
  "| Weights | Separate constant population weights and annual fishery/survey weights. | Original annual survey weights serve as population weights; catch conversion retains original fishery weights. | The current constant population vector was not recovered. No earlier assessment vector is silently substituted. |",
  "| Maturity and SSB | Adopted GOA female maturity, combined-sex fraction 0.5 and March spawning. | Same maturity and fraction; spawning time 2/12. | Native SSB is compared approximately because population weights differ. Shared-biology total biomass is labelled separately. |",
  "",
  "Official N is rounded to one million fish; native recruitment has finer precision and its own SE. Historical F-at-age is unavailable. The reported arithmetic mean F is not tinyAM’s population-weighted Fbar, and the terminal selectivity curve is not copied backwards.",
  "Default initialization runs first. Failed numerical checks permit one retry from aligned official N and constant M; historical F starts remain at package defaults. No observations are adjusted using fitted outputs."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ bh(ssb)),
  F_settings = list(process = "rw"),
  M_settings = list(process = "off", mu_form = ~ 1, mu_supplied = NULL),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(q_form = ~ q_block, sd_form = ~ 1),
  silent = silent, ...
)
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs, "population", "numbers_at_age",
                                     years, ages, 1e9)
    list(log_r0 = log(n[1, 1]), log_r = log(n[-1, 1]),
         log_n0 = log(n[1, -1]), mu_m = log(.206))
  }, attempt_file = file.path("analysis/comp_assessments/results/cache",
    source$assessment$assessment_id, "initialization_attempts.rds"))
}
