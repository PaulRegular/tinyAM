source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 1960:2024
ages <- 3:40
inputs <- source$inputs
# Keep the active biomass likelihood; omit length samples and inactive totals.
inputs <- inputs[!inputs$measure %in% c("frequency_at_length", "maturity_sample_size",
                                      "mature_sample_count") &
                   !(inputs$type == "index" & inputs$measure == "total_numbers"), ]
beta <- source$outputs$value[match(c("mat_beta1", "mat_beta2"), source$outputs$measure)]
mat <- source$inputs[source$inputs$type == "weight", ]
mat$type <- "maturity"
mat$measure <- "maturity_at_age"
mat$basis <- mat$unit <- "proportion"
mat$value <- plogis(beta[1] + beta[2] * mat$age)
mat$value[mat$age == 46] <- mean(c(tail(mat$value, 1), plogis(beta[1] + beta[2] * 100)))
mat$notes <- "Translation: fixed plug-in maturity from official fitted coefficients; not a supplied input."
inputs <- rbind(inputs, mat)
cw <- inputs[inputs$type == "catch_weight", ]
inputs <- rbind(inputs[inputs$type != "catch_weight", ],
                .translation_repeat_biology(cw, years, NULL, "Catch weights"))
obs <- database_to_tam_obs(source$assessment$assessment_id, inputs,
  years = years, ages = ages, maturity_multiplier = 0.5,
  sampling_times = c(AI_survey = 0.5, EBS_survey = 0.5),
  assumptions = source$assumptions)
obs$maturity$spawn_time <- 0.25
obs$pop_covariates$rec_era <- factor(ifelse(obs$pop_covariates$year < 1990,
                                          "before_1990", "from_1990"))

## Background and comparisons ----
age_plus_group <- 40
background <- c(
  "### BSAI Pacific ocean perch: official 2024 Model 24", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1960–2024, ages 3–46+; observations and exported N use 40+. | Same years; ages 3–40+. | Older ages are combined at the observation boundary. Source weights above 40 are retained for plus-group biology. |",
  "| N | Deterministic survival, equilibrium initial abundance; two recruitment means split at 1990; fixed log SD 0.75; last three recruitments fixed to a bias-corrected mean. | Deterministic survival, free initial ages, IID log recruitment with estimated SD around the same two mean periods. | Initial equilibrium, fixed recruitment SD, terminal recruitment fixing and the source bias correction are not reproduced. |",
  "| F | Annual fishing intensity and penalized bicubic selectivity. | Independent temporal random walks at each age. | The spline penalties and common annual intensity are not reproduced. |",
  "| M | One estimated constant under a tight prior (mean 0.05, CV 0.05). | Estimate one constant without that prior, using package initialization. | M and q may be weakly separated without the source priors. Source-estimated M is not fixed as an input. |",
  "| Catch | Biomass plus separate age and length samples. | Reconstruct numbers from biomass, number proportions and original catch weights. | Length samples and the source composition likelihood are omitted; age-specific lognormal errors replace them. |",
  "| Surveys | AI and EBS slope biomass indices and age samples; annual-average survival; constant q with logistic age selectivity. | Reconstructed age indices, separate rising logistic q curves and survey SDs; time 0.5. | Mid-year survival approximates an annual average. Source q priors, observation weights and ageing error are omitted. |",
  "| Biology and SSB | Supplied constant weights in grams; fitted female maturity; April 1 female SSB. | Original weights converted to kg, plug-in fitted maturity with female fraction 0.5, spawning time 0.25. | Maturity is held at the official point estimate. Plus-group age composition and coefficient rounding can affect SSB. |",
  "",
  "Comparisons use common historical years and reported 40+ groups. Native SSB and biomass remain available separately from quantities calculated under shared biology. Derived sums and F have no invented uncertainty.",
  "Default package starting values are tried first. Only failed numerical checks permit one retry from aligned official N and F; that retry changes initialization only."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ rec_era),
  F_settings = list(process = "rw"),
  M_settings = list(process = "off", mu_form = ~ 1, mu_supplied = NULL),
  catch_settings = list(sd_form = ~ age + I(age^2)),
  index_settings = list(q_form = ~ 0 + survey + logistic(age, by = survey),
                        sd_form = ~ 0 + survey),
  silent = silent, ...
)
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    n <- .translation_start_surface(source$outputs, "population", "numbers_at_age",
                                     years, ages, 1e6)
    f <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age",
                                     years, ages)
    list(log_r0 = log(n[1, 1]), log_r = log(n[-1, 1]),
         log_n0 = log(n[1, -1]), log_f = log(f))
  })
}
