## Observations ----
years <- 1972:2026
ages <- 0:8
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)

# The SAM configuration assigns one observation SD per key, in age order.
sd_key <- rbind(
  "Residual catch" = c(0, 1, 2, 2, 2, 2, 2, 2, 2),
  "delta-GAMNS-WCQ1" = c(NA, 3, 4, 4, 5, 5, 5, 6, 6),
  "delta-GAMNS-WCQ3+Q4" = c(7, 8, 9, 10, 10, 10, 11, 11, 11)
)
obs$catch$sd_key <- factor(sd_key["Residual catch", obs$catch$age + 1L])
survey_row <- match(obs$index$survey, rownames(sd_key))
obs$index$sd_key <- factor(sd_key[cbind(survey_row, obs$index$age + 1L)])
power_key <- c(
  "delta-GAMNS-WCQ1.1" = "Q1 age 1",
  "delta-GAMNS-WCQ3+Q4.0" = "Q3+Q4 age 0",
  "delta-GAMNS-WCQ3+Q4.1" = "Q3+Q4 age 1"
)
obs$index$q_power_group <- unname(power_key[as.character(obs$index$q_key)])
obs$pop_covariates$N_sd_key <- factor(ifelse(obs$pop_covariates$age == 8, "age8plus", "ages1to7"))

# Supply catch weights separately. There is no catch in the 2026 fit year.
catch_weights <- subset(source$inputs, type == "catch_weight" & measure == "weight_at_age")
stopifnot(nrow(catch_weights) == length(years[-length(years)]) * length(ages),
  all(catch_weights$unit == "kg"), !anyDuplicated(paste(catch_weights$year, catch_weights$age)))
obs$weight$type <- "stock"
catch_weight_grid <- transform(obs$weight, type = "catch")
catch_weight_grid$obs <- catch_weights$value[match(
  paste(pmin(catch_weight_grid$year, 2025L), catch_weight_grid$age),
  paste(catch_weights$year, catch_weights$age)
)]
obs$weight <- rbind(obs$weight, catch_weight_grid)

## Background and comparisons ----
age_plus_group <- 8
background <- c(
  "### Northern Shelf haddock: official 2026 SAM assessment",
  "",
  print_sources(source$assessment),
  "",
  "The fit uses ages 0–8+ and years 1972–2026. The catch series ends in 2025. Both surveys retain their reported timing, native units and relative precision.",
  "",
  "Recruitment and fishing mortality follow random walks. Annual log-F increments are correlated across ages. Both first-year recruitment and older first-year abundance are integrated without an initial density. Older-cohort survival deviations share one SD at ages 1–7 and have a separate SD at age 8+.",
  "",
  "Catch and survey observations are independent on the log scale with the official age-specific SD groups. Baseline catchability varies by survey and age. The abundance power is estimated separately for Q1 age 1 and Q3+Q4 ages 0 and 1.",
  "",
  "Natural mortality, stock weights, maturity and beginning-of-year spawning use the official inputs. The 2025 catch weights fill the reporting grid for the survey-only year 2026. There is no observed yield for 2026.",
  "",
  "Fbar comparisons use the official age range (2–4), with both models averaged using their own abundance at age. The published SAM Fbar is an arithmetic mean. [Native model checks](https://github.com/PaulRegular/tinyAM/tree/main/analysis/sam_equivalence#northern-shelf-haddock-2026) describe the likelihood and reporting comparisons."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "latent", sd_form = ~ 0 + N_sd_key),
    rec_settings = list(process = "rw", init = "latent"),
    F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:4),
    M_settings = list(process = "off", mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 0 + sd_key),
    index_settings = list(
      q_form = ~ 0 + q_key + q_power(N, by = q_power_group),
      sd_form = ~ 0 + sd_key, sd_supplied = ~ relative_sd
    ),
    silent = silent
  )
}
