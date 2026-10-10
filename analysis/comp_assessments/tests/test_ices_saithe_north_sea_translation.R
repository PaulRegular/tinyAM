source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_saithe_north_sea_2026"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
translated <- .test_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

inputs <- source_data$inputs
outputs <- source_data$outputs
m_expected <- c(.384, .335, .294, .259, .232, .212, .197, .177)
m_actual <- vapply(translated$ages, function(age) {
  unique(translated$obs$weight$M_assumption[
    translated$obs$weight$age == age
  ])
}, numeric(1))
accepted_f_9plus <- outputs[
  outputs$type == "mortality" &
    outputs$measure == "fishing_mortality_at_age" &
    !is.na(outputs$age_group) & outputs$age_group == "9+" &
    outputs$year %in% translated$years,
  c("year", "value")
]
comparison_f_10plus <- translated$comparison_outputs[
  translated$comparison_outputs$type == "mortality" &
    translated$comparison_outputs$measure == "fishing_mortality_at_age" &
    !is.na(translated$comparison_outputs$age_group) &
    translated$comparison_outputs$age_group == "10+" &
    translated$comparison_outputs$year %in% translated$years,
  c("year", "value")
]

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1967:2025),
  identical(translated$ages, 3:10),
  nrow(translated$obs$catch) == 59L * 8L,
  nrow(translated$obs$weight) == 59L * 8L,
  nrow(translated$obs$maturity) == 59L * 8L,
  nrow(translated$obs$index) == 34L * 6L,
  all(translated$obs$index$samp_time == .75),
  setequal(unique(translated$obs$index$survey), "NS-IBTS Q3-Q4"),
  sum(inputs$type == "index" &
        inputs$measure == "relative_biomass_index") == 26L,
  !any(translated$obs$index$survey == "Combined commercial trawl CPUE"),
  isTRUE(all.equal(m_actual, m_expected)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$F_settings$process, "ar1"),
  identical(dat$F_settings$mean_ages, 4:7),
  identical(dat$M_settings$process, "off"),
  identical(dim(translated$start_par$log_f), c(59L, 8L)),
  identical(dim(translated$start_par$log_f), dim(par$log_f)),
  all(is.finite(translated$start_par$log_f)),
  all(is.finite(translated$start_par$log_n)),
  nrow(accepted_f_9plus) == 59L,
  nrow(comparison_f_10plus) == 59L,
  identical(accepted_f_9plus$year, comparison_f_10plus$year),
  isTRUE(all.equal(accepted_f_9plus$value, comparison_f_10plus$value)),
  !anyNA(translated$comparison_outputs$year),
  grepl("omitted from the fit", paste(translated$background, collapse = " "))
)

cat("North Sea saithe translation structure passed.\n")
